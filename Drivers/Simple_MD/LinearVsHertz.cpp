//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "DPMBase.h"
#include <iostream>

/////////////////////////////////////////////////////////////////////////////////////////////////////
///This does the force computation
////////////////////////////////////////////////////////////////////////////////////////////////////
class Hertz : public DPMBase{
	
	void broadPhase(int i){
		for (int j=0;j<i;j++){
			computeInternalForces(i,j);
		}
	}

	void computeExternalForces(int CI)
	{
		//Add on gravity
		Particles[CI].Force += getGravity() * Particles[CI].get_mass();
		//Finally walls
		if (!Particles[CI].is_fixed()) computeWalls(CI);
		
	}
	
	void computeInternalForces(int P1, int P2);
	void computeWalls(int PI);

};

///fdotn computation is changed
void Hertz::computeInternalForces(int P1, int P2)
{
	Vec3D normal;
	Vec3D force;
	double dist;

	///Tangential spring information is always store in the real particle with highest index
	///When a Periodic contact is encountered it is always enstd::coutered twice, but only applied when the real particle has the highest index of both real indices
	///When a Particle is removed the tangential spring information has to be moved

	//Cases (start from P1 and P2 and go to PI and PJ
	//1 Normal-Normal		->PI=max(P1,P2), PJ=min(P1,P2)
	//2 Periodic-Normal		->if(P2>Real(P1)) (PI=P2 PJ=real(P1)) otherwise do nothing
	//3 Normal-Periodic		->if(P1>Real(P2)) (PI=P1 PJ=real(P2)) otherwise do nothing
	//4 Periodic-Periodic	->do nothing
	
	//Just some statements to handle the 4 cases
	int PI,PJ;
	int P1Per=Particles[P1].get_periodicFromParticle();
	int P2Per=Particles[P2].get_periodicFromParticle();
	int PJreal;
	if(P1Per==-1) {
		if(P2Per==-1) {
			//N-N
			PI=max(P1,P2);
			PJ=min(P1,P2);
			PJreal=PJ;
		} else {
			//N-P
			if(P1>P2Per) {
				PI=P1;
				PJ=P2;
				PJreal=P2Per;
			} else {
				return;
			}
		}
	} else {
		if(P2Per==-1) {
			//P-N
			if(P2>P1Per) {
				PI=P2;
				PJ=P1;
				PJreal=P1Per;
			} else {
				return;
			}
		} else {
			//P-P
			return;
		}
	}

	//this is because the rough bottom allows overlapping fixed particles
	if (Particles[PI].is_fixed()&&Particles[PJreal].is_fixed()) return;


	static Vec3D vrel, vrelt, forcet = Vec3D(0.0, 0.0, 0.0);
	static Vec3D deltat = Vec3D(0.0, 0.0, 0.0);
	
	#ifdef DEBUG_OUTPUT
		std::cerr << "In computing interal forces between particle "<<Particles[PI].Position<<" and "<<Particles[PJ].Position<<std::endl;
	#endif
	
	//Get the square of the distance between particle i and particle j
	double dist_squared=GetDistance2(Particles[PI].Position, Particles[PJ].Position);
	double radii_sum=Particles[PI].Radius+Particles[PJ].Radius;
	
	#ifdef DEBUG_OUTPUT_FULL
		std::cerr << "Square of distance between " << dist_squared << " square sum of radii " << radii_sum*radii_sum <<std::endl;
	#endif
	
	// True if the particles are in contact
	if (dist_squared<(radii_sum*radii_sum))
	{
		// For particles of the same species, set species vector to Species(PI);
		// for particles of different species, set species vector to MixedSpecies(PI,PJ)
		BaseSpecies* pSpecies = (Particles[PI].indSpecies==Particles[PJ].indSpecies)?&Species[Particles[PI].indSpecies]:speciesHandler.getMixedObject(Particles[PI].indSpecies,Particles[PJ].indSpecies);
		
		// Calculate distance between the particles
		dist=sqrt(dist_squared);
		
		// Compute normal vector
		normal=(Particles[PI].Position-Particles[PJ].Position)/dist;
		
		// Compute the overlap between the particles
		double deltan = radii_sum-dist;
		
		// Compute the relative velocity vector v_ij
		if (!pSpecies->mu) {
			vrel=(Particles[PI].Velocity-Particles[PJ].Velocity);
		} else {
			vrel=(Particles[PI].Velocity-Particles[PJ].Velocity) + Vec3D::cross(normal, Particles[PI].AngularVelocity * (Particles[PI].Radius - .5 * deltan) + Particles[PJ].AngularVelocity * (Particles[PJ].Radius - .5 * deltan) );
		}
		
		// Compute the projection of vrel onto the normal (can be negative)
		double vdotn=-Vec3D::dot(vrel,normal);
		
		// Compute normal force on particle i due to contact
		double fdotn = pSpecies->k*sqrt(deltan)*(deltan+pSpecies->disp*vdotn);
		force = normal * fdotn;
		//If tangential forces are present
		if (pSpecies->mu) { 
			
			//Compute the tangential component of vrel
			vrelt=vrel+normal*vdotn;

			//Compute norm of vrelt
			double vdott=vrelt.getLength();

			//Compute norm of normal force
			double norm_fn = abs(fdotn);
			
			if (!pSpecies->kt) { //if no tangential spring
				//tangential forces are modelled by a damper of viscosity dispt (sticking), 
				//but the force is limited by Coulomb friction (sliding):
				//f_t = -dispt*vrelt, if dispt*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
				if (vdott*pSpecies->dispt <= pSpecies->mus*norm_fn) { //sticking++;
					forcet = -pSpecies->dispt * vrelt; 
				} else { //sliding++;
					//set force to Coulomb limit
					forcet = -(pSpecies->mu * norm_fn / vdott) * vrelt;
				}
				
			} else { //if with tangential spring	
				//Retrieve the spring				
				Vec3D* delta = Particles[PI].TangentialSprings.select_particle(PJreal, getTime(), getTimeStep());
				
				//Integrate the spring
				//(*delta) += vrelt * dt;
				(*delta) += (vrelt - Vec3D::dot(*delta,Particles[PI].Velocity-Particles[PJ].Velocity)*normal/dist) *  getTimeStep();
				//std::cout << "check if tangential:" << Vec3D::Dot(normal,(*delta)/(delta->GetLength())) << std::endl;
				
				///\todo{The stored tangential spring is not made tangential, shouldn't this be the case?}
				//Retrieve tangential spring by subtracting normal component from total spring 
				deltat = (*delta);
				//To force the tangential spring into tangential direction use
				//deltat = (*delta) - normal*Vec3D::Dot(normal,(*delta));
				
				//Calculate test force including viscous force
				forcet = (-pSpecies->dispt) * vrelt - pSpecies->kt * deltat;
				double forcet2 = forcet.getLengthSquared();

				//tangential forces are modelled by a spring-damper of elastisity kt and viscosity dispt (sticking), 
				//but the force is limited by Coulomb friction (sliding):
				//f_t = -dispt*vrelt, if dispt*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
				if( forcet2 <= mathsFunc::square(pSpecies->mus*norm_fn) ) { 
					//sticking++;
				} else { 
					//sliding++;
					double norm_forcet = sqrt(forcet2);
					forcet *= pSpecies->mu * norm_fn / norm_forcet;
					///\todo{The spring should be cut back such that fdott=mu*fdotn. This is simple for dispt=0; we have to think about what happens in the sliding case with tang. dissipation; same for walls}
					//(*delta) = -(forcet + pSpecies->dispt * vrelt)/pSpecies->kt;
					(*delta) = forcet/(-pSpecies->kt);
				}
			} //end if tangential spring
			//Add tangential force to total force
			force += forcet;
			
			// Add torque due to tangential forces: t = Vec3D::Cross(l,f), l=dist*Wall.normal
			//~ add_torques(PI,PJreal, deltan);
			Vec3D Vec3D::cross = Vec3D::cross(normal, force);
			Particles[PI    ].Torque -= Vec3D::cross * (Particles[PI].Radius - .5 * deltan);
			Particles[PJreal].Torque -= Vec3D::cross * (Particles[PJ].Radius - .5 * deltan);

		} //end if tangential forces
		
		//~ add_forces(PI,PJreal);
		Particles[PI    ].Force+=force;
		Particles[PJreal].Force-=force;

		// output for ene and stat files:
		if (getEneFile().getSaveCurrentTimestep()) {
			addElasticEnergy(0.5 * (pSpecies->k  * sqrt(deltan) * mathsFunc::square(deltan) + pSpecies->kt * deltat.getLengthSquared()));
		}
		if (getFStatFile().getSaveCurrentTimestep()||getStatFile().getSaveCurrentTimestep()||getDoCGAlways()) {			
			double fdott = forcet.getLength();
			double deltat_norm = -deltat.getLength();
			
			///\todo{Define the center this way or are radii involved? Or maybe just use middle of overlap region?}
			Vec3D centre = 0.5 * (Particles[PI].Position + Particles[PJ].Position);
			
			///The second particle (i.e. the particle the force acts on) 
			///is always a flow particle 
			///\todo{Is it the first particle the force acts on?}
			
			if (!Particles[PI].is_fixed())
			{
				 gatherContactStatistics(PJreal,PI, centre, deltan, deltat_norm,fdotn,fdott,-normal,-(fdott?forcet/fdott:forcet));
				 fstat_file 
					<< getTime() << " "
					<< PJreal << " "
					<< PI << " "				
					<< centre << " "
					<< deltan << " "
					<< deltat_norm << " "
					<< fdotn << " "
					<< fdott << " "
					<< -normal << " "
					<< -(fdott?forcet/fdott:forcet) << std::endl;
			}
			if (!Particles[PJreal].is_fixed())
			{
				 gatherContactStatistics(PI,PJreal, centre, deltan, deltat_norm,fdotn,fdott,normal,(fdott?forcet/fdott:forcet));
				 fstat_file 
					<< getTime() << " "
					<< PI << " "
					<< PJreal << " "
					<< centre << " "
					<< deltan << " "
					<< deltat_norm << " "
					<< fdotn << " "
					<< fdott << " "
					<< normal << " "
					<< (fdott?forcet/fdott:forcet) << std::endl;
			}
		}
		
	} // end if particle i and j are overlapping
}

///fdotn computation is changed
void Hertz::computeWalls(int PI) 
{
	Vec3D normal;
	Vec3D force;
	double dist;

	//No need to compute interactions between periodic particle images and walls
	if(Particles[PI].get_periodicFromParticle()!=-1)
		return;
	
	static Vec3D vrel, vrelt, forcet = Vec3D(0.0, 0.0, 0.0);
	static Vec3D deltat = Vec3D(0.0, 0.0, 0.0);
		
	for (int i=0; i<get_NWall(); i++) 
	{
		
		// note: normal points away from particle i
		bool touch = Walls[i].getDistance_and_normal(Particles[PI], dist, normal);
		
		//If the wall is being touched (I think :Ant)
		if (touch) 
		{
			
			BaseSpecies* pSpecies = (Particles[PI].indSpecies==Walls[i].indSpecies)?&Species[Particles[PI].indSpecies]:speciesHandler.getMixedObject(Particles[PI].indSpecies,Walls[i].indSpecies);
			
			double deltan = Particles[PI].Radius-dist;

			// Compute the relative velocity vector v_ij
			if (!pSpecies->mu) {
				vrel = Particles[PI].Velocity;
			} else {
				vrel = Particles[PI].Velocity - Vec3D::cross(normal, Particles[PI].AngularVelocity) * dist;
			}
			
			// Compute the projection of vrel onto the normal (can be negative)
			double vdotn = -Vec3D::dot(vrel, normal);
			
			// Compute normal force on particle i due to contact
			double fdotn = pSpecies->k*sqrt(deltan)*(deltan+pSpecies->disp*vdotn);
			//double fdotn = pSpecies->k * deltan - pSpecies->disp * vdotn;
			force = normal * (-fdotn);

			if (pSpecies->mu) { // Tangential force limited by Coulomb force fdott <= pSpecies->mu * fdotn
				
 				// Compute the tangential component of vrel
				vrelt = vrel + normal * vdotn;
				
				double vdott=sqrt(vrelt.getLengthSquared());
				fdotn = abs(fdotn);
				
				if (!pSpecies->kt) { //if no tangential spring
					if (vdott*pSpecies->dispt <= pSpecies->mus*fdotn) { 
						//sticking, viscous force
						forcet = -pSpecies->dispt * vrelt; 
					}	else { 
						//sliding, set force to Coulomb limit
						forcet = -(pSpecies->mu * fdotn / vdott) * vrelt;
					}
					
				} else { //if with tangential spring
					
					//Retrieve the tangential spring				
					Vec3D* delta = Particles[PI].TangentialSprings.select_wall(i, getTime(), getTimeStep());

					//Integrate the spring
					(*delta) += vrelt * getTimeStep(); //no correction because the normal direction is constant

					//Retrieve tangengial spring by subtracting normal component from total spring 
					deltat = (*delta);
					//To force the tangential spring into tangential direction use
					//deltat = (*delta) - normal*Vec3D::Dot(normal,(*delta));
					
					//Calculate test force including viscous force
					forcet = (-pSpecies->dispt) * vrelt - pSpecies->kt * deltat;
					double forcet2 = forcet.getLengthSquared();
					
					//tangential forces are modelled by a spring-damper of elastisity kt and viscosity dispt (sticking), 
					//but the force is limited by Coulomb friction (sliding):
					//f_t = -dispt*vrelt, if dispt*vrelt<=mu_s*fdotn, f_t=mu+s*fdotn*t, else
					if( forcet2 <= mathsFunc::square(pSpecies->mus*fdotn) ) {
						//sticking++;
					} else { 
						//sliding++;
						double norm_forcet = sqrt(forcet2);
						forcet *= pSpecies->mu * fdotn / norm_forcet;
						//(*delta) = -(forcet + pSpecies->dispt * vrelt)/pSpecies->kt;
						(*delta) = forcet/(-pSpecies->kt);
					}
				} //end if (!pSpecies->kt)
				
				///\todo{How do you calculate Hertzian forces for walls and particles of different diameter?}
				//if (Hertzian) force *= sqrt(deltan/(2.0*Particles[PI].Radius));
				
				// Add torque due to tangential forces t = Vec3D::Cross(l,f), l=dist*Wall.normal
				force += forcet;
				Particles[PI].Torque += Vec3D::cross(normal, force) * dist;
				
			} //end if tangential forces
			
			// Add force due to contact
			Particles[PI].Force += force;
			
			if (getEneFile().getSaveCurrentTimestep()) {
				addElasticEnergy(0.5 * (pSpecies->k * sqrt(deltan) * mathsFunc::square(deltan) + pSpecies->kt * deltat.getLengthSquared()));
			}
			if (getFStatFile().getSaveCurrentTimestep()||getStatFile().getSaveCurrentTimestep()||getDoCGAlways()) {			
				double fdott = forcet.getLength();
				double deltat_norm = deltat.getLengthSquared();

				
				if (!Particles[PI].is_fixed())
				{
					 gatherContactStatistics(PI,-(i+1), Particles[PI].Position + normal*(Particles[PI].Radius-deltan), deltan, deltat_norm,fdotn,fdott,normal,-(fdott?forcet/fdott:forcet));
					 fstat_file 
						<< getTime() << " "
						<< PI << " "
						<< -(i+1) << " "
						<< Particles[PI].Position + normal*(Particles[PI].Radius-deltan) << " "
						<< deltan << " "
						<< deltat_norm << " "
						<< fdotn << " "
						<< fdott << " "
						<< normal << " "
						<< -(fdott?forcet/fdott:forcet) << std::endl;
				}
			} // end if touch
			
		}
	}
}					 

class collinearH : public Hertz{
public:
	void setupInitialConditions()
	{
		setSystemDimensions(3);
		setParticleDimensions(3);
		
		setXMax( 3);
		setYMax( 2);
		setZMax( 2);
		setGravity(Vec3D(0,0,0));
		setDensity(6/constants::pi);
		
		set_NWall(0);
		
		set_N(2);
		particleHandler.getObject(0)->Position=Vec3D( 1,1,1);
		particleHandler.getObject(0)->setVelocity(Vec3D( 1,0,0);)
		particleHandler.getObject(0)->getRadius()=0.5;
		particleHandler.getObject(1)->Position=Vec3D( 2,1,1);
		particleHandler.getObject(1)->setVelocity(Vec3D(-1,0,0));
		particleHandler.getObject(1)->getRadius()=0.5;
		
		double tc = 1e-1;
		setCollisionTimeAndRestitutionCoefficient(tc,1,1);
		setStiffness(getStiffness()*5); //equal to linear model at deltan=1/5;
				
		setTimeStep(tc/50.0);
		setTimeMax(tc*2.0);
		setSaveCount(1);
		
	}
	
};

class collinear : public DPMBase{
public:
	void setupInitialConditions()
	{
		setSystemDimensions(3);
		setParticleDimensions(3);
		
		setXMax( 3);
		setYMax( 2);
		setZMax( 2);
		setGravity(Vec3D(0,0,0));
		setDensity(6/constants::pi);
		
		set_NWall(0);
		
		set_N(2);
		particleHandler.getObject(0)->Position=Vec3D( 1,1,1);
		particleHandler.getObject(0)->setVelocity(Vec3D( 1,0,0));
		particleHandler.getObject(0)->getRadius()=0.5;
		particleHandler.getObject(1)->Position=Vec3D( 2,1,1);
		particleHandler.getObject(1)->setVelocity(Vec3D(-1,0,0));
		particleHandler.getObject(1)->getRadius()=0.5;
		
		double tc = 1e-1;
		setCollisionTimeAndRestitutionCoefficient(tc,1,1);		
		setTimeStep(tc/50.0);
		setTimeMax(tc*2.0);
		setSaveCount(1);
		
	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{	
	collinear md;
	md.set_name("collinear");
	md.solve();	
	md.write(std::cout,true);

	collinearH mdH;
	mdH.set_name("collinearH");
	mdH.solve();	
	mdH.write(std::cout,true);
}
