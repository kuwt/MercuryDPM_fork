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

#include "scr/Mercury3D.h"
 
class ScrewConveyor : public Mercury3D{
	private:
	
	void printTime() const {
		std::cout <<"t="<< getTime()<<", tmax=" << getTimeMax() << std::endl;
		std::cout.flush();
	}
	
	void setupInitialConditions()
	{
		std::cout<<"Setup_particles_initial_conditions()"<<std::endl;
		wallHandler.clear();
		InfiniteWall w0;
		w0.set(Vec3D( 0, 0,-1), -getZMin());
		wallHandler.copyAndAddObject(w0);
		
		CylindricalWall* cylindricalWall=new CylindricalWall(1.0);
		wallHandler.addObject(cylindricalWall);
		
		//Screw(Start,Length,Radius,Turn,Omega,Thickness)
		Screw* screw=new Screw(Vec3D(0,0,0),getZMax(),1.0,4,-1.0,0.5*ParticleRadius);
		wallHandler.addObject(screw);
		
		DeletionBoundary* deletionBoundary=new DeletionBoundary();
		deletionBoundary->set(Vec3D(0,0,1),getZMax());
		boundaryHandler.addObject(deletionBoundary);
		
		particleHandler.clear();
		BaseParticle p0;
		p0.setVelocity(Vec3D(0.0,0.0,0.0));
		p0.setRadius(BaseParticleRadius);
		p0.computeMass(Species);

		int Nx=20;
		int Ny=20;
		int Nz=20;

		Mdouble distance;
		Vec3D normal;		
		for(int i=0;i<Nx;i++)
			for(int j=0;j<Ny;j++)
				for(int k=0;k<Nz;k++)
				{
					p0.setPosition(Vec3D(getXMin()+(getXMax()-getXMin())*(0.5+i)/Nx,getYMin()+(getYMax()-getYMin())*(0.5+j)/Ny,getZMin()+0.4*(getZMax()-getZMin())*(0.5+k)/Nz));
					if(!dynamic_cast<CylindricalWall*>(wallHandler.getObject(1))->getDistance_and_normal(p0, distance, normal))
					if(!dynamic_cast<Screw*			 >(wallHandler.getObject(2))->getDistance_and_normal(p0, distance, normal))
						particleHandler.copyAndAddObject(p0);
				}
		/*p0.setPosition(Vec3D(0,-1.0+1.1*ParticleRadius,0.5));
		if(!dynamic_cast<CylindricalWall*>(wallHandler.getWall(1))->getDistance_and_normal(p0, distance, normal))
		{
			if(!dynamic_cast<Screw*			 >(wallHandler.getWall(2))->getDistance_and_normal(p0, distance, normal))
			{
				particleHandler.copyAndAddObject(p0);
			}
			else
			{
				std::cout<<"Particle hits screw"<<std::endl;
			}
		}
		else
		{
			std::cout<<"Particle hits cylinder"<<std::endl;
		}*/
	}
	
	void actionsBeforeTimeStep() {
	//	if(getTime()>2.0)
			wallHandler.getObject(2)->move_time(getTimeStep());			
	}
	
	public:
		double ParticleRadius;	

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	ScrewConveyor problem;
	
	
	//problem.set_name("Screw_test_ct_0.0001");
	problem.set_name("ScrewConveyor");
	problem.setDensity(1000);
	problem.setXMin(-1.0);
	problem.setXMax( 1.0);
	problem.setYMin(-1.0);
	problem.setYMax( 1.0);
	problem.setZMin( 0.0);
	problem.setZMax( 4.0);
	problem.setTimeMax( 5.0);
	//problem.ParticleRadius=0.05;	
	problem.ParticleRadius=0.05;	
	problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(0.001,0.4, pow(problem.ParticleRadius,3)*constants::pi*4.0/3.0*1000);
	problem.setGravity(Vec3D(0,-7,-7));
	//problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(0.0001,0.4, pow(problem.ParticleRadius,3)*constants::pi*4.0/3.0*1000);
	problem.setTimeStep();
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(1001,getTimeMax(),getTimeStep()));
	problem.set_HGRID_max_levels(1);
	problem.set_HGRID_sphere_to_cell_ratio(1.1);
	problem.set_UpdateEachTimeStep(false);
	problem.write(std::cout,false);
	problem.solve();
}
