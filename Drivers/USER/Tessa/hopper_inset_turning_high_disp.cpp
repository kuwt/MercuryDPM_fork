//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "ChuteWithHopperAndInset.h"
#include "StatisticsVector.h"
#include "Boundaries/HopperInsertionBoundary.h"

using namespace std;

class ChuteWithHopperAndInset_time : public ChuteWithHopperAndInset{
public:
	
	int num_created_triangle;
	int max_failed_triangle;
	int create_in_hopper;
	int max_created_triangle;
	int disp_old;
	
	void set_disp_old(Mdouble new_) {disp_old=new_;}
	Mdouble getDissipation_old() {return disp_old;}
	
	void set_create_in_hopper(unsigned int new_) {create_in_hopper = new_;}
	int get_create_in_hopper() {return create_in_hopper;}
	
	void set_max_created_triangle(unsigned int new_) {max_created_triangle = new_;}
	int get_max_created_triangle() {return max_created_triangle;}
	
	void set_num_created_triangle(unsigned int new_) {num_created_triangle = new_;}
	int get_num_created_triangle() {return num_created_triangle;}
	
	void setMaxFailed_triangle(unsigned int new_) {max_failed_triangle = new_;}
	
	void setupInitialConditions()
    {
        setupSideWalls();
        
        createBottom();
        
        addHopper();
        
        add_Inset();
    }

	void create_inflow_particle_triangle()
	{
		//the following formula yields polydispersed particle radii:
		P0.setRadius(random.getRandomNumber(getMinInflowParticleRadius(),getMaxInflowParticleRadius()));

		double insetHeight = get_InsetHeight();
		double insetWidth = get_InsetWidth();
		double insetAngle = get_InsetAngle();
		double a = insetHeight-insetWidth*(sin(insetAngle)/cos(insetAngle));
		double c= 0.5*(getHopperLength()-getHopperExitLength()); 
		
		Vec3D g = getGravity();
		double theta = atan(g.Z/g.X);
		
		double l=a/cos(theta);
		double h=(getMinInflowParticleRadius()+getMaxInflowParticleRadius())*5;
		
		double x1=random.getRandomNumber(0.0,1.0)*0.5*l+l*0.25;
		double z1=random.getRandomNumber(0.0,1.0)*h;
		
		P0.setPosition( Vec3D(
		x1*sin(theta)+z1*cos(theta) + insetWidth + c , 
		getMaxInflowParticleRadius() , 
		-x1*cos(theta)+z1*sin(theta) + a
		));
		Vec3D position;
		position = Vec3D(0.0,0.0,getHopperLift()) + P0.getPosition();
		P0.setPosition(position);
	}
	
	//change add_particles so there is a phase where no particles are added 
	void add_particles() 
	{
		//cout<<"Starting function add_particles"<<endl;
		int failed = 0;
		int failed_triangle=0;
		
		Mdouble time= getTime();
		Mdouble tmax=getTimeMax();
		
		Mdouble t1 = tmax*0.25; // end filling time keep constant gravityangle and let the system relax
		Mdouble t2 = tmax*0.5; // end relaxation time, start rotating
		Mdouble t3 = tmax*0.75; // end rotating, start adding flowing particles
		Mdouble t4 = tmax; // end simulation time
		
		bool adding_particles;
		
		if (time < t1){ /// filling
			//during 0<t<t1: constant gravity under angle and add particles 
			//to put particles in initial position
			adding_particles=true;
			set_create_in_hopper(0);
		} else {
			if (time < t2) { /// relaxation
				//during t1<t<t2: relaxing no added particles
				adding_particles=false;
				set_create_in_hopper(0);
			} else {
				if (time < t3) { /// turning
					// during t2<t<t3: turning phase no added particles
					adding_particles=false;
					set_create_in_hopper(0);
				} else { ///flowing particles
					// during t3<t<t4: normal gravity and add particles
					adding_particles=true;
					set_create_in_hopper(1);
				}
			}
		}
		
		
		//try max_failed times to find new insertable particle (only if I want to add particles)
		if (adding_particles) {
			//cout<<"Adding particle"<<endl;
			if (get_create_in_hopper()) {
                //cout<<"Adding particles in hopper"<<endl;
                HopperInsertionBoundary b1;
                b1.set(new SphericalParticle, getMaxFailed(), getYMin(), getYMax(), getChuteAngle(),
                       getFixedParticleRadius(), getIsHopperCentred(), getHopperDimension(), getHopperAngle(),
                       getHopperLength(), getHopperExitLength(), getHopperHeight(), getHopperLift(),
                       getHopperFillingPercentage());
                boundaryHandler.copyAndAddObject(b1);
                setInsertionBoundary(dynamic_cast<InsertionBoundary*>(boundaryHandler.getLastObject()));

//				while (failed<=max_failed){
//					create_inflow_particle();
//					if (IsInsertable(P0)) {
//						failed = 0; 
//						num_created++;
//					} else failed++;
            };
        }
        else
        {
            //cout<<"Adding particles in triangle"<<endl;
            while (failed_triangle <= max_failed_triangle && num_created_triangle < max_created_triangle)
            {
                create_inflow_particle_triangle();
                if (checkParticleForInteraction(P0))
                {
                    failed_triangle = 0;
                    num_created_triangle++;
                }
                else failed_triangle++;
            };
        }
		}
	
	void actionsBeforeTimeStep(){
		
		
		// Add the option to change the gravity with time
		Mdouble time= getTime();
		Mdouble gravity = getGravity().getLength();;
		Mdouble GravityAngle;
		Mdouble tmax=getTimeMax();
		Mdouble dt = getTimeStep();
		
		Mdouble t1 = tmax*0.25; // end filling time keep constant gravityangle and let the system relax
		Mdouble angle1 = -30; // angle of gravity during filling
		Mdouble t2 = tmax*0.5; // end relaxation time, start rotating
		Mdouble angle2 = 0; // final angle 
		Mdouble t3 = tmax*0.75; // end rotating, start adding flowing particles
		Mdouble t4 = tmax; // end simulation time
		
		Mdouble omega = (angle2-angle1)/(t3-t2); // constant rotational velocity during t2<t<t3 [deg/s]
		
		
		
		if (time < t1){ /// filling
			//during 0<t<t1: constant gravity under angle and add particles 
			//to put particles in initial position
			GravityAngle=angle1*constants::pi/180;
			setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
		} else {
			if (time < t2) { /// relaxation
				//during t1<t<t2: relaxing no added particles
				GravityAngle=angle1*constants::pi/180;
				setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
			} else {
				if (time < t3) { /// turning
					// during t2<t<t3: turning phase no added particles
					GravityAngle=+dt*omega*constants::pi/180;
					setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
					species->setDissipation(45*disp_old);
					species->setSlidingDissipation(45*disp_old);
				} else { ///flowing particles
					// during t3<t<t4: normal gravity and add particles
					GravityAngle=angle2*constants::pi/180;
					setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
					species->setDissipation(disp_old);
					species->setSlidingDissipation(disp_old);
				}
			}
		}
		
		// These two lines were already part of the original function Chute::actionsBeforeTimeStep
		//if (get_create_in_hopper()) {
			add_particles();
			cleanChute();
		//}
	}
	private:
        	SphericalParticle P0;
};
int main(int argc UNUSED, char *argv[] UNUSED)
{
	ChuteWithHopperAndInset_time problem;
	///Initialisatie
	problem.set_num_created_triangle(0);
	
	
	///Problem parameters
	problem.setName("hopper_inset_turning_high_disp");
	problem.setTimeMax(10.0);
	
	///Particle properties
	problem.species->setDensity(2500.0); //uit data sheet M type beads
	problem.setInflowParticleRadius(5.7e-3,6.3e-3);
	//problem.setInflowParticleRadius(5.7e-3/2.0,6.3e-3/2.0);
    Mdouble r =problem.getInflowParticleRadius();
	Mdouble mass=4/3*constants::pi*r*r*r*problem.species->getDensity();
	problem.species->setCollisionTimeAndRestitutionCoefficient(4.0e-4, 0.8,mass);
	problem.set_disp_old(problem.species->getDissipation());
	// set a rough particle bottom
	problem.setFixedParticleRadius(problem.getInflowParticleRadius()); // if FixedParticleRadius < 1e-12, the function Chute::createBottom (called by Chute.cc) creates a smooth instead of rough bottom
	problem.setRoughBottomType(MONOLAYER_ORDERED);
	problem.species->setSlidingFrictionCoefficient(0.8);
	
	int max_created_triangle=300; //300 is een best aantal voor problem.setInflowParticleRadius(5.7e-3,6.3e-3);
	problem.set_max_created_triangle(max_created_triangle);
	
	//Mdouble eps = problem.getRestitutionCoefficient();
	//Mdouble r =problem.getInflowParticleRadius();
	//Mdouble mass=4/3*constants::pi*r*r*r*problem.getDensity();
	//Mdouble tc=problem.getCollisionTime();
	//Mdouble disp = - mass / tc * log(eps);
	//Mdouble k = .5 * mass * (mathsFunc::square(constants::pi/tc) + mathsFunc::square(disp/mass));
	//cout << "STIJFHEID: "<<k<<endl;
	
	///Chute properties
	problem.setChuteAngle(0.0);
	problem.setChuteLength(1.005);
	problem.setChuteWidth(7.5e-3);
	problem.setMaxFailed(1);
	problem.setMaxFailed_triangle(1);
	
	///Hopper & Inset properties
	double Height = 0.535+0.20, ExitHeight = 0.55, ExitLength = 60.0e-3, hopperAngle_ = 60.0, hopperLength_ = 8.0 * ExitLength;
	problem.setHopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_, Height);
	//set_Hopper(0.01, 0.01, 60.0, 0.08, 0.02);
	problem.setHopperLift(0.0);
	double insetHeight=0.332,insetWidth=0.21,insetAngle=30;
	problem.set_Inset(insetHeight,insetWidth,insetAngle);
	
	///Rolling resistance
	//problem.species->setRollingStiffness(2./5.*problem.species->getStiffness());
    //problem.species->setRollingFrictionCoefficient(10);
    //problem.species->setRollingDissipation(2./5.*problem.species->getDissipation());
	
	///time settings
	//speed allowed before particles move through each other!
	//cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; 
	//cout << "Maximum possible speed of particles: " << problem.getMaximumVelocityInducedByGravity() << endl; 
	problem.setTimeStep(4e-4/50.0);
	problem.DPMBase::setTimeStep(problem.getTimeStep()*5.0);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(100,problem.getTimeMax(),problem.getTimeStep()));
	//problem.setSaveCount(1);
	cout << "dt=" << problem.getTimeStep() << endl;
	problem.solve();
	
	problem.writeRestartFile();
}
