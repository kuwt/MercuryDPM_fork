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
#include "Statistics.h"
#include "Boundaries/HopperInsertionBoundary.h"

using namespace std;

class ChuteWithHopperAndInset_time : public ChuteWithHopperAndInset{
public:
	
	int numCreatedHopper;
	int numCreatedTriangle;
	int maxFailedTriangle;
	int createInHopper;
	int maxCreatedHopper;
	int maxCreatedTriangle;
	Mdouble angle1;
	Mdouble angle2;
	Mdouble t1;
	Mdouble t2;
	Mdouble t3;
	Mdouble t4;
	
	void set_CreateInHopper(unsigned int new_) {createInHopper = new_;}
	int get_CreateInHopper() {return createInHopper;}
	
	void set_MaxCreatedHopper(unsigned int new_) {maxCreatedHopper = new_;}
	int get_MaxCreatedHopper() {return maxCreatedHopper;}
	
	void set_NumCreatedHopper(unsigned int new_) {numCreatedHopper = new_;}
	int get_NumCreatedHopper() {return numCreatedHopper;}
	
	void set_MaxCreatedTriangle(unsigned int new_) {maxCreatedTriangle = new_;}
	int get_MaxCreatedTriangle() {return maxCreatedTriangle;}
	
	void set_NumCreatedTriangle(unsigned int new_) {numCreatedTriangle = new_;}
	int get_NumCreatedTriangle() {return numCreatedTriangle;}
    
    int getNumCreated()
    {return numCreatedTriangle+numCreatedHopper;}
	
	void set_MaxFailedTriangle(unsigned int new_) {maxFailedTriangle = new_;}
	
	void set_Angle1(Mdouble new_) {angle1=new_;} //all input of angle is in degrees
	Mdouble get_Angle1() {return angle1;}
	
	void set_Angle2(Mdouble new_) {angle2=new_;} //all input of angle is in degrees
	Mdouble get_Angle2() {return angle2;}
	
	void set_T1(Mdouble new_) {t1=new_;}
	Mdouble get_T1() {return t1;}
	void set_T2(Mdouble new_) {t2=new_;}
	Mdouble get_T2() {return t2;}
	void set_T3(Mdouble new_) {t3=new_;}
	Mdouble get_T3() {return t3;}
	void set_T4(Mdouble new_) {t4=new_;}
	Mdouble get_T4() {return t4;}
    
    void setupInitialConditions() override
    {
        setupSideWalls();
        
        createBottom();
        
        addHopper();
        
        add_Inset();
    }
	
	
	void create_inflow_particle_triangle()
	{
		//the following formula yields polydispersed particle radii:
        inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));

		double insetHeight = get_InsetHeight();
		double insetWidth = get_InsetWidth();
		double insetAngle = get_InsetAngle(); //output is in radians
		//cout<<"cppfile:create_inflow_particle_triangle, insetAngle="<<insetAngle<<endl;
		double a = insetHeight-insetWidth*(sin(insetAngle)/cos(insetAngle)); //a=triangleHeight;
		double c= 0.5*(getHopperLength()-getHopperExitLength());
		
		//Vec3D g = getGravity();
		//double theta = atan(g.Z/g.X); //theta=90graden-insetAngle(grad) OR theta=1/2pi rad-insetAngle(rad)
		
		double theta=0.5*constants::pi-insetAngle;
		/*
		cout<<"---------------"<<endl;
		cout<<"theta "<<theta<<endl;
		cout<<"alpha "<<insetAngle<<endl;
		cout<<"1/2pi-alpha "<<0.5*constants::pi-insetAngle<<endl;
		*/
		
		double l=a/cos(theta);
		double h=sin(theta)*a; //a=triangleHeight;
		double b=cos(theta)*a;
		
		double x1,z1;
		do
		{
			//double h=(MinInflowParticleRadius+MaxInflowParticleRadius)*10;
			
			//double x1=random.get_RN(0.0,1.0)*0.5*l+l*0.25;
			x1=random.getRandomNumber(0.0,1.0)*l;

			//l*0.8+0.05*l: levert met hoek van 30 graden en 1400 deeltjes van 3mm radius wel goede driehoek op
			z1=random.getRandomNumber(0.0,1.0)*h;
		} while (z1>h/b*x1 || z1>h/(l-b)*(l-x1));
		
		//
//		if (z1>h/b*x1) {
			//deeltje ligt boven de lijn van (0,0) naar (b,h): z=h/b x
			//doe niks (oftewel dit deeltje wordt weggegooid)
	//	} else {
		//	if (z1>h/(l-b)*(l-x1)) {
				//deeltje ligt boven de lijn van (b,h) naar (l,0): z=h/(l-b)*(l-x)
				//doe ook niks (oftewel dit deeltje wordt weggegooid)
			//} else {
				//het deeltje valt binnen de gewenste driehoek en mag blijven
				inflowParticle_.setPosition( Vec3D(
				x1*sin(theta)+z1*cos(theta) + insetWidth + c , 
				getMaxInflowParticleRadius() , 
				-x1*cos(theta)+z1*sin(theta) + a
				));
                Vec3D position;
                position=inflowParticle_.getPosition();
                position+=Vec3D(0.0,0.0,getHopperLift());
				inflowParticle_.setPosition(position);
			//};
		//};
	}
	
	//change add_particles so there is a phase where no particles are added 
	void add_particles() 
	{
		//cout<<"Starting function add_particles"<<endl;
		int failed = 0;
		int failedTriangle=0;
		
		Mdouble time= getTime();
		Mdouble tmax=getTimeMax();
		
		Mdouble t1 = get_T1(); // end filling time keep constant gravityangle and let the system relax
		Mdouble t2 = get_T2(); // end relaxation time, start rotating
		Mdouble t3 = get_T3(); // end rotating, start adding flowing particles
		Mdouble t4 = get_T4(); // end simulation time
		
		bool addingParticles;
		
		if (time < t1){ /// filling
			//during 0<t<t1: constant gravity under angle and add particles 
			//to put particles in initial position
			addingParticles=true;
			set_CreateInHopper(0);
		} else {
			if (time < t2) { /// relaxation
				//during t1<t<t2: relaxing no added particles
				addingParticles=false;
				set_CreateInHopper(0);
			} else {
				if (time < t3) { /// turning
					// during t2<t<t3: turning phase no added particles
					addingParticles=false;
					set_CreateInHopper(0);
				} else { ///flowing particles
					// during t3<t<t4: normal gravity and add particles
					addingParticles=true;
					set_CreateInHopper(1);
				}
			}
		}
		
		
		//try max_failed times to find new insertable particle (only if I want to add particles)
		if (addingParticles) {
			//cout<<"Adding particle"<<endl;
			if (get_CreateInHopper()) {
                //cout<<"Adding particles in hopper"<<endl;
                HopperInsertionBoundary b1;
                b1.set(new SphericalParticle, getMaxFailed(), getYMin(), getYMax(),
                       getChuteAngle(), getFixedParticleRadius(), getIsHopperCentred(), getHopperDimension(),
                       getHopperAngle(), getHopperLength(),
                       getHopperExitLength(), getHopperHeight(), getHopperLift(), getHopperFillingPercentage());
                boundaryHandler.copyAndAddObject(b1);
                setInsertionBoundary(dynamic_cast<InsertionBoundary*>(boundaryHandler.getLastObject()));
                
                /*
                while (failed<=max_failed){
                    create_inflow_particle();
                    if (IsInsertable(P0)) {
                        failed = 0;
                        num_created++;
                    } else failed++;
                };
                */
            } else {
				//cout<<"Adding particles in triangle"<<endl;
				while (failedTriangle <= maxFailedTriangle && numCreatedTriangle < maxCreatedTriangle) {
					create_inflow_particle_triangle();
					if (checkParticleForInteraction(inflowParticle_)) {
						failedTriangle = 0;
						numCreatedTriangle++;
					} else failedTriangle++;
				};	
			}
		}
	}
	
	void actionsBeforeTimeStep() override {
		
		
		// Add the option to change the gravity with time
		Mdouble time= getTime();
		Mdouble gravity = getGravity().getLength();;
		Mdouble GravityAngle;
		Mdouble tmax=getTimeMax();
		Mdouble dt = getTimeStep();
		
		Mdouble t1 = get_T1(); // end filling time keep constant gravityangle and let the system relax
		Mdouble t2 = get_T2(); // end relaxation time, start rotating
		Mdouble t3 = get_T3(); // end rotating, start adding flowing particles
		Mdouble t4 = get_T4(); // end simulation time
		
		
		Mdouble angle1=get_Angle1();
		Mdouble angle2=get_Angle2();
		
		//Mdouble omega = (angle2-angle1)/(t3-t2); // constant rotational velocity during t2<t<t3 [deg/s]
		
		
		
		if (time < t1){ /// filling
			//during 0<t<t1: constant gravity under angle and add particles 
			//to put particles in initial position
			GravityAngle=angle1*constants::pi/180.0;
			setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
		} else {
			if (time < t2) { /// relaxation
				//during t1<t<t2: relaxing no added particles
				GravityAngle=angle1*constants::pi/180.0;
				setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
			} else {
				if (time < t3) { /// turning
					// during t2<t<t3: turning phase no added particles
					//cout<<"WAARDES STAAN HIER: "<<"angle1:"<<angle1<<", angle2:"<<angle2<<", time:"<<time<<", t2:"<<t2<<", t3:"<<t3<<endl;
					GravityAngle=((angle2-angle1)*(1-cos((time-t2)/(t3-t2)*constants::pi))*0.5+angle1)*constants::pi/180.0;
					setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
				} else { ///flowing particles
					// during t3<t<t4: normal gravity and add particles
					GravityAngle=angle2*constants::pi/180.0;
					setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
				}
			}
		}
		
		// These two lines were already part of the original function Chute::actionsBeforeTimeStep
		//if (get_CreateInHopper()) {
			add_particles();
			cleanChute();
		//}
	}
	
	void actionsBeforeTimeLoop() override {
		
		//automatically sets dt if dt is not specified by the user
        //if (!getTimeStep()) setTimeStep(); // change "dt" to "getTimeStep()", otherwise error: scr/MD.h:643: error: 'Mdouble MD::dt' is private
		
		
		Mdouble angle1 = get_Angle1();
		Mdouble angle2 = get_Angle2();
		Mdouble t1 = get_T1();
		Mdouble t2 = get_T2();
		Mdouble t3 = get_T3();
		Mdouble t4 = get_T4();
				
		ofstream myfile;
		myfile.open ("gravity_file.txt");
		myfile<<"Constants: "<<"angle1:"<<angle1<<", angle2:"<<angle2<<", t1:"<<t1<<", t2:"<<t2<<", t3:"<<t3<<", t4:"<<t4<<endl;
		myfile.close();
		
		
	};
	
	
	void actionsAfterTimeStep() override {
		//save the gravity information in "gravity_file.txt" after every 100th time step
		Mdouble factor = 10; //100
		
		Mdouble time= getTime();
		Mdouble dt= getTimeStep();
		
		if (abs(fmod(time+0.5*dt,factor*dt)-0.5*dt)<dt*1e-2) {
			//ugly fix to be able to compare two doubles
			//(als je door een afrondingsfout net onder factor*dt zit, komt de modulus ervan niet net onder nul te liggen, als je er eerst 0.5*dt bij optelt, dan modulus doet en daarna weer 0.5*dt van aftrekt wel)
			Vec3D gravity = getGravity();
			fstream myfile;
			myfile.open ("gravity_file.txt",fstream::in | fstream::out | fstream::app);
			//myfile << "links "<< abs(fmod(time,factor*dt))<<" rechts "<<dt*1e-2<<endl;
			myfile <<"time:"<<fixed<<setprecision(5)<<setw(9)<<time<<", g:"<<setw(17)<<setprecision(2)<<gravity<<endl;
			myfile.close();
		} else {
			//fstream myfile;
			//myfile.open ("gravity_file.txt",fstream::in | fstream::out | fstream::app);
			//myfile <<"time:"<<fixed<<setprecision(5)<<setw(9)<<time<<" this time is skipped"<<endl;
			//myfile.close();
		}
	
	}
    
private:
     SphericalParticle inflowParticle_;
};
int main(int argc UNUSED, char *argv[] UNUSED)
{
	ChuteWithHopperAndInset_time problem;
	///Initialisatie
	problem.set_NumCreatedTriangle(0);
	
	
	///Problem parameters
	problem.setName("hopper_inset_turning_cos_rolling_resistance5");
	Mdouble tmax=5.0;
	problem.setTimeMax(tmax);
	
	///Particle properties
	problem.species->setDensity(2500.0); //uit data sheet M type beads
	//problem.setInflowParticleRadius(5.7e-3,6.3e-3);
	problem.setInflowParticleRadius(5.7e-3/2.0,6.3e-3/2.0);
	
	// set a bottom made of regular placed particles
	problem.setFixedParticleRadius(problem.getInflowParticleRadius()); // if FixedParticleRadius < 1e-12, the function Chute::createBottom (called by Chute.cc) creates a smooth instead of rough bottom
	problem.setRoughBottomType(MONOLAYER_ORDERED); //bottom_type=0:grid wise bottom
	
	///values Thomas used in BegerSetup(ThomasValues_BegerSetup.cpp):
	//problem.setGravity(Vec3D(0,0,-9.8));
	//collision time 5 ms
	Mdouble tc = 50e-4;
	Mdouble eps = 0.57; //gemiddelde van run3 en run4 in air collision //thomas: 0.97
	Mdouble beta = 0.44;
	Mdouble mass = problem.species->getDensity()*mathsFunc::cubic(10e-3)*constants::pi/6.;
	problem.species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, eps, beta, mass);
	problem.species->setSlidingFrictionCoefficient(0.1);
	problem.species->setRollingStiffness(2./5.*problem.species->getStiffness());
	problem.species->setRollingDissipation(2./5.*problem.species->getDissipation());
	problem.species->setRollingFrictionCoefficient(0.);
	
	
	
	
	///Chute properties
	Mdouble chuteAngle = 0.0; // in degrees
	Mdouble chuteLength = 1.005;
	problem.setChuteAngle(chuteAngle);
	problem.setChuteLength(chuteLength);
	//if the particles are small enough to fit in, set the width to the real-life width (7.5 mm), else width=2.5*radius
	if (problem.getInflowParticleRadius()<7.5e-3*0.5) {
		problem.setChuteWidth(7.5e-3);
	} else {
		problem.setChuteWidth(problem.getInflowParticleRadius()*2.5); //7.5/3=2.5, dus zelfde verhouding als bij 3mm radius 7.5mm width
	}
	problem.setMaxFailed(100);
	problem.set_MaxFailedTriangle(1);
	
	///Hopper & Inset properties
	double hopperHeight, exitHeight, exitLength =0.1, hopperAngle = 60.0, hopperLength = 0.57; //hopperHeight = 0.535+0.13, exitHeight = 0.4
	//set_Hopper(0.01, 0.01, 60.0, 0.08, 0.02);
	problem.setHopperLift(0.0);
	
	//Referentie grootte inset en driehoek met deeltjes // DIT NIET AANPASSEN, zijn referentiewaarden
	Mdouble insetHeightRef=0.332;
	Mdouble insetWidthRef=0.21;
	Mdouble insetAngleRef=30*constants::pi/180.0;
	Mdouble tanInsetAngleRef=tan(insetAngleRef);
	Mdouble triangleHeightRef=insetHeightRef-insetWidthRef*tanInsetAngleRef;
	Mdouble triangleWidthRef=triangleHeightRef/tanInsetAngleRef;
	Mdouble aTriangleRef=0.5*triangleHeightRef*triangleWidthRef;
	
	//grootte van nieuwe inset en driehoek met deeltjes, de eis is dat de oppervlakte hetzelfde blijft
	Mdouble insetAngleDeg=10;
	Mdouble insetAngle=insetAngleDeg*constants::pi/180.0; //in radians //insetAngle=30*constants::pi/180.0;
	Mdouble insetWidth=insetWidthRef; 	//take insetWidth constant
	Mdouble aTriangle=aTriangleRef;		//take aTriangle constant
	Mdouble tanInsetAngle=tan(insetAngle);
	Mdouble triangleWidth=sqrt(2*aTriangleRef/tanInsetAngle); 
	Mdouble triangleHeight=triangleWidth*tanInsetAngle;
	
	//Adjust the insetHeight depending on the insetAngle, triangleHeight
	Mdouble insetHeight=triangleHeight+insetWidth*tanInsetAngle;
	
	Mdouble fallingHeight=0.3; //fallingHeight horende bij mijn opstelling is ongeveer 0.21
	hopperHeight=(hopperLength-exitLength)*0.5/tan(hopperAngle*constants::pi/180.0)+fallingHeight+insetHeight; //hopperHeight = 0.535+0.13
	exitHeight = hopperHeight-(hopperLength-0.5*(hopperLength-exitLength))/tan(hopperAngle/180*constants::pi)+tan(chuteAngle/180*constants::pi)*exitLength;//ik verwacht exitHeight >= 0.4
	cout<<"exitHeight "<<exitHeight<<endl;
	
	//Set the hopper and inset (input of angle in degrees!)
	problem.setHopper(exitLength,exitHeight,hopperAngle,hopperLength,hopperHeight);
	problem.set_Inset(insetHeight,insetWidth,insetAngleDeg);
	
	//Nu de chute langer maken als dat nodig is
	if (insetWidth+triangleWidth>chuteLength) {
		// nu vallen er een aantal zojuist gecreerde deeltjes al direct van de chute af
		Mdouble oversteek=(insetWidth+triangleWidth)-chuteLength;
		chuteLength=chuteLength+floor(oversteek*2.0)+1;
		problem.setChuteLength(chuteLength);
	} else {
		// nu vallen de gecreeerde deeltjes gewoon op de chute
		// hoeft niets te veranderen aan chuteLength
	}
	
	
	
	
	
	//Omdat de oppervlakte van de driehoek met deeltje hetzelfde blijft, blijft het aantal benodigde deeltjes ook gelijk (hoop ik!)
	Mdouble maxCreatedTriangle=1300;
	
	
	//Mdouble triangleHeight=insetHeight-insetWidth*tan(insetAngle);
	//Mdouble aTriangle=0.5*triangleHeight*triangleHeight/tan(insetAngle);
	//Mdouble insetHeightRef=0.332,insetWidthRef=0.21,insetAngleRef=30*constants::pi/180.0;
	
	//1400 deeltjes in triangle is best aardig voor bij inset: height=0.332m, width=0.21m, angle=30 degr.
	
	
	
	
	//double insetWidth=0.21, insetHeight=triangleHeight+insetWidth*tan(insetAngle); //insetHeight=0.332
	
	/*
	cout<<"triangleWidth "<<triangleWidth<<endl;
	cout<<"triangleHeight "<<triangleHeight<<endl;
	cout<<"aTriangle "<<aTriangle<<endl;
	cout<<"insetAngle "<<insetAngle/constants::pi*180.0<<" [degrees]"<<endl;
	cout<<"insetAngle "<<insetAngle<<" [radians]"<<endl;
	cout<<"tan(insetAngle) "<<tanInsetAngle<<endl;
	cout<<"insetWidth "<<insetWidth<<endl;
	cout<<"insetHeight "<<insetHeight<<endl;
	cout<<"maxCreatedTriangle "<<maxCreatedTriangle<<endl;
	cout<<"triangleHeightRef "<<triangleHeightRef<<endl;
	cout<<"triangleWidthRef "<<triangleWidthRef<<endl;
	cout<<"tan(insetAngleRef) "<<tanInsetAngleRef<<endl;
	cout<<"aTriangleRef "<<aTriangleRef<<endl;
	*/
	
	
	
	//int maxCreatedTriangle=1400; //350 is een best aantal voor problem.setInflowParticleRadius(5.7e-3,6.3e-3);
	int maxCreatedHopper=maxCreatedTriangle;
	problem.set_MaxCreatedTriangle(maxCreatedTriangle);
	problem.set_MaxCreatedHopper(maxCreatedHopper);
	
	
	///Turning settings
	problem.set_T1(tmax*0.1); // end filling time keep constant gravityangle and let the system relax
	problem.set_Angle1(-insetAngle/constants::pi*180); // angle of gravity during filling [degrees]
	problem.set_T2(tmax*0.2); // end relaxation time, start rotating
	problem.set_Angle2(0); // final angle 
	problem.set_T3(tmax*0.8); // end rotating, start adding flowing particles
	problem.set_T4(tmax); // end simulation time
	
	
	
	///time settings
	//speed allowed before particles move through each other!
	//cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl;
	//cout << "Maximum possible speed of particles: " << problem.getMaximumVelocityInducedByGravity() << endl;
	problem.setTimeStep(tc/50.0);
	problem.DPMBase::setTimeStep(problem.getTimeStep()*5.0);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(100,problem.getTimeMax(),problem.getTimeStep()));
	//problem.setSaveCount(1);
	cout << "dt=" << problem.getTimeStep() << endl;
	problem.solve();
	
	problem.writeRestartFile();
}
