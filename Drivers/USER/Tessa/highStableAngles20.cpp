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
#include "Walls/IntersectionOfWalls.h"

//based on hopper_inset_turning_cos_rolling_resistance_values_thomas5.cpp
//fallingHeight and inset Angle are the two variables the other dimensions (of inset and hopper) ares based on

class ChuteWithHopperAndInset_time : public ChuteWithHopperAndInset{
public:
	
	void setEKinThresholdFilling(Mdouble new_) {eKinThresholdFilling=new_;}
	Mdouble getEKinThresholdFilling() const {return eKinThresholdFilling;}
	
	void setEKinThresholdTurning(Mdouble new_) {eKinThresholdTurning=new_;}
	Mdouble getEKinThresholdTurning() const {return eKinThresholdTurning;}
	
	void setEKinThresholdHopper(Mdouble new_) {eKinThresholdHopper=new_;}
	Mdouble getEKinThresholdHopper() const {return eKinThresholdHopper;}
	
	void setEKinThresholdEnd(Mdouble new_) {eKinThresholdEnd=new_;}
	Mdouble getEKinThresholdEnd() const {return eKinThresholdEnd;}
	
	void setTurningTime(Mdouble new_) {turningTime=new_;}
	Mdouble getTurningTime() const {return turningTime;}
	

	void setStartWait1Time(Mdouble new_) {startWait1Time=new_;}
	Mdouble getStartWait1Time() const {return startWait1Time;}
	
	void setStartTurningTime(Mdouble new_) {startTurningTime=new_;}
	Mdouble getStartTurningTime() const {return startTurningTime;}
	
	void setStartWait2Time(Mdouble new_) {startWait2Time=new_;}
	Mdouble getStartWait2Time() const {return startWait2Time;}
	
	void setStartFillHopperTime(Mdouble new_) {startFillHopperTime=new_;}
	Mdouble getStartFillHopperTime() const {return startFillHopperTime;}
	
	void setStartWait3Time(Mdouble new_) {startWait3Time=new_;}
	Mdouble getStartWait3Time() const {return startWait3Time;}
	
	void setStartFinishedTime(Mdouble new_) {startFinishedTime=new_;}
	Mdouble getStartFinishedTime() const {return startFinishedTime;}
		
	void setCreateInHopper(unsigned int new_) {createInHopper = new_;}
	int getCreateInHopper() const {return createInHopper;}
	
	void setMaxCreatedHopper(unsigned int new_) {maxCreatedHopper = new_;}
	int getMaxCreatedHopper() const {return maxCreatedHopper;}
	
	void setNumCreatedHopper(unsigned int new_) {numCreatedHopper = new_;}
        //int getNumCreatedHopper() const {return numCreatedHopper;}

        int getNumCreatedHopper() const
        {
              const InsertionBoundary* b0= dynamic_cast<const InsertionBoundary*>(boundaryHandler.getLastObject());
              if(b0)
              {
                    return b0->getNumberOfParticlesInserted();
              }
              else
              {
                    return 0;
              }
        }

	
	void setMaxCreatedTriangle(unsigned int new_) {maxCreatedTriangle = new_;}
	int getMaxCreatedTriangle() const {return maxCreatedTriangle;}
	
	void setNumCreatedTriangle(unsigned int new_) {numCreatedTriangle = new_;}
	int getNumCreatedTriangle() const {return numCreatedTriangle;}
	
	void setMaxFailedTriangle(unsigned int new_) {maxFailedTriangle = new_;}
	
	void setAngle1(Mdouble new_) {angle1=new_;} //all input of angle is in degrees
	Mdouble getAngle1() const {return angle1;}
	
	void setAngle2(Mdouble new_) {angle2=new_;} //all input of angle is in degrees
	Mdouble getAngle2() const {return angle2;}
			
	void createInflowParticleTriangle()
	{
		//the following formula yields polydispersed particle radii:
	  inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(),getMaxInflowParticleRadius()));

		double insetHeight = get_InsetHeight();
		double insetWidth = get_InsetWidth();
		double insetAngle = get_InsetAngle(); //output is in radians
		double a = insetHeight-insetWidth*(sin(insetAngle)/cos(insetAngle)); //a=triangleHeight;
		double c= 0.5*(getHopperLength()-getHopperExitLength()); 
		
		double theta=0.5*constants::pi-insetAngle;
		
		double l=a/cos(theta);
		double h=sin(theta)*a; //a=triangleHeight;
		double b=cos(theta)*a;
		
		double x1,z1;
		do
		{
			x1=random.getRandomNumber(0.0,1.0)*l;
			z1=random.getRandomNumber(0.0,1.0)*h;
		} while (z1>h/b*x1 || z1>h/(l-b)*(l-x1));
		//whilevoorwaarde: deeltje ligt boven de lijn van (b,h) naar (l,0): z=h/(l-b)*(l-x) of
		//deeltje ligt boven de lijn van (0,0) naar (b,h): z=h/b x, dan zit hij in het verkeerde
		//gebied en moet een andere 'random' positie gezocht worden
		
		//als hij hier aankomt zit hij in de driehoek en mag het deeltje gecreerd worden (gebeurt in IsInsertable(P0))
        Vec3D position;
        position.X = x1 * sin(theta) + z1 * cos(theta) + insetWidth + c;
        position.Y = getMaxInflowParticleRadius();
        position.Z = -x1 * cos(theta) + z1 * sin(theta) + a;
        // shifting Z
        position.Z = position.Z + getHopperLift();
        
        inflowParticle_.setPosition(position);

	//		P0.setPosition( Vec3D(
	//			 x1*sin(theta)+z1*cos(theta) + insetWidth + c , 
	//		getMaxInflowParticleRadius() , 
	//		-x1*cos(theta)+z1*sin(theta) + a));
	//	
	//	P0.getPosition().Z+=getHopperLift();
	}
	
	//change add_particles so there is a phase where no particles are added 
  //	void addParticlesHopper() 
  //	{
  //		int failed = 0;
  //		
  //		//try max_failed times to find new insertable particle in the hopper
  //		while (failed<=max_failed && numCreatedHopper < maxCreatedHopper){
  //			createInflowParticle();
  //			if (IsInsertable(P0)) {
  //				failed = 0; 
  //				num_created++;
  //				numCreatedHopper++;
  //			} else failed++;
  //		}
  //	}
	
	void addParticlesTriangle() 
	{
		int failedTriangle=0;
		
		//try maxFailedTriangle times to find new insertable particle in triangle
		while (failedTriangle <= maxFailedTriangle && numCreatedTriangle < maxCreatedTriangle)                {
			createInflowParticleTriangle();
                        if (checkParticleForInteraction(inflowParticle_))
                        {
                             particleHandler.copyAndAddObject(inflowParticle_);
                             failedTriangle = 0;
                             ++numCreatedTriangle;
                        }
                        else
                        {
                             failedTriangle++;
                        } 

			//			if (IsInsertable(P0)) {
			//	failedTriangle = 0;
			//	num_created++;
			//	numCreatedTriangle++;
			//} else failedTriangle++;
		}
	}
	
	double calcEKin() const
	{
		double eKin=0;
		for (unsigned int i=0;i<particleHandler.getNumberOfObjects();i++)
		{
			if (!particleHandler.getObject(i)->isFixed())
				eKin += .5 * particleHandler.getObject(i)->getMass() * particleHandler.getObject(i)->getVelocity().getLengthSquared();
		}
		return eKin;
	}
	
	void actionsBeforeTimeStep(){
		Mdouble gravity = getGravity().getLength();;
		Mdouble GravityAngle;	
		
		///Switch with all six statuses
		switch (currStatus) {
			case fillTriangle:
				{
					GravityAngle=angle1*constants::pi/180.0;
					setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
					addParticlesTriangle();
					if (numCreatedTriangle >= maxCreatedTriangle) {
						setStartWait1Time(getTime());
						currStatus=wait1;
					}
				}
				break;
			case wait1:
				if (calcEKin()<getEKinThresholdFilling()) 
				{
					currStatus=turning;
					setStartTurningTime(getTime());
				}
				break;
			case turning:
				{
					//turning with a cosine
					//check time elapsed 
					Mdouble time,t1,t2;
					time =getTime();
					t1 = getStartTurningTime();
					t2=t1+getTurningTime();
					if (time<=t2)
					{ //turn
						GravityAngle=((angle2-angle1)*(1-cos((time-t1)/(t2-t1)*constants::pi))*0.5+angle1)*constants::pi/180.0;
						setGravity(Vec3D(sin(GravityAngle), 0.0, -cos(GravityAngle))*gravity);
					}
					else
					{
						setStartWait2Time(t2);
						currStatus=wait2;
					}
				}
				break;
			case wait2:
				if (calcEKin()<getEKinThresholdTurning()) 
				{
					setStartFillHopperTime(getTime());
					currStatus=fillHopper;
				}
				break;
			case fillHopper:
				{
					if (calcEKin()<getEKinThresholdFilling()&&getNumCreatedHopper()>0)
					{
						setStartWait3Time(getTime());
						currStatus=wait3;
                                                boundaryHandler.removeLastObject();
                                                removeExtraHopperWall();

						//removeExtraHopperWall();
						
						//added to finish the simulation after the turning
						//currStatus=finished;
					}
					//if(getNumCreatedHopper() < getMaxCreatedHopper())
					//{
					//	addParticlesHopper();
					//}
				}
				break;
			case wait3:
				if (calcEKin()<getEKinThresholdEnd()) 
				{
					setStartFinishedTime(getTime());
					currStatus=finished;
				}
				break;
		       case finished:
		                break;
		}
		cleanChute();
	}	
	
	void actionsBeforeTimeLoop() {
		
		//automatically sets dt if dt is not specified by the user
	        //if (!getTimeStep()) setTimeStep(); 
				
            	std::ofstream myfile;
		myfile.open ("gravity_file.txt");
		myfile<<"Constants: "<<"angle1:"<<getAngle1()<<", angle2:"<<getAngle2()<<", turningTime:"<<getTurningTime()<<std::endl;
		myfile.close();
		
		//create extra finite wall to keep the particles in the hopper
		createExtraHopperWall();
	};
	
	void createExtraHopperWall() {
	       //int n = getNWall();
	       //setNWall(n+1);
		
		//double s = sin(getChuteAngle());
		//double c = cos(getChuteAngle());
		
		Vec3D A, B, C, temp, normal;
		
		///The points A, B, C should be the corners of the triangular 
		///finite wall in clockwise order. 
		
		///Define the three points between which the finite extra hopper wall is created
		A = Vec3D(0,0,getHopperHeight()-0.5*(getHopperLength()-getHopperExitLength())/tan(getHopperAngle()));
		B = Vec3D(getHopperExitLength(),0.0,getHopperHeight()-0.5*(getHopperLength()-getHopperExitLength())/tan(getHopperAngle()));
		C = Vec3D(0.5*getHopperExitLength(),0.0,getHopperHeight()-0.5*(getHopperLength()-getHopperExitLength())/tan(getHopperAngle())-0.5*getHopperExitLength());
		
		///Move points A,B,C a distance 'shift' down the chute, so the extra hopper wall starts at the beginning of the chute
		A.X+= getHopperShift();
		B.X+= getHopperShift();
		C.X+= getHopperShift();
                IntersectionOfWalls w0;
                ///create a finite wall from B to A, from C to B and from A to C
                temp = B - A;
                normal = Vec3D(temp.Z, 0.0, -temp.X) / std::sqrt(temp.getLengthSquared());
                w0.addObject(normal, A); //Walls[n].addObject(normal, Dot(normal,A));
                temp = C - B;
                normal = Vec3D(temp.Z, 0.0, -temp.X) / std::sqrt(temp.getLengthSquared());
                w0.addObject(normal, B); //Walls[n].addObject(normal, Dot(normal,B));
                temp = A - C;
                normal = Vec3D(temp.Z, 0.0, -temp.X) / std::sqrt(temp.getLengthSquared());
                w0.addObject(normal, C); //Walls[n].addObject(normal,Dot(normal,C));
                wallHandler.copyAndAddObject(w0);
		
		///create a finite wall from B to A, from C to B and from A to C
		//		temp = B-A;
		//normal  = Vec3D(temp.Z,0.0,-temp.X) / sqrt(temp.getLengthSquared());
		//Walls[n].addObject(normal, Dot(normal,A));
		//temp = C-B;
		//normal  = Vec3D(temp.Z,0.0,-temp.X) / sqrt(temp.getLengthSquared());
		//Walls[n].addObject(normal, Dot(normal,B));
		//temp = A-C;
		//normal = Vec3D(temp.Z,0.0,-temp.X)/sqrt(temp.getLengthSquared());
		//Walls[n].addObject(normal,Dot(normal,C));
		
	}
	void removeExtraHopperWall() {
                wallHandler.removeLastObject();
//		int n = getNWall();
//		setNWall(n-1);
	}
	
	void actionsAfterTimeStep(){
		//save the gravity information in "gravity_file.txt" after every factor'th time step
		Mdouble factor = 10; //100
		
		Mdouble time= getTime();
		Mdouble dt= getTimeStep();
		
		if (std::abs(fmod(time+0.5*dt,factor*dt)-0.5*dt)<dt*1e-2) {
			//ugly fix to be able to compare two doubles
			//(als je door een afrondingsfout net onder factor*dt zit, komt de modulus ervan niet net onder nul te liggen, als je er eerst 0.5*dt bij optelt, dan modulus doet en daarna weer 0.5*dt van aftrekt wel)
			Vec3D gravity = getGravity();
			std::fstream myfile;
			myfile.open ("gravity_file.txt",std::fstream::in | std::fstream::out | std::fstream::app);
			myfile <<"time:"<<std::fixed<<std::setprecision(5)<<std::setw(9)<<time<<", g:"<<std::setw(17)<<std::setprecision(2)<<gravity<<std::endl;
			myfile.close();
		} 	
	}
	
	virtual void actionsAfterSolve(){
	        std::ofstream myfile;
		myfile.open ("times.txt");
		myfile<<"This file includes all starting/end times of all six stages"<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startFillTriangleTime: "<<0.0<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startWait1Time: "<<getStartWait1Time()<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startTurningTime: "<<getStartTurningTime()<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startWait2Time: "<<getStartWait2Time()<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startFillHopperTime: "<<getStartFillHopperTime()<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startWait3Time: "<<getStartWait3Time()<<std::endl;
		myfile<<std::fixed<<std::setw(25)<<"startFinishedTime: "<<getStartFinishedTime()<<std::endl;
		myfile.close();				
	};
	
	bool continueSolve() const {
		if (currStatus ==finished) {
			return false;
		} else {
			return true;
		}
	}
	
	void printTime() const {
	  std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() 
		    << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
		    << ", N=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects();
			switch (currStatus) {
				case fillTriangle:
				  std::cout<<"Filling triangle, number triangle particles = "<<std::setw(4)<<numCreatedTriangle <<" (maxCreatedTriangle = "<<maxCreatedTriangle<<")";
					break;
				case wait1:
				  std::cout<<"wait for eKin ("<<std::setw(9)<<calcEKin()<<") to drop below eKinThresholdFilling ("<<getEKinThresholdFilling()<<")";
					break;
				case turning:
				  std::cout<<"turning, end turning time is "<<getStartTurningTime()+getTurningTime();
					break;
				case wait2:
				  std::cout<<"wait for eKin ("<<std::setw(9)<<calcEKin()<<") to drop below eKinThresholdTurning ("<<getEKinThresholdTurning()<<")";
					break;
				case fillHopper:
				  std::cout<<"Filling hopper, number hopper particles = "<<std::setw(4)<<numCreatedHopper <<" (maxCreatedHopper = "<<maxCreatedHopper<<")";
					break;
				case wait3:
				  std::cout<<"wait for eKin ("<<std::setw(9)<<calcEKin()<<") to drop below eKinThresholdEnd ("<<getEKinThresholdEnd()<<")";
					break;
              			case finished:
			                break;
			}
			std::cout<< std::endl;
	}
	
	virtual double getInfo(BaseParticle& P UNUSED) {
		return P.getId();
	}
	
	private:
	int numCreatedHopper;
	int numCreatedTriangle;
	int maxFailedTriangle;
	int createInHopper;
	int maxCreatedHopper;
	int maxCreatedTriangle;
	enum status
	{
		fillTriangle,
		wait1,
		turning,
		wait2,
		fillHopper,
		wait3,
		finished
	};
	status currStatus;
	
	Mdouble angle1,angle2;
		
	Mdouble eKinThresholdFilling,eKinThresholdTurning,eKinThresholdHopper,eKinThresholdEnd;
	Mdouble turningTime;
	
	Mdouble startWait1Time,startTurningTime,startWait2Time;
	Mdouble startFillHopperTime,startWait3Time,startFinishedTime;
  SphericalParticle inflowParticle_;
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	ChuteWithHopperAndInset_time problem;
	///Initialisatie
	problem.setNumCreatedTriangle(0);
	
	///Problem parameters
	problem.setName("highStableAngles");
	Mdouble tmax=50.0;
	problem.setTimeMax(tmax);
	
	///Particle properties
	problem.species->setDensity(2500.0); //uit data sheet M type beads
	//problem.setInflowParticleRadius(5.7e-3,6.3e-3);
	problem.setInflowParticleRadius(5.7e-3/2.0,6.3e-3/2.0);
	
	// set a bottom made of regular placed particles
	problem.setFixedParticleRadius(problem.getInflowParticleRadius()); // if FixedParticleRadius < 1e-12, the function Chute::createBottom (called by Chute.cc) creates a smooth instead of rough bottom
	problem.setRoughBottomType(MULTILAYER); //bottom_type=0:grid wise bottom, 1=random (thin), 2=random (thick)
	
	///values Thomas used in BegerSetup(ThomasValues_BegerSetup.cpp):
	//problem.setGravity(Vec3D(0,0,-9.8));
	//collision time 5 ms
	Mdouble tc = 50e-4;
	Mdouble eps = 0.57; //gemiddelde van run3 en run4 in air collision //thomas: 0.97
	Mdouble beta = 0.44;
	Mdouble mass = problem.species->getDensity()*mathsFunc::cubic(10e-3)*constants::pi/6.;
	problem.species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, eps, beta, mass);
	problem.species->setSlidingFrictionCoefficient(0.1); //T: 0.1
	problem.species->setRollingStiffness(2./5.*problem.species->getStiffness());//T:2./5.*problem.species->getStiffness()
	problem.species->setRollingDissipation(2./5.*problem.species->getDissipation());//T:2./5.*problem.species->getDissipation()
	problem.species->setRollingFrictionCoefficient(0.0); //T: 0.
	
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
	problem.setMaxFailedTriangle(1);
	
	///Hopper & Inset properties
	double hopperHeight, exitHeight, exitLength =0.1, hopperAngle = 60.0, hopperLength = 0.57; //hopperHeight = 0.535+0.13, exitHeight = 0.4
	//setHopper(0.01, 0.01, 60.0, 0.08, 0.02);
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
	Mdouble insetAngleDeg=20;
	Mdouble insetAngle=insetAngleDeg*constants::pi/180.0; //in radians //insetAngle=30*constants::pi/180.0;
	Mdouble insetWidth=insetWidthRef; 	//take insetWidth constant
	Mdouble aTriangle=aTriangleRef;		//take aTriangle constant
	Mdouble tanInsetAngle=tan(insetAngle);
	Mdouble triangleWidth=sqrt(2*aTriangleRef/tanInsetAngle); 
	Mdouble triangleHeight=triangleWidth*tanInsetAngle;
	
	//Adjust the insetHeight depending on the insetAngle, triangleHeight
	Mdouble insetHeight=triangleHeight+insetWidth*tanInsetAngle;
	
	Mdouble fallingHeight=0.3; //fallingHeight horende bij mijn opstelling is ongeveer 0.21
	hopperHeight=(hopperLength-exitLength)*0.5/tan(hopperAngle*constants::pi/180.0)+fallingHeight+insetHeight-0.5*insetWidth*tan(insetAngle); //hopperHeight = 0.535+0.13
	exitHeight = hopperHeight-(hopperLength-0.5*(hopperLength-exitLength))/tan(hopperAngle/180*constants::pi)+tan(chuteAngle/180*constants::pi)*exitLength;//ik verwacht exitHeight >= 0.4
	std::cout<<"exitHeight "<<exitHeight<<std::endl;
	
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
	//int maxCreatedTriangle=1400; //350 is een best aantal voor problem.setInflowParticleRadius(5.7e-3,6.3e-3);
	int maxCreatedHopper=maxCreatedTriangle;
	problem.setMaxCreatedTriangle(maxCreatedTriangle);
	problem.setMaxCreatedHopper(maxCreatedHopper);
	
	
	///Turning settings
	problem.setEKinThresholdFilling(1e-12); //1e-12 is the value thomas used in BegerSetup
	problem.setEKinThresholdTurning(1e-12); //1e-12 is the value thomas used in BegerSetup
	problem.setEKinThresholdHopper(1e-12);  //1e-12 is the value thomas used in BegerSetup
	problem.setEKinThresholdEnd(1e-12);     //1e-12 is the value thomas used in BegerSetup
	problem.setTurningTime(10.0); //eenheid in seconden real time
	problem.setStartTurningTime(0.0); //initialiseren op nul, wordt vanzelf verandert
	problem.setAngle1(-insetAngle/constants::pi*180); // angle of gravity during filling [degrees]
	problem.setAngle2(0); // final angle 
	
	///time settings
	//speed allowed before particles move through each other!
	//	std::cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << std::endl; 
	//std::cout << "Maximum possible speed of particles: " << problem.getMaximumVelocityInducedByGravity() << std::endl; 
	problem.setTimeStep(tc/50.0);
	problem.DPMBase::setTimeStep(problem.getTimeStep()*5.0);
	//problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(500,problem.getTimeMax(),problem.getTimeStep()));
	problem.setSaveCount(100);
	std::cout << "dt=" << problem.getTimeStep() << std::endl;
	problem.solve();
	
	problem.writeRestartFile();
}
