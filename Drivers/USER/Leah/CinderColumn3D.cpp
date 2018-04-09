#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "File.h"
#include "Walls/InfiniteWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/IntersectionOfWalls.h"
#include "Particles/BaseParticle.h"

//NOTE: power function for squared numbers can be replaced with: mathsFunc::square(x), mathsFunc::cubic(x) 

using namespace std;	

class CinderDriver : public Mercury3D {


public:

double getInfo(BaseParticle& P) {return P.getIndex();}


// ******************** SET DEFAULTS, overwritten by set function (invisible to user)************//
CinderDriver()						
	{
	verif_file.open("verification_file.dat");
	NumEjected=0;                                           // Number of particles ejected
	count=0;                                                // counts nuber of timesteps to determine when to delete 'landed' particles
	nburst=0;                                               // initial number of bursts that have ejected
	
	// set seed for RNG
	double seed = (time(NULL));                             // set seed for RNG
	random.setRandomSeed(time(NULL));                       // random integers used to determine particle ejection positions and velocities		
    	}
    	
//********************************end default values***************************//



// **************** Define smallest and largest BaseParticles possible based on user-set initial conditions ********************//
BaseParticle* getSmallestParticle(){
	BaseParticle P0;
	P0.setRadius(getMinParticleRadius());
	BaseParticle* Pp=&P0;
	return Pp;
}


BaseParticle* getLargestParticle(){
	BaseParticle P0;
	P0.setRadius(getMaxParticleRadius());
	BaseParticle* Pp=&P0;
	return Pp;
}



// ************************* SET INITIAL CONDITIONS, FILES TO WRITE *****************************//

void setupInitialConditions(){	

	// Create a cylindrical wall:
    AxisymmetricIntersectionOfWalls w0;
    w0.setPosition(Vec3D(0,0,0)); 
    w0.setOrientation(Vec3D(0, 0, 1));
            
    std::vector<Vec3D> Points(3);
    Points[0] = Vec3D(getVentDiameter()*0.5001,  0,  getVentElevation()+getConduitDepth());
    Points[1] = Vec3D(getVentDiameter()/2,0,  getVentElevation()+getConduitDepth());
    Points[2] = Vec3D(getVentDiameter()/2,0,  getVentElevation());

    w0.createOpenPrism(Points);
    wallHandler.copyAndAddObject(w0);   
      
       
    // Create file to store landing positions of particles: CinderColumn3D.N.fpos   <--- does't work
    std::string name=getName();     
    int number=getRunNumber();
       
    fpos.setName(getName()+".fpos");
    fpos.setOpenMode(std::fstream::out | std::fstream::app);
        
     

        
    // write input variables to file CinderDriver.X.restart
    writeRestartFile();		
	}	
	

    // Write collisions information to fstat files
    void writeFstatHeader(std::ostream& os) const 
    {
        for (std::vector<BaseInteraction*>::const_iterator it = interactionHandler.begin(); it != interactionHandler.end(); ++it)
            {
                //(*it)->writeToFStat(os);
                os << (*it)->getTimeStamp()
                    << " " << (*it)->getP()->getIndex()
                    << " " << (*it)->getI()->getIndex()
                    << " " << (*it)->getOverlap()
                    << " " << (*it)->getContactPoint()
                    << " " << (*it)->getForce()
                    << " " << (*it)->getRelativeVelocity()
                    << std::endl;
                }
             }
    
    // Add particle number to data file
    double getInfo(const BaseParticle& P) const
    {
        return P.getIndex();
    }

    // Print time to screen
    void printTime() const
    {
        std::cout << "\rt=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
            << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax() << "\r";
        std::cout.flush();
        
        
    //   actionsBeforeTimeStep();          UNCOMMENT TO RUN XBALLS!!
    //   BaseParticle P0;
    //   P0.setPosition(Vec3D(45,45,700));
    //   P0.setRadius(getMinParticleRadius());
    //   particleHandler.copyAndAddObject(P0);
    }
    
// ***************************end set initial conditions***************************//



// *************************** INSERT PARTICLES PARAMETERS ***********************//

	
	void actionsBeforeTimeStep()
	{
	
	 DPMBase::actionsBeforeTimeStep();
	
	
// ***************************************************** Insert Particles  **********************************************

		double t=getTime();                         // time
		double bt=getBurstTime();                   // burst duration
		double rt=getRestTime();                    // time between bursts

	    if (nburst < getMaxNburst()) {              // if number of bursts is less than total number of bursts to be ejected
	   	        if (nburst == 0.0 || bt == 0 || (t > nburst*(bt+rt) && nburst < t/(bt+rt)) ) {
	   	            nburst=nburst+1;
			        int i=0;
			
			        while (i<getEjectionRate()*getBurstTime()) {
        	            BaseParticle P0;
			            P0.setRadius(calcParticleRadiusPower());                                                // set particle radius
  	    	            P0.setPosition(calcParticlePositionBubble(getVentElevation()-P0.getRadius()));  // set position in 3 dimensions  	
    		            P0.setVelocity(calcVelocityBurstDecay(P0.getRadius(),P0.getPosition()));    			// set 3D velocity

			            //Check if last inserted particle is overlaping and insert particle
			            P0.setHandler(&particleHandler);

			            if (checkParticleForInteraction(P0))
				            {
			                particleHandler.copyAndAddObject(P0);
				            NumEjected++;
				            i++;
				        }
		            }
		            
		            std::cout << "\nNumEjected: " << NumEjected <<  "\n";
	            }	
        }
        
        
}



	void actionsAfterTimeStep()
		{
			
			
		// every 1000 timesteps, delete particles that fall below the vent elevation //(should be landing surface, problem with walls) 
		// save particle positions to file CinderColumn.X.fpos
		double MaxN = getMaxNburst()*getBurstTime()*getEjectionRate();
		count=count+1;
		if (count == 1000){
				count = 0;
				
				for (int i=0; i < particleHandler.getNumberOfObjects();i++) { 
			//	for (int i=0; i < 100; i++) { 
                    if (particleHandler.getObject(i)) {
				    double X=particleHandler.getObject(i)->getPosition().X;
				    double Y=particleHandler.getObject(i)->getPosition().Y;
				    double Xa=abs(X);
				    double Ya=abs(Y);
				    double Z=particleHandler.getObject(i)->getPosition().Z;
				    double Vz=particleHandler.getObject(i)->getVelocity().Z;
			
			
			
			
				    
				    // Remove particles falling into the conduit
				/*	    if (Z < getVentElevation()+getMaxParticleRadius() && Vz < 0.0  && Y < getVentDiameter()/2  && X < getVentDiameter()/2) {		
						particleHandler.removeObject(i);	
						 }
                */


					 
			/*		 // Make fixed particles that have landed
					 if (Z < getVentElevation()+getConduitDepth()+getMaxParticleRadius()*2 && Vz < 0 && (Ya > getVentDiameter()/2 || Xa > getVentDiameter()/2)) {		
	                    particleHandler.getObject(i)->setPosition(Vec3D(X,Y,getVentElevation()+getConduitDepth()+getMaxParticleRadius()*2));
						//particleHandler.getObject(i)->setVelocity(Vec3D(0,0,0);
						particleHandler.getObject(i)->fixParticle();
					 }
					 
			*/		 
					// Remove particles falling to the ground
				//    if (Z < getVentElevation()+getConduitDepth()+getMaxParticleRadius()*2+3 && (Vz < 0.0) && (Y > getVentDiameter()/2 || X > getVentDiameter()/2)) {		
		            if (Z < getVentElevation()+getConduitDepth()+getMaxParticleRadius()*2+3 && (Vz < 0.0) && (Ya > getVentDiameter()/2 || Xa > getVentDiameter()/2)) {		
	
					/*	fpos.getFstream() << particleHandler.getObject(i)->getPosition().X
						<< " " << particleHandler.getObject(i)->getPosition().Y
						<< " " << particleHandler.getObject(i)->getPosition().Z
						<< " " << particleHandler.getObject(i)->getVelocity().X
						<< " " << particleHandler.getObject(i)->getVelocity().Y
						<< " " << particleHandler.getObject(i)->getVelocity().Z
						<<std::endl;		*/
						
						particleHandler.removeObject(i);	}
					    } 
					} 
			}
		}
		
		
// ****************************end particle parammeters ***********************//	




void actionAfterSolve()
		{
							
fpos.close();

}



//****************** EXTERNAL FORCES i.e. drags, gravity etc...***************//

	void computeExternalForces(BaseParticle* P0)
	
{	
	// Call the MD compute_external_forces function (turns on gravity)
	DPMBase::computeExternalForces(P0);     // delete this line to turn off mercury gravity - turns off more than that!!

	// Drag Force - constant Cd
	double Po;              // part of pressure calculation
	double Temp;            // Temperature in the atmosphere
	double Pressure;        // Atmospheric pressure
	double density_air;     // Air density
	double Cd;              // Drag coefficient (constant)
	double A;               // Cross sectional area of particle
	double Fx;              // Force in the x dimension
	double Fy;              // Force in the y dimension
	double Fz;              // Force in the z dimension: note that this includes gravity
	double v;               // Particle velocity
	double mass;            // Particle mass

    // Drag equations, from Eject!
	A=3.1415*pow(P0->getRadius(),2);
	Cd=getCdconst();
	v=pow(pow((P0->getVelocity().X-getVwind()),2)+pow(P0->getVelocity().Y,2)+pow(P0->getVelocity().Z,2),0.5);	
	Po=101300.0*pow(((getTzeroK()-(getVentElevation())*getLapse()/1000.0)/getTzeroK()),(-9.81/(getRair()*getLapse()/1000.0)));
	Temp=getTzeroK()-(P0->getPosition().Z+getVentElevation())*getLapse()/1000.0;
	Pressure=Po*pow(Temp/getTzeroK(),9.81/(getRair()*getLapse()/1000.0));
	density_air=Pressure/(getRair()*Temp);
	mass=P0->getMass();
	
	// IS THIS THE DRAG FORCE OR THE TOTAL FORCE?
	Fx=-((P0->getVelocity().X-getVwind())*density_air*v*A*Cd)/(2);        // Wind is in the x direction
	Fy=-(P0->getVelocity().Y*density_air*v*A*Cd)/(2);
	Fz=-(P0->getVelocity().Z*density_air*v*A*Cd)/(2)+getLGravity().Z*((getDensity()-density_air)/getDensity())*P0->getMass();    //+ because gravity is negative                             

    // Force within the conduit: Equation 6.17
    // NOTE: NEED TO DECAY WITH DISTANCE!!!!  -- working only inside conduit?
    // NOTE: now particles will exit too quickly!!  (experience force within the conduit but already had their final velocities)
     if (P0->getPosition().Z < getVentElevation()+getConduitDepth() && abs(P0->getPosition().X) < getVentDiameter()/2 && abs(P0->getPosition().Y) < getVentDiameter()/2 ) {
        Fz=Fz+getGasDensity()*getCdconst()*A/2*pow((getGasVelocity()-P0->getVelocity().Z),2)-mass*((getDensity()-getGasDensity())/getDensity())*getLGravity().Z;
        }
    
    // Commented out to verify code against eject !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	//P0->addForce(Vec3D(Fx, Fy, Fz));           // REPLACE GRAVITY. . . otherwise we have gravity twice??
	

   
   
    /* Variable Drag (Not coded, requires tables, etc.)	
	//c = 20.116*pow(Temp,0.5)                                       
	//visc=(0.000172*(390/(Temp+117))*pow(Temp/273,1.5))/10;   
	//vwind = pow(pow(getVelocity().Z,2)+pow(getVelocity().X-w),2),0.5);     
	//reynolds = rhoa*vwind*2*radius/visc;
    //mach = vwind/c;                       
    */
    
      
}
	

//**************************end external forces*******************************//



//***************************CALCULATION FUNCTIONS****************************//

//Calculate initial particle position in the vent, Gaussian distibution (more particles in the center)
Vec3D calcParticlePositionBubble(double ZposA){

	double negpos, np, np1, np2, np3, np4, RN1, RN2, RN3, RN4, RN5, RN6, RN7, RN8, RN9, RN10, RN11, RN12, RN;   // random number declarations
	negpos=random.getRandomNumber(-1,1);        //randomly (+1) or (-1)
	np1=negpos/abs(negpos);
	negpos=random.getRandomNumber(-1,1);
	np2=negpos/abs(negpos);
	double Xpos, Ypos, Zpos;
				
	//Gaussian Distribution RNG: Random number between 0 and 0.5???
	RN1=random.getRandomNumber(0,1);
	RN2=random.getRandomNumber(0,1);
	RN3=random.getRandomNumber(0,1);
	RN4=random.getRandomNumber(0,1);
    RN5=random.getRandomNumber(0,1);
	RN6=random.getRandomNumber(0,1);
	RN7=random.getRandomNumber(0,1);
	RN8=random.getRandomNumber(0,1);
	RN9=random.getRandomNumber(0,1);
	RN10=random.getRandomNumber(0,1);
	RN11=random.getRandomNumber(0,1);
	RN12=random.getRandomNumber(0,1);
    RN=((RN1+RN2+RN3+RN4+RN5+RN6+RN7+RN8+RN9+RN10+RN11+RN12)-6)/6;  
       
	
	Xpos = RN1*0.5*getVentDiameter()*np1;
	Ypos = RN2*pow(abs(Xpos*Xpos-0.5*getVentDiameter()*0.5*getVentDiameter()),0.5)*np2;
	Zpos = ZposA+pow(0.5*getVentDiameter()*0.5*getVentDiameter()-Xpos*Xpos-Ypos*Ypos,0.5);	
	
	//Xpos=0;
	//Ypos=1;
	//Zpos=0;

	return Vec3D(Xpos,Ypos,Zpos);
	            }



//Calculate initial particle velocity: Bubble Burst
Vec3D calcVelocityBurstDecay(double radius, Vec3D position){

		double p=getPSDn()+1;                     // Particle size distribution exponent
		double temp;                              // parameter in particle size distribution calculation
		double DistFromCenter=pow((position.X*position.X)+(position.Y*position.Y)+((position.Z-getVentElevation())*(position.Z-getVentElevation())),0.5);
		
		// set 3D velocity
		// speed coupled to the gas jet (max speed for particle size)
		double Speed = getGasVelocity()-pow(4*9.81/(3*getCdconst())*(getDensity()/getGasDensity()-1)*(2*radius),0.5); 	
        double SpeedX=Speed*position.X/DistFromCenter;
        double SpeedY=Speed*position.Y/DistFromCenter;
        double SpeedZ=abs(Speed*(position.Z-getVentElevation())/DistFromCenter);

		Vec3D Speed3D = Vec3D(SpeedX,SpeedY,SpeedZ);    // actually velocity not speed
		
		//Single Speed
		//Speed3D = Vec3D(cos(45*constants::pi/180.0),0,sin(45)*constants::pi/180.0);
		
		
		return Speed3D;
		}
	
//Calculate random velocities within user-suplpied range
Vec3D calcVelocityRandom(){
			double negpos, np1;

			//randomly (+1) or (-1)
			//random.setRandomSeed(time(NULL));
			negpos=random.getRandomNumber(-1,1);
			np1=negpos/abs(negpos);
	
			//calculate angle in 3 dimensions
			double angfromvert=(90.0-getEjectionMedAngle())*np1+90.0;;
			double Angle = (angfromvert+getEjectionAngleVariance()*random.getRandomNumber(-1,1))*constants::pi/180.0;	
			double AngleXY=random.getRandomNumber(0.0,2.0*constants::pi);

			//calculate 3D velocity
			double Speed = getEjectionMedVelocity()+getEjectionVelocityVariance()*random.getRandomNumber(-1,1);	
			double SpeedX = Speed*cos(Angle)*cos(AngleXY);
			double SpeedY = Speed*cos(Angle)*sin(AngleXY);	
			return Vec3D(SpeedX,SpeedY,Speed*sin(Angle));
			}
			
			
//Calculate particle size based on a power distribution			
double calcParticleRadiusPower(){
		double rmin=getMinParticleRadius();        // min possible radius
		double rmax=getMaxParticleRadius();        // max possible radius
		double p=getPSDn()+1;                      // Particle size distribution exponent
		double temp;                               // parameter in particle size distribution calculation

		temp=(pow(rmax,p)-pow(rmin,p))*random.getRandomNumber(0.05,1)+pow(rmin,p);      // particle radius: PSD distrubution
		return pow(temp,(1/p));	
		//return 10;
		}
		
//Calculate particle size based on a linear distribution
double calcParticleRadiusLinear(){
        return random.getRandomNumber(MinParticleRadius,MaxParticleRadius);	
       }


void actionsAfterSolve()	// After model completes running
	{
	verif_file.close();
	}






// ********************** VARIABLE FUNCITONS ****************//

// Set Functions

    // interaction parameters
    void setTc(double t){Tc=t;}
	void setRest(double r){Rest=r;}
    
    // particle parameters
	void setDensity(double R){ speciesHandler.getObject(0)->setDensity(R); }
	void setParticleRadius(double rmin,double rmax){
		if (rmin>rmax)	{
			cerr << "Your maximum radius is less than your minimum. Particle radius has not been set" << endl;
			return; }
		if (rmin<0) {
			cerr << "Your minimum radius is less than 0, radius has not been set"<< endl;		
			return;	}
		if (rmax<0) {
			cerr << "Your minimum radius is less than 0, radius has not been set"<< endl;
			return;	}
		MaxParticleRadius=rmax;
		MinParticleRadius=rmin;
		}
		void setCdconst(double cconst){ Cdconst=cconst;}
	    void setPSDn(double n){PSDn=n;}
	
	void setEjectionMedVelocity(double vel){ EjectionMedVelocity=vel; }
	void setEjectionVelocityVariance(double vel){ EjectionVelocityVariance=vel; }
    void setEjectionVelocityMax(double vel){ EjectionVelocityMax=vel;}
    void setEjectionVelocityMin(double vel){ EjectionVelocityMin=vel;}
    void setEjectionVelocityFactor(double f){EjectionVelocityFactor=f;}
   	void setEjectionMedAngle(double new_vel_ang){ EjectionMedAngle=new_vel_ang; }
	void setEjectionAngleVariance(double ang_var){ EjectionAngleVariance=ang_var; }
	void setEjectionRate(double rate){EjectionRate=rate;}
	void setMinAngle(double angle){MinAngle=angle;}

	void setBurstTime(double t){BurstTime=t;}
	void setRestTime(double rt){RestTime=rt;}
	void setnburstT(int n) {nburstT=n; nburst=0;} 
	void setMaxNburst(int n){MaxNburst=n;}

    void setSpeedMode(int speed){SpeedMode=speed;}
    void setAngleMode(int ang){AngleMode=ang;}
    void setVelocityMode(int v){VelocityMode=v;}
	void setGasDensity(double gd) {GasDensity=gd;}
	void setGasVelocity(double gv) {GasVelocity=getEjectionMedVelocity()+pow(4/3*9.81*getCdconst()*(1250/getGasDensity()-1)*(2*0.1),0.5);}	
	
	// vent dimensions
	void setVentElevation(double event){ VentElevation=event-getConduitDepth();}
	void setVentDiameter(double new_dim){ VentDiameter=new_dim;}
	void setConduitDirection(Vec3D dir){ ConduitDirection=dir; }
	void setConduitDepth(double depth){ ConduitDepth=depth; }
	
    // particle world
	void setWorldLength(double length) {
		WorldLength=length;
		if (WorldLength<0) {
			cerr << "The x dimension of your world is less than zero. WorldLength has not been set"<< endl;		
			return;	}
		setXMax(WorldLength/2.0);
		setXMin(-WorldLength/2.0);
		setYMax(WorldLength/2.0);
		setYMin(-WorldLength/2.0);
		}
	void setVwind(double v){Vwind=v;}
	void setLandingElevation(double elev){LandingElevation=elev;}
	void setLGravity(Vec3D LG){LGravity=LG;}

    // atmospheric parameters
	void setTzeroC(double To){ TzeroC=To; setTzeroK(TzeroC);}
	void setLapse(double l){ Lapse=l; }
	void setRair(double R){ Rair=R; }
	
    // parameters not set by user
	void setN(double n){ N=n; }
	void setTzeroK(double TzeroC){ TzeroK=TzeroC+273.15;}


// Get Functions	

    // interaction parameters
    double getRest(){return Rest;}
	double getTc(){return Tc;}
	
	//particle parameters
	double getDensity(){return speciesHandler.getObject(0)->getDensity();}
	double getMaxParticleRadius(){return MaxParticleRadius;}
	double getMinParticleRadius(){return MinParticleRadius;}
	double getCdconst(){return Cdconst;}
	double getPSDn(){return PSDn;}
	
	double getEjectionMedVelocity(){ return EjectionMedVelocity; }
	double getEjectionVelocityVariance(){ return EjectionVelocityVariance; }
	double getEjectionVelocityMax(){ return EjectionVelocityMax; }
	double getEjectionVelocityMin(){ return EjectionVelocityMin; }
	double getEjectionVelocityFactor(){ return EjectionVelocityFactor;}
	double getEjectionMedAngle(){return EjectionMedAngle;}
	double getEjectionAngleVariance(){return EjectionAngleVariance;}
	double getEjectionRate(){return EjectionRate;}
	double getMinAngle(){return MinAngle;}

	double getBurstTime(){return BurstTime;}
	double getRestTime(){return RestTime;}
	int getnburstT(){return nburstT;}
	int getMaxNburst(){return MaxNburst;}
	
	int getSpeedMode(){return SpeedMode;}
	int getAngleMode(){return AngleMode;}
	int getVelocityMode(){return VelocityMode;}
	double getGasDensity() {return GasDensity;}
	double getGasVelocity() {return GasVelocity;}
	
    // vent dimensions
	double getVentDiameter(){return VentDiameter;}
	double getVentElevation(){return VentElevation;}
	Vec3D getConduitDirection(){return ConduitDirection;}
	double getConduitDepth(){return ConduitDepth;}
	
    //particle world
    double getWorldLength(){ return WorldLength; }
    double getVwind(){return Vwind;}
    double getLandingElevation(){return LandingElevation;}
    Vec3D getLGravity(){return LGravity;}
    	    
    // atmospheric parameters
	double getTzeroC(){return TzeroC;}
	double getLapse(){return Lapse;}
	double getRair(){return Rair;}
	    
    // parameters not set by user
	double getN(){return N;}
	double getTzeroK(){return TzeroK;}



// ***********************end set / get functions**************************//



// ******************** DEFINE PRIVATE VARIABLES / FUNCTIONS *****************//
private:

    // interaction parameters
	double Rest;
	double Tc;
	
	// particle parameters
	double MaxParticleRadius;
	double MinParticleRadius;
	double Cdconst;
	double PSDn;
		
	double EjectionMedVelocity;
	double EjectionVelocityVariance;
	double EjectionVelocityMax;
	double EjectionVelocityMin;
	double EjectionVelocityFactor;
	double EjectionMedAngle; 	
	double EjectionAngleVariance;
	double EjectionRate;
	double MinAngle;
	
	double BurstTime;
	double RestTime;
	int nburstT;
	int nburst;
	int MaxNburst;
	
	int SpeedMode;
	int AngleMode;
	int VelocityMode;
	double GasDensity;
	double GasVelocity;
		
	// vent dimenstions
	double VentDiameter; 
	double VentElevation;
	double ConduitDepth;		
	Vec3D ConduitDirection;	
	
	//particle world
	double WorldLength;
	double Vwind;
	double LandingElevation;
	Vec3D LGravity;
	
	//atmospheric parameters	
	double TzeroC;
	double Lapse;
	double Rair;

    //parameters not set by user
	double N;
	int count;
	int NumEjected;
	int NewSeed;
	double TzeroK;
	
	
	//File for saving x location of collision with ground
	File fpos;      
	ofstream verif_file;
};

// ***********************end private variable definitions********************//


//**************************end CinderDriver**********************************//






// ****************************************************************************//
// *********************** MAIN: SET USER PARAMETERS***************************//
// ****************************************************************************//

int main(int argc, char *argv[])
{
	CinderDriver problem;
	problem.setName("CinderColumn3D");
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    
	// set interaction parameters
	problem.setRest(0.5); 	   //0.5                        // restitution coefficient
	problem.setTc(2e-5); 				                    // contact time
	//problem.setTc(1e-5); 

	problem.restartFile.setFileType(FileType::ONE_FILE);	// turns on/off creation of restart file
	problem.fStatFile.setFileType(FileType::ONE_FILE);	    // turn on/off (force info) stress file
	problem.eneFile.setFileType(FileType::ONE_FILE);		// turn on/off (force info) stress file

	// particle parameters
	problem.setParticleRadius(.10,0.5);	                        // particle radius (min, max)
//	problem.setParticleRadius(.064,1.0);	                    // particle radius (min, max)
	problem.setDensity(1250); 		        	                // density (kg/m3)
	problem.setCdconst(1);		        	                    // drag coefficient
	problem.setPSDn(-3);                                        // particle size power distribution exponent
	
	
	// ejection parameters
	problem.setEjectionMedVelocity(30);		                    // median ejection velocity (m/s)
	problem.setEjectionVelocityVariance(20);	                // velocity variance
	problem.setEjectionVelocityMax(200);                        // maximum ejection velocity   <------USED?
	problem.setEjectionVelocityMin(10);                         // mininum ejection velocity   <------USED?
	problem.setEjectionVelocityFactor(0.75);                    // this % of particles will be calculated with gas relation, the rest randomly from a linear distribution  <----- NOT ENABLED
	problem.setEjectionMedAngle(80);			                // max ejection angle (degrees from horizontal)
	problem.setEjectionAngleVariance(30);		                // ejection angle variance
	problem.setEjectionRate(100);//3000		                // ejection rate (particles/s)
    problem.setMinAngle(45);


	problem.setBurstTime(0.1);				                    // number of seconds over which to eject particles. set to zero to erupt indefinitely
	problem.setRestTime(0.3);				                    // time between bursts of particles
	problem.setMaxNburst(1);                                    // total number of bursts
	
//	problem.setSpeedMode(1);                                    // Speed mode: 0=randomly selected, 1=coupled to gas jet
//	problem.setAngleMode(1);                                    // Angle mode: 0=randomly selected, 1=coupled to position
//    problem.setVelocityMode(1);                                 // 0=randomply selected, 1=spherical, coupled to gas jet
	problem.setGasDensity(286.98);                              // Density of gas jet gas
	problem.setGasVelocity(30);                                   // Gas jet velocity


	// vent dimensions
	problem.setVentDiameter(2.5);                               // vent diameter
	problem.setVentElevation(700.0);                            // vent elevation
//	problem.setConduitDirection(Vec3D(0,0,1));                  // conduit angle
	problem.setConduitDepth(2);					            // conduit depth: used to calculate acceleration at the vent


	// particle world
	problem.setWorldLength(100);	                            // length of landing surface (m)
	problem.setZMax(800);			                            // max height of display box
	problem.setZMin(700);			                            // min height of display box
	problem.setLandingElevation(690);
	problem.setGravity(Vec3D(0,0,0)); 	                    // gravity
	problem.setVwind(0);  
	problem.setLGravity(Vec3D(0,0,-9.81));                                      // wind velocity


	// atmospheric parameters for Earth
	problem.setTzeroC(25);			                            // temperature at sea level
	problem.setLapse(6.5);                                      // adiabatic lapse rate
	problem.setRair(286.98);                                    // atmospheric density


	// set normal spring coefficient, and contact dissipation
	double mass;
	mass=pow(problem.getMinParticleRadius(),3)*problem.getDensity()*4.0/3.0*constants::pi;          // mass of smallest particle, used to set collision time, restitution coefficient
	species->setCollisionTimeAndRestitutionCoefficient(problem.getTc(),problem.getRest(),mass);	    // set collision time and restitution coefficient
	problem.setTimeStep(problem.getTc()/50.0);		                                                // time step for calculation

	// time step for calculation
	species->setSlidingStiffness(species->getStiffness()*2.0/7.0);	// tangential spring coefficient	
	species->setSlidingDissipation(species->getDissipation());		// tangential dissipation
	species->setSlidingFrictionCoefficient(0.5);				    // friction coefficent tan fric_ang = mu(static friction=dynamic friction by default)

	//turns on rolling and torsional friction
	species->setTorsionStiffness(2./7.*species->getStiffness());    // spring coefficient
	species->setTorsionFrictionCoefficient(0.5);                    // friction coefficient
	species->setTorsionDissipation(0);                              // dissipation

	species->setRollingStiffness(2./7.*species->getStiffness());    //spring coefficient
	species->setRollingFrictionCoefficient(0.5);                    //friction coefficient
	species->setRollingDissipation(0.0);                            // dissipation


	// timing
	problem.setTimeMax(40);                         	            // max time to consider in calculation (0=infinity?)
	//problem.setSaveCount(problem.getTimeMax()/problem.getTimeStep()/9000.0);	        // total number of saves ENOUGH?
    problem.fStatFile.getFstream().precision(4);
    problem.fStatFile.setSaveCount(50000);  // 500 frames per second = 5000

/* ------------------------------------ Finished Setting Parameters ------------------------------------------------------*/

	// simulation
//	problem.setRNtypeLCG();
	problem.setXBallsAdditionalArguments(" -v0 -solid -s .01 -mo 2700 -moh 200");	// parameters for plotting	
	problem.autoNumber();			// auto number output files
	
	// write to screen
	cout << "Run number:" << problem.getRunNumber() << endl;
	//cout << "max speed allowed: " << problem.getMaximumVelocity() << endl;
	cout << "time step:" << problem.getTimeStep() << endl;
	cout << "mass: " << mass << endl;
	cout << "Tc: " << problem.getTc() << endl;
	cout << "K: " << species->getStiffness() << endl;
//	problem.readArguments(argc,argv);	// allows passing of variables at command line (w/out compiling)
	problem.solve();			// calls function to setup particles initial conditions

	
}
// ******************************end main***************************************//


