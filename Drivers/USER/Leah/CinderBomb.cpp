#include<iostream>

#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionIrreversibleAdhesiveSpecies.h"
#include <iomanip>

// NOTES: August 28, 2012
/* This code is driving me crazy. All I want is to randomly arrange particles in a box 1 m3. 
 These particles need to fall to the floor of the box and come to rest.
 Particles should fill the bottom 1/2 of the box after coming to rest.

 When both box and particles were smaller, this appeared to work.
 Currently particles fall through the bottom of the box, giving me a giant headache!!! */

class my_problem : public Mercury3D
{
public:

    my_problem()
    {
    	settled=false;
        flag = true;
    }

    void setupInitialConditions()
    {        
       //impactorSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionIrreversibleAdhesiveSpecies());
       //sandSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionIrreversibleAdhesiveSpecies());

     //   setParticleDimensions(3);                   // NECESARRY??? CORRECT??
      
        
        InfiniteWall w0;
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0,0,-getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0,0,getZMax()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(0,0,-getXMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(0,0,getXMax()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0,0,-getYMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0,0,getYMax()));
        wallHandler.copyAndAddObject(w0);

        int N1 = static_cast<int>(std::pow(N, 0.33)) + 1;

        BaseParticle P0;
        P0.setSpecies(sandSpecies);
        for (int i = 0; i < N; i++)
        {

            int ix = static_cast<int>(i % N1);
            int iz = static_cast<int>(i / N1 / N1);

            int iy = (i - ix - N1 * N1 * iz) / N1;

            double x = (getXMax() - getXMin()) * (ix + 1) / (N1 + 1);
            double y = (getYMax() - getYMin()) * (iy + 1) / (N1 + 1);
            double z = (getZMax() - getZMin()) * (iz + 1) / (N1 + 1);

            P0.setPosition(Vec3D(x, y, z));
            P0.setVelocity(Vec3D(random.getRandomNumber(-1.0, 1.0), random.getRandomNumber(-1.0, 1.0), random.getRandomNumber(-1.0, 1.0)));
            P0.setRadius(0.025);

            particleHandler.copyAndAddObject(P0);
        }

        double tc = 50e-5; //collision time
        double rc = 0.88; //restitution coefficient // (speed after collision)/(speed before) or sqrt(energy loss)
        sandSpecies->setCollisionTimeAndRestitutionCoefficient(tc, rc, P0.getRadius()*P0.getRadius()*P0.getRadius()* sandSpecies->getDensity() * 4.0/3.0 * constants::pi); //////// CORRECT?
        
        setTimeStep(tc/50.0);
    }

    int N;
    bool flag;
    
    void computeExternalForces(BaseParticle* P0){	
	// Call the MD compute_external_forces function (turns on gravity)
	DPMBase::computeExternalForces(P0);
	
	// delete above line to turn off mercury gravity

	// Drag Force - constant Cd
	double Po;              // part of pressure calculation
	double Temp;            // Temperature in the atmosphere
	double Pressure;        // Atmospheric pressure
	double density_air;     // Air density
	double Cd;              // Drag coefficient (constant)
	double A;               // Cross sectional area of particle
	double Fdragx;          // Force in the x dimension
	double Fdragy;          // Force in the y dimension
	double Fdragz;          // Force in the z dimension: note that this includes gravity
	double v;               // Particle velocity
	double mass;            // Particle mass

	A=3.1415*pow(P0->getRadius(),2);
	Cd=getCDConst();
	v=pow(pow((P0->getVelocity().X-getVWind()),2)+pow(P0->getVelocity().Y,2)+pow(P0->getVelocity().Z,2),0.5);	
	Po=101300.0*pow(((getTZeroK()-(getVentElevation())*getLapse()/1000.0)/getTZeroK()),(-9.81/(getRAir()*getLapse()/1000.0)));
	Temp=getTZeroK()-(P0->getPosition().Z+getVentElevation())*getLapse()/1000.0;
	Pressure=Po*pow(Temp/getTZeroK(),9.81/(getRAir()*getLapse()/1000.0));
	density_air=Pressure/(getRAir()*Temp);
	mass=getDensity()*(4.0/3.0*3.14159*pow(P0->getRadius(),3)); 
    
    }
    
    

// Set Functions
    void set_Depth(double depth)    {Depth = depth;    }
    void setVentElevation(double vent) {VentElevation=vent;}
    void setVWind(double wind){VWind=wind;}
    void setLapse(double L){Lapse=L;}
    void setTZeroK(double T){TZeroK=T;}
    void setRAir(double r){RAir=r;}
    void setRadius(double r){Radius=r;}
    void setDensity(double d){Density=d;}
    void setMass(double m){Mass=m;}
    void setCDConst(double cd){CDConst=cd;}
    void setRImpact(double r){RImpact=r;}
    void set_Length(double length)
    {
        Length = length;
//		setXMin(-length/2.0);
//		setXMax(length);
//		setYMin(-length/2.0);
//		setYMax(length);
//		setZMax(length);
    }

// Get Functions
    double get_Depth()    {return Depth;    }
    double get_Length()    {return Length;    }
    double getVentElevation(){return VentElevation;}
    double getVWind(){return VWind;}
    double getLapse(){return Lapse;}
    double getTZeroK(){return TZeroK;}
    double getRAir(){return RAir;}
    double getRadius(){return Radius;}
    double getDensity(){return Density;}
    double getCDConst(){return CDConst;}
    double getMass(){return Mass;}
    double getRImpact(){return RImpact;}

   
private:
    //void computeExternalForces(BaseParticle *P){computeWalls(P);}
    double Length;
    double Depth;
    bool settled;
    
    double VentElevation;
    double VWind;
    double Lapse;
    double TZeroK;
    double RAir;
    double Radius;
    double Density;
    double Mass;
    double CDConst;
    double RImpact;



	void actions_before_time_step(){
		double t=getTime();

		//insert the impacting particle
		if ((t>10) && (flag))
			{
			// delete particle above z=0.3 m
		//	for (unsigned int i=0;i<getParticleHandler().getNumberOfParticles();){
			//if (getParticleHandler().getParticle(i)->getPosition().Z>0.4)
				//removeParticle(i);
				//else i++;}

			//insert particle
			double x=0.3;	//location where impactor appears
			double y=0.3;
			double z=0.9;
			
		    BaseParticle P0;
            P0.setSpecies(sandSpecies);
			P0.setPosition(Vec3D(x,y,z));
			P0.setVelocity(Vec3D(0.0,0.0,-200.0));
			P0.setRadius(getRImpact());

		
			// Insert particle in grid
			P0.setHandler(&particleHandler);
			particleHandler.copyAndAddObject(P0);
			std::cout << "New particle created " << std::endl;

			flag=false;	//flag used so that only one particle is created

			} 

	}

    // Print time to screen
    void printTime() const
    {
        std::cout << "\rt=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
            << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax() << "\r";
        std::cout.flush();
    }
 
public:
     LinearViscoelasticFrictionIrreversibleAdhesiveSpecies* sandSpecies;  
     LinearViscoelasticFrictionIrreversibleAdhesiveSpecies* impactorSpecies;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc UNUSED, char *argv[] UNUSED)
{

    ///Start off my solving the default problem
    my_problem problem;

    //problem.N=1200;
    problem.N = 500;
    problem.setName("CinderBomb");
    //	problem.set_dissipation(100.0);
    
     //problem.impactorSpecies = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionIrreversibleAdhesiveSpecies());
     //problem.sandSpecies = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionIrreversibleAdhesiveSpecies());
     LinearViscoelasticFrictionIrreversibleAdhesiveSpecies species;
     problem.impactorSpecies = problem.speciesHandler.copyAndAddObject(species);
     problem.sandSpecies = problem.speciesHandler.copyAndAddObject(species);
     
     auto mixedSpecies = problem.speciesHandler.getMixedObject(problem.sandSpecies, problem.impactorSpecies);
//     auto mixedSpecies = dynamic_cast<LinearViscoelasticFrictionIrreversibleAdhesiveSpecies*>(problem.speciesHandler.getMixedObject(0,1));
    
    problem.setVentElevation(698);
    problem.setVWind(0);
    problem.setLapse(2);
    problem.setTZeroK(270);
    problem.setRAir(160);
    problem.setRadius(0.6);
    problem.setDensity(1250);
    problem.setMass(1);
    
    
    problem.setGravity(Vec3D(0.0, 0.0, -100.0));
    problem.sandSpecies->setStiffness(1e3);
    problem.sandSpecies->setAdhesionStiffness(1e3);
    problem.sandSpecies->setAdhesionForceMax(1e3);
    problem.impactorSpecies->setStiffness(1e3);
    problem.impactorSpecies->setAdhesionStiffness(1e3);
    problem.impactorSpecies->setAdhesionForceMax(1e3);
    problem.setSaveCount(1000); //set_number_of_saves_all(1000);
    problem.setTimeMax(1.5);
    problem.setSystemDimensions(3);
    problem.set_Length(1);
    problem.autoNumber();   
    
    problem.impactorSpecies->setDensity(1000.0);
    problem.sandSpecies->setDensity(1000.0);
    
    //Turning On Friction: sand
    problem.sandSpecies->setSlidingStiffness(problem.sandSpecies->getStiffness()*2.0/7.0);	// tangential spring coefficient	
	problem.sandSpecies->setSlidingDissipation(problem.sandSpecies->getDissipation());		// tangential dissipation
	problem.sandSpecies->setSlidingFrictionCoefficient(0.5);				    // friction coefficent tan fric_ang = mu(static friction=dynamic friction by default)

	problem.sandSpecies->setTorsionStiffness(2./7.*problem.sandSpecies->getStiffness());        // spring coefficient
	problem.sandSpecies->setTorsionFrictionCoefficient(0.5);                    // friction coefficient
	problem.sandSpecies->setTorsionDissipation(0);                              // dissipation

	problem.sandSpecies->setRollingStiffness(2./7.*problem.sandSpecies->getStiffness());    //spring coefficient
	problem.sandSpecies->setRollingFrictionCoefficient(0.5);                    //friction coefficient
	problem.sandSpecies->setRollingDissipation(0.0);                            // dissipation
	
	    //Turning On Friction: impactor
    problem.impactorSpecies->setSlidingStiffness(problem.impactorSpecies->getStiffness()*2.0/7.0);	// tangential spring coefficient	
	problem.impactorSpecies->setSlidingDissipation(problem.impactorSpecies->getDissipation());		// tangential dissipation
	problem.impactorSpecies->setSlidingFrictionCoefficient(0.5);				    // friction coefficent tan fric_ang = mu(static friction=dynamic friction by default)

	problem.impactorSpecies->setTorsionStiffness(2./7.*problem.impactorSpecies->getStiffness());        // spring coefficient
	problem.impactorSpecies->setTorsionFrictionCoefficient(0.5);                    // friction coefficient
	problem.impactorSpecies->setTorsionDissipation(0);                              // dissipation

	problem.impactorSpecies->setRollingStiffness(2./7.*problem.impactorSpecies->getStiffness());    //spring coefficient
	problem.impactorSpecies->setRollingFrictionCoefficient(0.5);                    //friction coefficient
	problem.impactorSpecies->setRollingDissipation(0.0);                            // dissipation

    //Turning On Friction: mixed
    mixedSpecies->setSlidingStiffness(mixedSpecies->getStiffness()*2.0/7.0);	// tangential spring coefficient	
	mixedSpecies->setSlidingDissipation(mixedSpecies->getDissipation());		// tangential dissipation
	mixedSpecies->setSlidingFrictionCoefficient(0.5);				    // friction coefficent tan fric_ang = mu(static friction=dynamic friction by default)

	mixedSpecies->setTorsionStiffness(2./7.*mixedSpecies->getStiffness());        // spring coefficient
	mixedSpecies->setTorsionFrictionCoefficient(0.5);                    // friction coefficient
	mixedSpecies->setTorsionDissipation(0);                              // dissipation

	mixedSpecies->setRollingStiffness(2./7.*mixedSpecies->getStiffness());    //spring coefficient
	mixedSpecies->setRollingFrictionCoefficient(0.5);                    //friction coefficient
	mixedSpecies->setRollingDissipation(0.0);                            // dissipation

    problem.setZMin(0);
    problem.setZMax(1);
    problem.setYMin(0);
    problem.setYMax(1);
    problem.setXMin(0);
    problem.setXMax(1);
    
    	//Notifications
//	std::cout << "number of particles in bed: " << problem.N << std::endl;
//	std::cout << "radius of particles in bed: " << problem.get_Rbed() << std::endl;	
//	std::cout << "collision time: " << problem.get_tc() << std::endl;
//	std::cout << "stiffness: " << problem.get_k() << std::endl;
//	std::cout << "particle massI: " << problem.get_MassI() << std::endl;

    std::cout << "Solving by checking all collisions - 3D version" << std::endl;
    problem.solve();

}

