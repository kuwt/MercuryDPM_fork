#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Species/LinearViscoelasticFrictionSpecies.h>

#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
// NOTE: r is restitution coefficient not radius!! Use getMinParticleRadius when setting collision time and stiffness

class CinderDriver : public Mercury3D
{
    
public:

// ******************** SET DEFAULTS, overwritten by set function (invisible to user)************//
    CinderDriver()
    {
        species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

//	EjectionMaxAngle=0.1;	
        ConduitAMag = 1;
        ConduitDirection = Vec3D(0, 0, 1);
        flag = true;
        verif_file.open("verification_file.dat");
        NumEjected = 0;
    }
//********************************end default values***************************//
    
// Use to change contact law. NOT RECOMENDED!!!
// Might want to run once with a different contact model, almost certainly will do my runs with current contact model (spring-dashpot)
//void computeInternalForces(BaseParticle *PI, Particle *PJ) {
//		computePlasticInternalForces(PI, PJ);	// Turn on plastic contact law. under redevelopment. 
    // require different parameters than restitution coefficient, etc.
//		visco-elastic and sintering contact laws are under development. 
//	}
    
///Defines the smallest Particle to be the smallest possible Ejected particle.
    BaseParticle* getSmallestParticle()
    {
        
        //P0 is the smallest possible ejected particle
        SphericalParticle P0;
        P0.setSpecies(species);
        P0.setRadius(MinParticleRadius);
        BaseParticle* Pp = &P0;
        return Pp;
        
//	for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) { 
//		if (particleHandler.getObject(i)->getRadius()<pP->getRadius()) pP = particleHandler.getObject(i); 
//	}
//	return pP;
    }
    
///////////////////////////////////////
    void write(std::ostream& os, bool print_all) const override
            {
        DPMBase::write(os, print_all);
        os << "-----------------------------" << std::endl;
        os << "Problem Name: " << dataFile.getName() << std::endl; //problem_name.str().c_str() << endl;
        os << " " << std::endl;
        
        os << "Time Step: " << getTimeStep() << std::endl;
        os << "End time: " << getTimeMax() << std::endl;
        os << "Dimensions of problem:" << getSystemDimensions() << std::endl;
        os << "  x :[" << getXMin() << "," << getXMax() << "]    y:[" << getYMin() << "," << getYMax() << "]    z:[" << getZMin() << "," << getZMax() << "]" << std::endl;
        os << " " << std::endl;
        
        //os << " Number of particles ejected:" << get_N() << ", Memory alloted for a maximum of:" << get_Nmax() << " particles" << endl;
        //if (print_all) for (unsigned int i=0; i<Particles.size(); i++) { os << "  "; Particles[i].write(os); os << endl; }
        os << "Walls: " << std::endl;
        os << "   Number of walls:" << wallHandler.getNumberOfObjects() << std::endl; //get_NWall() << endl;
        os << "   Number of Periodic walls:" << boundaryHandler.getNumberOfObjects() << std::endl; //get_NWallPeriodic() << endl;
        /*for (int i = 0; i < wallHandler.getNumberOfObjects(); i++)
         {
         os << "     ";
         wallHandler.getObject(i)->write(os);
         os << std::endl;
         }*/
        for (int i = 0; i < boundaryHandler.getNumberOfObjects(); i++)
        {
            os << "     ";
            boundaryHandler.getObject(i)->write(os);
            os << std::endl;
        }
        //	os << "   Number of walls:" << get_NWall() << endl;
        //	os << "   Number of Periodic walls:" << get_NWallPeriodic() << endl;
        //	for (int i=0; i<get_NWall(); i++) { os << "     "; Walls[i].write(os); os << endl; }
        //	for (int i=0; i<get_NWallPeriodic(); i++) { os << "     "; WallsPeriodic[i].write(os,print_all); os << endl; }
        os << " " << std::endl;
        os << " " << std::endl;
        
        os << "Interaction Properties:" << std::endl;
        os << "   stiffness (k): " << species->getStiffness() << std::endl;
        os << "   normal dissipation (disp): " << species->getDissipation() << std::endl;
        os << "   tangential stiffness: " << species->getSlidingStiffness() << std::endl;
        os << "   tangential dissipation: " << species->getSlidingDissipation() << std::endl;
	    os << "   dynamic coulomb friction coefficient (mu): " << species->getSlidingFrictionCoefficient() << "    static coulomb friction coefficient (mus): " << species->getSlidingFrictionCoefficientStatic() << std::endl;
        //os << "   cohesion: " << species->getCohesionStiffness() << std::endl;
        os << " " << std::endl;
        //k2max:0
        
        os << " Gravity:" << getGravity() << std::endl;
        os << " Density of particles: " << species->getDensity() << std::endl;
        //os << " Maximum Ejection Angle (deg):" << get_EjectionMaxAngle() << endl;
        os << " Vent Diameter:" << get_VentDiameter() << std::endl;
        os << " Conduit Acceleration:" << get_ConduitAcceleration() << std::endl;
        os << " " << std::endl;
        os << " " << std::endl;
        
        os << "File Saving: " << std::endl;
        os << "   force statistics: " << static_cast<unsigned int>(fStatFile.getFileType()) //get_options_fstat()
                << "   data: " << static_cast<unsigned int>(dataFile.getFileType()) //get_options_data()
                << "   energy: " << static_cast<unsigned int>(eneFile.getFileType()) //get_options_ene()
                << "   restart: " << static_cast<unsigned int>(restartFile.getFileType()) << std::endl; //get_options_restart() << endl;
        os << "Frequency of Saves: " << std::endl;
        os << "   data: " << dataFile.getSaveCount() << std::endl;
        os << "   energy statistics: " << eneFile.getSaveCount() << std::endl;
        os << "   force statistics: " << fStatFile.getSaveCount() << std::endl;
        os << " " << std::endl;
    }
    
// ************************* SET PROBLEM GEOMETRY *****************************//
    
    void setupInitialConditions() override
    {
        
        setYMax(MaxParticleRadius); // y dimension restricted to width of particle
        setYMin(-MaxParticleRadius);
        
        PeriodicBoundary b0;
        b0.set(Vec3D(0., 1., 0.), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        //set_NWallPeriodic(1);				// periodic walls bounding y dqimension (no force)
        //WallsPeriodic[0].set(Vec3D(0,1,0), getYMin(), getYMax());
        //set_NWall(0);
        
        // infinitely deep bottom wall
        //wallHandler.getObject(0)->set(Vec3D(0.0, 0.0, -1.0),-getZMin());//+getMaxInflowParticleRadius()*2.0);	
        
        //create base
        int N_quarter = (get_WorldLength() - get_VentDiameter()) / (2.0 * getFixedParticleRadius()) / 2;
        set_N(N_quarter * 4); // number of fixed particles on bottom surface
        for (int i = 0; i < N_quarter; i++)
        {
            
            //top left layer
            SphericalParticle p1;
            p1.setSpecies(species);
            p1.setPosition(Vec3D(getXMin() + (2.0 * i + 1.0) * getFixedParticleRadius(), 1, 0.0));
            p1.setRadius(getFixedParticleRadius());
            p1.fixParticle();
            particleHandler.copyAndAddObject(p1);
            
            //bottom left layer
            SphericalParticle p2;
            p2.setSpecies(species);
            p2.setPosition(Vec3D(getXMin() + (2.0 * i + 1.0) * getFixedParticleRadius(), 1, -getFixedParticleRadius() * 2.0));
            p2.setRadius(getFixedParticleRadius());
            p2.fixParticle();
            particleHandler.copyAndAddObject(p2);
            
            //top right layer
            SphericalParticle p3;
            p3.setSpecies(species);
            p3.setPosition(Vec3D(getXMax() - (2.0 * i + 1.0) * getFixedParticleRadius(), 1, 0.0));
            p3.setRadius(getFixedParticleRadius());
            p3.fixParticle();
            particleHandler.copyAndAddObject(p3);
            
            //bottom right layer
            SphericalParticle p4;
            p4.setSpecies(species);
            p4.setPosition(Vec3D(getXMax() - (2.0 * i + 1.0) * getFixedParticleRadius(), 1, -getFixedParticleRadius() * 2.0));
            p4.setRadius(getFixedParticleRadius());
            p4.fixParticle();
            particleHandler.copyAndAddObject(p4);
            
        }
        
        writeRestartFile(); // saves input variables to file CinderDriver.X.inf
        
    }
    
// ***************************end problem geometry***************************//
    
// *************************** SET PARTICLE PARAMETERS ***********************//
    
    void actionsBeforeTimeStep() override
    {
        //if (flag)		// hack to do one particle
        
        double t, CurrentRate;
        t = getTime();
        CurrentRate = NumEjected / t;
        
        if (CurrentRate < EjectionRate)
        {
            SphericalParticle P0;
            P0.setSpecies(species);
            P0.setRadius(random.getRandomNumber(MinParticleRadius, MaxParticleRadius)); // particle radius
//			P0.computeMass(Species);	// particle mass  (Species)?????
//			can ask anthony for another way to get mass if using particles of different densities
                    
            // particles form in middle of calculation box, at the bottom
            P0.setPosition(Vec3D(0.5 * get_VentDiameter() * random.getRandomNumber(-1, 1), 0.0, getZMin() + P0.getRadius()));
            double Angle = (get_EjectionMedAngle() + get_EjectionAngleVariance() * random.getRandomNumber(-1, 1)) * constants::pi / 180; // ejection angle
            double Speed = get_EjectionVelocity() + get_EjectionVelocityVariance() * random.getRandomNumber(-1, 1); // ejection speed
                    
            //	double Angle = 25*constants::pi/180;	// comment out to return control to user
            //	double Speed = 20;	
                    
            P0.setVelocity(Vec3D(Speed * cos(Angle), 0.0, Speed * sin(Angle)));
            
            if (checkParticleForInteraction(P0))
            {
                NumEjected++;
                // Insert particle in grid
                particleHandler.copyAndAddObject(P0);
            }
            
        }
    }
    
// ****************************end particle parammeters ***********************//	
    
    // TALK ABOUT THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
//	void write(std::ostream& os, bool print_all);	// declare print function
//	leads to errors!!
    
    void actionsAfterTimeStep() override
    {
        // Writes out certain particle info to screen 
        //		double z=particleHandler.getObject(996)->getPosition().Z;
//		//cout << z << endl;
//
//		//cout << ((z<0) && (flag));
//
//		if (z<0 && flag)
//			{
//			double x=particleHandler.getObject(996)->getPosition().X;
//			verif_file << x <<endl;
//			flag=false;
//			}
    }
    
//****************** EXTERNAL FORCES i.e. drags, gravity etc...***************//
    
    void computeExternalForces(BaseParticle* P0) override

    {
        // Call the MD computeExternalForces function (turns on gravity)
        DPMBase::computeExternalForces(P0);

        //double g = getGravity().getLength();
        
        // Drag Force
        /*

         double Po;
         double Temp;
         double Pressure;
         double density_air;
         double Cd;
         double A;
         double Fdmag;
         double Fdragx;
         double Fdragz;
         double v2;
         double theta;
         
         A=3.1415*pow(P0->getRadius(),2);
         Cd=get_Cdconst();
         
         v2=(pow(P0->getVelocity().X,2)+pow(P0->getVelocity().Z,2));	//velocity squared
         theta=atan(P0->getPosition().X/P0->getPosition().Z);
         
         Po=101300.0*pow(((get_TzeroK()-get_ventElevation()*get_lapse()/1000)/get_TzeroK()),(-9.81/(get_rair()*get_lapse()/1000)));	
         Temp=get_TzeroK()-(P0->getPosition().Z+get_ventElevation())*get_lapse()/1000;

         Pressure=Po*pow(Temp/get_TzeroK(),getGravity().Z/(get_rair()*get_lapse()/1000));
         density_air=Pressure/(get_rair()*Temp);
         Fdmag=1/2*density_air*v2*Cd*A;
         Fdragx=Fdmag*sin(theta);
         Fdragz=Fdmag*cos(theta);

         
         P0->addForce(Vec3D(Fdragx,0.0,Fdragz));

         */

    }
//**************************end external forces*******************************//
    
// ********************** VARIABLE FUNCITONS ****************//
    
// Set Functions
    void set_EjectionMedAngle(double new_vel_ang)
    {
        EjectionMedAngle = new_vel_ang;
    }
    void set_EjectionAngleVariance(double ang_var)
    {
        EjectionAngleVariance = ang_var;
    }
    void set_VentDiameter(double new_dim)
    {
        VentDiameter = new_dim;
    }
    void set_ConduitDirection(Vec3D dir)
    {
        ConduitDirection = dir;
        set_ConduitAcceleration(ConduitAMag, ConduitDirection);
    }
    void set_ConduitAcceleration(double mag, Vec3D dir)
    { // normalized acceleration direction
        ConduitAcceleration = 1 / (pow((pow(dir.X, 2) + pow(dir.Y, 2) + pow(dir.Z, 2)), (1 / 3))) * mag * dir;
    }
    void set_ConduitAMag(double mag)
    {
        ConduitAMag = mag;
        set_ConduitAcceleration(ConduitAMag, ConduitDirection);
    }
    void set_N(double n)
    {
        N = n;
    }
    void set_TzeroC(double To)
    {
        TzeroC = To;
        set_TzeroK(TzeroC);
    }
    void set_lapse(double l)
    {
        lapse = l;
    }
    void set_rair(double R)
    {
        rair = R;
    }
    void set_TzeroK(double TzeroC)
    {
        TzeroK = TzeroC + 273.15;
    }
    void set_ventElevation(double event)
    {
        ventElevation = event;
    }
    void set_Cdconst(double cconst)
    {
        Cdconst = cconst;
    }
    
    void set_ParticleRadius(double rmin, double rmax)
    {
        if (rmin > rmax)
        {
            std::cerr << "Your maximum radius is less than your minimum. Particle radius has not been set" << std::endl;
            return;
        }
        if (rmin < 0)
        {
            std::cerr << "Your minimim radius is less than 0, radius has not been set" << std::endl;
            return;
        }
        if (rmax < 0)
        {
            std::cerr << "Your minimim radius is less than 0, radius has not been set" << std::endl;
            return;
        }
        MaxParticleRadius = rmax;
        MinParticleRadius = rmin;
    }
    
    void set_WorldLength(double length)
    {
        WorldLength = length;
        if (WorldLength < 0)
        {
            std::cerr << "The x dimension of your world is less than zero. WorldLength has not been set" << std::endl;
            return;
        }
        setXMax(WorldLength / 2.0);
        setXMin(-WorldLength / 2.0);
    }
    
    void setFixedParticleRadius(double rad)
    {
        FixedParticleRadius = rad;
    }
    void set_EjectionVelocity(double vel)
    {
        EjectionVelocity = vel;
    }
    void set_EjectionVelocityVariance(double vel)
    {
        EjectionVelocityVariance = vel;
    }
    void set_EjectionRate(double rate)
    {
        EjectionRate = rate;
    }
    
// Get Functions	
    double get_EjectionMedAngle()
    {
        return EjectionMedAngle;
    }
    double get_EjectionAngleVariance()
    {
        return EjectionAngleVariance;
    }
    double get_VentDiameter() const
    {
        return VentDiameter;
    }
    double get_ConduitAMag()
    {
        return ConduitAMag;
    }
    Vec3D get_ConduitDirection()
    {
        return ConduitDirection;
    }
    Vec3D get_ConduitAcceleration() const
    {
        return ConduitAcceleration;
    }
    double get_N()
    {
        return N;
    }
    double get_TzeroC()
    {
        return TzeroC;
    }
    double get_lapse()
    {
        return lapse;
    }
    double get_rair()
    {
        return rair;
    }
    double get_TzeroK()
    {
        return TzeroK;
    }
    double get_ventElevation()
    {
        return ventElevation;
    }
    double get_Cdconst()
    {
        return Cdconst;
    }
    double get_MaxParticleRadius()
    {
        return MaxParticleRadius;
    }
    double get_MinParticleRadius()
    {
        return MinParticleRadius;
    }
    double get_WorldLength()
    {
        return WorldLength;
    }
    double getFixedParticleRadius()
    {
        return FixedParticleRadius;
    }
    double get_EjectionVelocity()
    {
        return EjectionVelocity;
    }
    double get_EjectionVelocityVariance()
    {
        return EjectionVelocityVariance;
    }
    double get_EjectionRate()
    {
        return EjectionRate;
    }
    
// ***********************end set / get functions**************************//
    
    void actionsAfterSolve() override // After model completes running
    {
        verif_file.close();
    }
    
// ******************** DEFINE PRIVATE VARIABLES / FUNCTIONS *****************//
private:
    double EjectionMedVelocity;

    double EjectionAngleVariance;
    double EjectionVelocity;
    double EjectionVelocityVariance;
    double EjectionMedAngle; // in deg
    double VentDiameter; // meters
    double ConduitAMag;
    double N;
    double EjectionRate;
    Vec3D ConduitDirection;
    Vec3D ConduitAcceleration; // acceleration in the conduit
    
    double TzeroC;
    double lapse;
    double rair;
    double TzeroK;
    double ventElevation;
    double Cdconst;

    double MaxParticleRadius;
    double MinParticleRadius;
    double WorldLength;
    double FixedParticleRadius;

    int NumEjected;

    //flag if particle has gone below z<0
    bool flag;
    //File for saving x location of collision with ground.
    std::ofstream verif_file;
public:
    LinearViscoelasticFrictionSpecies* species;
};
// ***********************end private variable definitions********************//

//**************************end CinderDriver**********************************//

// ****************************************************************************//
// *********************** MAIN: SET USER PARAMETERS***************************//
// ****************************************************************************//

int main(int argc, char *argv[])
{
    CinderDriver problem;
    problem.setName("CinderDriver");
    
    // set interaction parameters
    double tc = 1e-5; //collision time
    double r = 0.1; //restitution coefficient	// (speed after collision)/(speed before) or sqrt(energy loss)
    double density = 1000;
    
    problem.restartFile.setFileType(FileType::ONE_FILE); //	problem.restartFile.setFileType(FileType::ONE_FILE);	// turns on/off creation of restart file
    problem.fStatFile.setFileType(FileType::ONE_FILE); //	problem.fStatFile.setFileType(FileType::ONE_FILE);	// turn on/off (force info) stress file
    problem.eneFile.setFileType(FileType::ONE_FILE); //	problem.eneFile.setFileType(FileType::ONE_FILE);	// turn on/off (force info) stress file
            
    // particle parameters
    problem.set_ParticleRadius(0.1, 0.1); // particle radius (min, max)
    problem.species->setDensity(density); // density
    problem.set_Cdconst(1);
    
//	cout << "angle " << problem.get_Angle() << endl;
//	cout << "velocity " << problem.getVelocity() << endl;
    
    // ejection parameters
    problem.set_EjectionVelocity(25); // velocity at bootom of conduit
    problem.set_EjectionVelocityVariance(1); // velocity variance
    problem.set_EjectionMedAngle(60); // max ejection angle (from vertical)
    problem.set_EjectionAngleVariance(20);
    problem.set_EjectionRate(1000);
//	problem.set_ConduitAMag(1125);				// acceleration within the conduit
    
    // conduit dimensions
    problem.set_ConduitDirection(Vec3D(0, 0, 1));
    problem.set_VentDiameter(10.0);
    problem.setFixedParticleRadius(problem.get_MaxParticleRadius()); // radius of particles fixed to ground
    problem.set_ventElevation(0);
    
    // particle world
    problem.set_WorldLength(100); // length of landing surface
    problem.setZMax(50); // dimensions of the problem
    problem.setZMin(0);
    problem.setGravity(Vec3D(0, 0, -9.8)); // gravity
            
    // Atmospheric parameters
    problem.set_TzeroC(25);
    problem.set_lapse(6.5);
    problem.set_rair(286.98);
    
    // sets normal spring coefficient, and contact dissipation
    problem.species->setCollisionTimeAndRestitutionCoefficient(tc, r, r * r * r * density * 4.0 / 3.0 * constants::pi);
    problem.setTimeStep(tc / 50.); // time step for calculation
            
    // time step for calculation
    problem.species->setSlidingStiffness(problem.species->getStiffness() * 2.0 / 7.0); // tangential spring coefficient
    problem.species->setSlidingDissipation(problem.species->getDissipation()); // tangential dissipation
    problem.species->setSlidingFrictionCoefficient(0.5); // friction coefficent tan fric_ang = mu(static friction=dynamic friction by default)

    //turns on rolling and torsional friction
    problem.species->setTorsionStiffness(2. / 7. * problem.species->getStiffness());
    problem.species->setTorsionFrictionCoefficient(10);
    problem.species->setTorsionDissipation(0);
    
    problem.species->setRollingStiffness(2. / 7. * problem.species->getStiffness());
    problem.species->setRollingFrictionCoefficient(10);
    problem.species->setRollingDissipation(0.0);
    
    problem.setTimeMax(100.0); // max time problem will run (=infinity)
            
    //std::cout << "max speed allowed: " << problem.particleHandler.getSmallestParticle()->calculateMaximumVelocity(problem.speciesHandler.getObject()) << std::endl;
    problem.setSaveCount(10000);
    
    // simulation
    problem.setXBallsAdditionalArguments(" -v0 -solid "); // parameters for plotting	
    problem.autoNumber(); // auto number output files
//	problem.readArguments(argc,argv);	// allows passing of variables at command line (w/out compiling)
    problem.solve(); // calls function to setup particles initial conditions
    
}
// ******************************end main***************************************//

