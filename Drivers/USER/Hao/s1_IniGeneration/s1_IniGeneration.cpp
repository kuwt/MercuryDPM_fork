//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
//#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class initial_poly : public Mercury3D{
public:
/*! s1_IniGeneration creates the initial configuration before the 
 * isotropic compresssion process. The particles, with given 
 * volume fracion and polydispersity, are placed in a cubic box 
 * of dimensions Xmax, Ymax, Zmax and period boundaries are applied
 * in all directions.
 * The particles are elastoplastic (en = 0.804) and frictionless.
 * Due to the overlap, particles bounce in all directions.
 * The relax stage after the generation of this configuration will
 * remove the overlap between the particles to get it homogeneous.
 * In 's2_IsotropicCompression', the last particle positions configuration (restart)
 * will be loaded and then particles will be compressed isotropically.
 * */
 
	initial_poly()
	{
		setSystemDimensions(3);
        particleSpecies = speciesHandler.copyAndAddObject(LinearPlasticViscoelasticFrictionSpecies());
        
	}
	
	LinearPlasticViscoelasticFrictionSpecies* particleSpecies;
    Mdouble particleDiameter, rhop, en, K1, K2, Kc, Phic, mu_slid, mu_roll,mu_tor, poly;
    Mdouble volumeFraction, N, tmax, H;
    Mdouble dampingCoefficient = 0.0;
    
    void computeExternalForces (BaseParticle * CI) override	
	{
		DPMBase::computeExternalForces (CI);
		if (!CI->isFixed())
        {
          // Background damping
          CI->addForce(-dampingCoefficient * CI->getVelocity());
        }
	}
	
    void setupInitialConditions()
    {
		
		double mass = rhop*constants::pi*mathsFunc::cubic(2*particleDiameter/(poly+1))/6.0;  //mass or effective mass in Mercury =/frac{2}{1/m1+1/m2}
		double tc = std::sqrt(mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) )); //contact duration for smallest pair of particles
		double Lm = std::pow(N*constants::pi*mathsFunc::cubic(particleDiameter)/6.0/volumeFraction,1.0/3.0)*1.000;   //estimated box length
		Mdouble velocity = 0.0;
		
		setTimeStep(tc/50);
        setTimeMax(tmax);
        
		 // particleSpecies
		particleSpecies->setDensity(rhop);
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(tc,en,mass);
        particleSpecies->setPlasticParameters(K1,K2,Kc,Phic);
		if ( mu_slid == 0){
		particleSpecies->setSlidingStiffness(0.0);
		}
		else{
		particleSpecies->setSlidingStiffness(2.0/10.0*particleSpecies->getLoadingStiffness());
		}
		if ( mu_roll == 0){
		particleSpecies->setRollingStiffness(0.0);
		}
		else{
		particleSpecies->setRollingStiffness(2.0/10.0*particleSpecies->getLoadingStiffness());
		}
		if ( mu_tor == 0){
		particleSpecies->setTorsionStiffness(0.0);
		}
		else{
		particleSpecies->setTorsionStiffness(2.0/10.0*particleSpecies->getLoadingStiffness());
		}
		particleSpecies->setSlidingFrictionCoefficient(mu_slid);
		particleSpecies->setSlidingFrictionCoefficientStatic(mu_slid);
		particleSpecies->setRollingFrictionCoefficient(mu_roll);
		particleSpecies->setRollingFrictionCoefficientStatic(mu_roll);
		particleSpecies->setTorsionFrictionCoefficient(mu_tor);
		particleSpecies->setTorsionFrictionCoefficientStatic(mu_tor);
		particleSpecies->setSlidingDissipation(2.0/10.0*particleSpecies->getDissipation());
		particleSpecies->setRollingDissipation(2.0/10.0*particleSpecies->getDissipation());
		particleSpecies->setTorsionDissipation(2.0/10.0*particleSpecies->getDissipation());
		dampingCoefficient =  0.6*particleSpecies->getDissipation();
        setGravity(Vec3D(0.0,0.0,0.0));
		
		BaseParticle p0;
        p0.setSpecies(particleSpecies);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        Mdouble Vp = 0.0;
         
         for (int i=0; i < N; i++) {
			p0.setRadius(random.getRandomNumber(particleDiameter/(poly+1),particleDiameter*poly/(poly+1)));
			double posX = random.getRandomNumber(0,Lm);
			double posY = random.getRandomNumber(0,Lm);
			double posZ = random.getRandomNumber(0,Lm);
			p0.setPosition(Vec3D(posX,posY,posZ));
			particleHandler.copyAndAddObject(p0);
			Vp = Vp + p0.getVolume();
		}
		double L = std::pow(Vp/volumeFraction,1.0/3.0);
        setXMax(L);
        setYMax(L);
        setZMax(L);
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);

        
        // Lees Edwards bc in y direction & periodic boundary in x direction
        //LeesEdwardsBoundary leesEdwardsBoundary;
        //leesEdwardsBoundary.set(
            //[velocity] (double time) { return time*velocity; },
            //[velocity] (double time UNUSED) { return velocity; },
            //getXMin(),getXMax(),getYMin(),getYMax());
        //boundaryHandler.copyAndAddObject(leesEdwardsBoundary);
        
        // periodic boundary 
		PeriodicBoundary normWall;
        normWall.set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(normWall);
        

        
		std::cout << "N = " << N << std::endl;
        std::cout << "nu = " << volumeFraction << std::endl;
        
        std::cout << "Vp = " << Vp << std::endl;
        std::cout << "Vt = " << L*L*L << std::endl;
        std::cout << "Initial box L = " << Lm << std::endl;
        std::cout << "Corrected box Lx = " << (getXMax()-getXMin()) << ", Ly = " << (getYMax()-getYMin()) << ", Lz = " << (getZMax()-getZMin()) << std::endl;
        std::cout << "output = " << getName() << std::endl;
        std::cout << "mass = " << mass << std::endl;
        std::cout << "tc = " << tc << std::endl;
        std::cout << "delta t = " << getTimeStep() << "," << particleSpecies->computeTimeStep(mass) << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;

        
    }
    
	
    
    
    
};

int main(int argc , char *argv[] )
{
	initial_poly problem;
	
	problem.particleDiameter = 2.0;		//set particle diameter
    problem.rhop = 2000.0;				//set particle density
    problem.en = 0.804;					//set restitution coefficient
    problem.K1 = 100000;				//set loading stiffness
    problem.K2 = 100000;				//set unloading stiffness
    problem.Kc = 0.0;					//set cohesive stiffness
    
    problem.mu_slid = 0.01;				//set sliding friction coefficient
    problem.mu_roll = 0.0;				//set rolling friction coefficient
    problem.mu_tor = 0.0;				//set torsional friction coefficient
    problem.Phic = 0.5;					// penetration DepthMax, the maximum depth of linear plastic-viscoelastic normal force
    problem.N = 4096;					// number of particles
    problem.volumeFraction = 0.68;		// initial volume fraction
    problem.poly = 3.0;					//polydispersity dmax/dmin
    // ----------------------------------------------------------------
    
    // assign problem configuration: volume fraction and name
    
	problem.setName("ini_nu_w1");
	problem.tmax = 2000;
    problem.setSaveCount(500);
    problem.eneFile.setSaveCount(500);
    
    problem.dataFile.setFileType(FileType::NO_FILE);
	//problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);

    
    problem.solve(argc,argv);
    
}
