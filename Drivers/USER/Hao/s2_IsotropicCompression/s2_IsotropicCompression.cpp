//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
//#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class strain_iso: public Mercury3D{
public:
/*! s2_IsotropicCompression will load the configuration from s1_IniGeneration
 * The particles, with given volume fracion, are
 * placed in a cubic box of dimensions Xmax, Ymax, Zmax
 * and period boundaries are applied in all directions.
 * With a given strain-rate, all 6 period boundaries are moved inwards
 * isotropically to compress the particle gas to target volume fraction.
 * */    
    strain_iso(std::string restartName)
    {
        setName(restartName);
        readRestartFile();
        setRestarted(false);
        particleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
        std::cout << "input = " << getName() << std::endl;
    }
    
    LinearPlasticViscoelasticFrictionSpecies* particleSpecies;
    Mdouble tmax, N, Xp, Yp, Zp, Lx, Ly, Lz, Px, Py, Pz;
    Mdouble particleDiameter, rhop, en, K1, K2, Kc, Phic, mu_slid, mu_roll,mu_tor, poly;
    Mdouble dot_strain_xx, dot_strain_yy, dot_strain_zz;
    Mdouble nu_ini, nu_final, dot_strain_iso;
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
        double Rmin = particleHandler.getObject(0)->getRadius();
		double Rmax = particleHandler.getObject(0)->getRadius();
		double N = particleHandler.getNumberOfObjects();
		Mdouble Vp = 0;
		// particles properties and initial positions
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
        double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = std::sqrt( mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));
        

        setTimeStep(tc/50);
		setTimeMax(tmax);
        setGravity(Vec3D(0.0,0.0,0.0));		
        
        
		
        // particleSpecies (d average = 2)
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
		dampingCoefficient =  0.1*particleSpecies->getDissipation();
		
		//boundaryHandler.removeObject(0);
		
        // assign velocity to Lees Edward y-boundary
        //LeesEdwardsBoundary leesEdwardsBoundary;
        //leesEdwardsBoundary.set(
            //[velocity] (double time) { return time*velocity; },
            //[velocity] (double time UNUSED) { return velocity; },
            //getXMin(),getXMax(),getYMin(),getYMax());
        //boundaryHandler.copyAndAddObject(leesEdwardsBoundary);
        
		
		
        for (int i=0; i < N; i++) {
			//particleHandler.getObject(i)->setVelocity(Vec3D(0.0, 0.0, 0.0));
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
		dot_strain_iso = std::log(nu_ini/nu_final)/(3*tmax);   //calculate the isotropic strain rate
		dot_strain_xx = dot_strain_iso;
		dot_strain_yy = dot_strain_iso;
		dot_strain_zz = dot_strain_iso;
		std::cout << "dot_strain_xx " << dot_strain_xx << std::endl;
		
		
		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;
		
        std::cout << "Lx = " << getXMax()-getXMin() << ", Ly = " << getYMax()-getYMin() << ", Lz = " << getZMax()-getZMin() << std::endl;
        std::cout << "nu = " << Vp/((getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())) << std::endl;
        std::cout << "k1 = " << K1 << std::endl;
        std::cout << "output = " << getName() << std::endl;
        
        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		 
    }
    
     // Strain-rate application before integration
    void actionsBeforeTimeStep() override
    {
		Lx = (getXMax()-getXMin());					// Box length Lx, Lyand Lz in current time step
		Ly = (getYMax()-getYMin());
		Lz = (getZMax()-getZMin());
		Px = (getXMax()+getXMin())/2.0;
		Py = (getYMax()+getYMin())/2.0;
		Pz = (getZMax()+getZMin())/2.0;
		//std::cout << "Lx = " << Lx << std::endl;
		//std::cout << "Ly = " << Ly << std::endl;
		//std::cout << "Lz = " << Lz << std::endl;
		
        //  Change the system size according to next time step
		setXMax(getXMax()+0.5*Lx*dot_strain_xx*getTimeStep());
        setXMin(getXMin()-0.5*Lx*dot_strain_xx*getTimeStep());
        setYMax(getYMax()+0.5*Ly*dot_strain_yy*getTimeStep());
        setYMin(getYMin()-0.5*Ly*dot_strain_yy*getTimeStep());
        setZMax(getZMax()+0.5*Lz*dot_strain_zz*getTimeStep());
        setZMin(getZMin()-0.5*Lz*dot_strain_zz*getTimeStep());
        
        //  Give the strain-rate for all particles and move them to next timestep before integration
        N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
			Xp = particleHandler.getObject(i)-> getPosition().X - Px;
			Xp = particleHandler.getObject(i)-> getPosition().Y - Py;
			Zp = particleHandler.getObject(i)-> getPosition().Z - Pz;
			particleHandler.getObject(i)->move(Vec3D(dot_strain_xx*getTimeStep()*Xp,dot_strain_yy*getTimeStep()*Yp, dot_strain_zz*getTimeStep()*Zp));				
        
        
        //  Move the boundary in x direction to next time step
        PeriodicBoundary* normWall;
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(0));
        normWall->set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        
         //  Move the boundary in y direction to next time step
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(1));
        normWall->set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
  
        //  Move the boundary in z direction to next time step
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(2));
        normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
		
		
		}
		
	}
    
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::string restartName ("ini_nu_w1"); //the Prefix of your restart file from stage 1
	
	strain_iso problem(restartName);
	
	//  --------------------------------------------------
   
	problem.particleDiameter = 2.0;		//set particle diameter
    problem.rhop = 2000.0;				//set particle density
    problem.en = 0.804;					//set restitution coefficient
    problem.K1 = 100000;				//set loading stiffness
    problem.K2 = 100000;				//set unloading stiffness
    problem.Kc = 0.0;					//set cohesive stiffness
    
    problem.mu_slid = 0.0;				//set sliding friction coefficient
    problem.mu_roll = 0.0;				//set rolling friction coefficient
    problem.mu_tor = 0.0;				//set torsional friction coefficient
    problem.Phic = 0.5;					// penetration DepthMax, the maximum depth of linear plastic-viscoelastic normal force
    problem.poly = 1;					//set polydispersity
    problem.nu_ini = 0.40;				//set initial volume fraction
    problem.nu_final = 0.90;			//set final volume fraction
    problem.tmax = 3000;				//set simulation time
    //----------------------------------------------------------------

    problem.setName("mu0-0.4-0.9");
    
    
    problem.setSaveCount(500);
    problem.eneFile.setSaveCount(500);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.eneFile.setFileType(FileType::ONE_FILE);
 
    problem.solve();
}
