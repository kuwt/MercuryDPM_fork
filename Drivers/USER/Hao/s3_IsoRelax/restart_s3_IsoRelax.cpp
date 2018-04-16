//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
//#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class relax: public Mercury3D{
public:
/* s3_IsoRelax will load the 1 or more configurations from s2_IsotropicCompression
 * and relax those samples in order to create the initial homogeneous jammed samples 
 * for the further shear process. Here we basiclly doing nothing but just let all the
 * particles cool down in the cubic box. Be aware that if the sample is compressed with
 * friciton on, the following process should be also with frition on, same for cohesion.
 * */    
    relax()
    {
        setName("mu0-phi-w1-relax");
        readRestartFile();
        setRestarted(true);
        particleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
        write(std::cout);
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
		//! particles properties and initial positions
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
        
        
		
        //! particleSpecies (d average = 2)
        particleSpecies->setDensity(rhop);
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(tc,en,mass);
        particleSpecies->setPlasticParameters(K1,K2,Kc,Phic);
		particleSpecies->setSlidingStiffness(2.0/10.0*particleSpecies->getLoadingStiffness());
		particleSpecies->setRollingStiffness(2.0/10.0*particleSpecies->getLoadingStiffness());
		particleSpecies->setTorsionStiffness(2.0/10.0*particleSpecies->getLoadingStiffness());
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
		
		
		
        for (int i=0; i < N; i++) {
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
		
		
		
		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;
		
        std::cout << "Lx = " << getXMax() << ", Ly = " << getYMax() << ", Lz = " << getZMax() << std::endl;
        std::cout << "nu = " << Vp/(getXMax()*getYMax()*getZMax()) << std::endl;
        std::cout << "output = " << getName() << std::endl;
        
        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		 
    }
    
   
		
	
    
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	relax problem;
	
	//!  --------------------------------------------------
   
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
    problem.poly = 1.0;					//set polydispersity
    problem.tmax = 10000;				//set simulation time
    //! ----------------------------------------------------------------

    problem.setName("mu0-phi-w1-relax");
    
    
    problem.setSaveCount(5000);
    problem.eneFile.setSaveCount(5000);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.eneFile.setFileType(FileType::ONE_FILE);
 
    problem.solve();
}
