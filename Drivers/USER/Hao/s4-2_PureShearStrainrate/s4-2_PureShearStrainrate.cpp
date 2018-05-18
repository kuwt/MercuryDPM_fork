
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
//#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class strain_pureshear: public Mercury3D{
public:
/* 
 * In 's4-2_PureShearStrainRate', it will read the restart file for a relaxed
 * configuration from s3_IsoRelax and start doing the pure shear with volume of box = const,
 * shear will be controlled by strain-rate tensor with dot_epsilon_xx = - dot_epsilon_yy, dot_epsilon_zz = 0.
 * Then shear is applied to all the individual particles (according to position) 
 * and the 4 periodic boundaries: x-y-directions. 
 * Be aware the box can not be sheared forever as in simple shear mode,
 * because with large shear strain, the whole box will be streched into a
 * rectangle and thus less particles  will consist in x-direction (compression)
 * and finally collapse when particles are having self contacts.
 * */    
    strain_pureshear(std::string restartName)
    {
        setName(restartName);								//  read last configuration from restart file
        readRestartFile();
        setRestarted(false);									//  set particle species
        particleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
        std::cout << "input = " << getName() << std::endl;
    }
    
     LinearPlasticViscoelasticFrictionSpecies* particleSpecies;
     
    Mdouble tmax, N, Xp, Yp, Zp, Lx, Ly, Lz, Px, Py, Pz;				//  Initialize all relaive parameters
    Mdouble particleDiameter, rhop, en, K1, K2, Kc, Phic, mu_slid, mu_roll,mu_tor, poly, dissipation, volumeFraction;
    Mdouble dot_strain_xx, dot_strain_yy, dot_strain_zz;
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

    void setupInitialConditions()								//  set up initial conditions
    {		
		
		//  Initialize particle related parameters : Radius: Rmin Rmax, Diameter, mass, collision time tc
		double Rmin = particleHandler.getObject(0)->getRadius();	
		double Rmax = particleHandler.getObject(0)->getRadius();
		double N = particleHandler.getNumberOfObjects();	//  number of particles
		Mdouble Vp = 0;										//  volume of particles
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());			//  get the minimum and maximum particle size in the system
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
		
        double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = std::sqrt( mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));
       
       //  set the time step tc for integration and tmax for total simulation time
        setTimeStep(tc/50);
        //setTimeStep(0.004);
		setTimeMax(tmax);
		
       
       
		
        //  particleSpecies    set the species parameters
        particleSpecies->setDensity(rhop);
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(tc,en,mass);
        particleSpecies->setPlasticParameters(K1,K2,Kc,Phic);
		//particleSpecies->setDissipation(dissipation);
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
		setGravity(Vec3D(0.0,0.0,0.0));
		
		
		//  set particles properties and initial positions
	
		
        for (int i=0; i < N; i++) {
			//Xp = particleHandler.getObject(i)-> getPosition().X - (getXMax()+getXMin())/2.0;    //  relative particle position in x direction
			//Zp = particleHandler.getObject(i)-> getPosition().Z - (getZMax()+getZMin())/2.0;	//  relative particle position in z direction
			//particleHandler.getObject(i)->setVelocity(Vec3D(dot_strain*Xp, 0.0, dot_strain*Zp));  //  set new initial velocity from relative position
			//particleHandler.getObject(i)->setVelocity(Vec3D(0.0, 0.0, 0.0));
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();						//  calculate the total particle volume
		}
		
		
		
		//double L = std::pow(Vp/volumeFraction,1.0/3.0);			//  Get the box size from volume fraction						
		//setXMax(L);												//  Set the initial system size
        //setYMax(L);
        //setZMax(L);
        //setXMin(0.0);
        //setYMin(0.0);
        //setZMin(0.0);
        
        //interactionHandler.clear();
        boundaryHandler.clear();								//  Delete all exist boundaries
        
        
        //  Set up new box periodic boundaries
        PeriodicBoundary normWall;
        normWall.set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(normWall);
		
			
		//  Print the important infos for the further check
		
		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;	
        std::cout << "Lx = " << getXMax()-getXMin() << ", Ly = " << getYMax()-getYMin() << ", Lz = " << getZMax()-getZMin() << std::endl;
        std::cout << "nu = " << Vp/((getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())) << std::endl;
        std::cout << "k = " << K1 << std::endl;
        std::cout << "output = " << getName() << std::endl;
        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		
    }
    
    
    //  Strain-rate application before integration
    void actionsAfterTimeStep() override
    {
		Lx = (getXMax()-getXMin());					//  Box length Lx, Lyand Lz in current time step
		Ly = (getYMax()-getYMin());
		Lz = (getZMax()-getZMin());
		Px = (getXMax()+getXMin())/2.0;
		Py = (getYMax()+getYMin())/2.0;
		Pz = (getZMax()+getZMin())/2.0;
		//std::cout << "Lx = " << Lx << std::endl;
		//std::cout << "Ly = " << Ly << std::endl;
		//std::cout << "Lz = " << Lz << std::endl;

       	
		
        //   Change the system size according to next time step
		setXMax(getXMax()+0.5*Lx*dot_strain_xx*getTimeStep());
        setXMin(getXMin()-0.5*Lx*dot_strain_xx*getTimeStep());
        setYMax(getYMax()+0.5*Ly*dot_strain_yy*getTimeStep());
        setYMin(getYMin()-0.5*Ly*dot_strain_yy*getTimeStep());
        setZMax(getZMax()+0.5*Lz*dot_strain_zz*getTimeStep());
        setZMin(getZMin()-0.5*Lz*dot_strain_zz*getTimeStep());
        
        //   Move the boundary in x direction to next time step
        PeriodicBoundary* normWall;
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(0));
        normWall->set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        
         //   Move the boundary in y direction to next time step
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(1));
        normWall->set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
  
        //   Move the boundary in z direction to next time step
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(2));
        normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
		
		//   Give the strain-rate for all particles and move them to next time step before integration
        N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
			Xp = particleHandler.getObject(i)-> getPosition().X - Px;
			Xp = particleHandler.getObject(i)-> getPosition().Y - Py;
			Zp = particleHandler.getObject(i)-> getPosition().Z - Pz;
			//particleHandler.getObject(i)->setVelocity(Vec3D(dot_strain*Xp+velocity_x, velocity_y, dot_strain*Zp+velocity_z));
			particleHandler.getObject(i)->move(Vec3D(dot_strain_xx*getTimeStep()*Xp,dot_strain_yy*getTimeStep()*Yp, dot_strain_zz*getTimeStep()*Zp));				
        
		}
		
	}
	
   
    
   
};

int main(int argc UNUSED, char *argv[] UNUSED)
{	
	std::string restartName ("mu1-250-relax");
	strain_pureshear problem(restartName);
	
	//   --------------------------------------------------
	// 	particle properties
    problem.rhop = 2000.0;
    problem.en = 0.804;					//set restitution coefficient
    problem.K1 = 100000;				//set loading stiffness
    problem.K2 = 100000;				//set unloading stiffness
    problem.Kc = 0.0;					//set cohesive stiffness
    
    problem.mu_slid = 0.5;				//set sliding friction coefficient
    problem.mu_roll = 0.0;				//set rolling friction coefficient
    problem.mu_tor = 0.0;				//set torsional friction coefficient
    problem.Phic = 0.5;					// penetration DepthMax, the maximum depth of linear plastic-viscoelastic normal force
  
    
    //  ----------------------------------------------------------------
    
    //  Simulation related settings
    problem.setName("mu1-nu0.6047");
    problem.tmax = 10000;
    problem.dot_strain_xx = 5e-4; //strain rate of pure shear
    problem.dot_strain_yy = -5e-4;  //strain rate of pure shear
    problem.dot_strain_zz = 0.0;
    
    problem.setSaveCount(4000);
    problem.eneFile.setSaveCount(2000);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.eneFile.setFileType(FileType::ONE_FILE);
 
    problem.solve();
}
