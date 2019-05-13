
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
//#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Particles/LiquidFilmParticle.h"
class poly_simpleshear: public Mercury3D{
public:
/* 
 * In 's4-1_SimpleShearStrainRate', it will read the restart file for a relaxed
 * configuration from s3_IsoRelax and start doing the shear with stress control in yy direction,
 * shear will be controlled by strain-rate input which is du_x/dy = dot_strain_xy =gamma_dot
 * The shear is applied by strainrate tensor to all the individual particles (according to position) 
 * and the Lees Edwards boundary. the stress value will be set by stressYYGoal and kept constant after reaching
 * the steady state.
 * */    
    poly_simpleshear()
    {
        setName("Plane_shear-relax");
        readRestartFile();
        setRestarted(false);
        particleSpecies = dynamic_cast<LinearViscoelasticFrictionLiquidMigrationWilletSpecies*>(speciesHandler.getObject(0));
        double N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
			particleHandler.getObject(i)->setSpecies(particleSpecies);
		}
        //std::cout << "input = " << getName() << std::endl;
    }
    
    LinearViscoelasticFrictionLiquidMigrationWilletSpecies* particleSpecies;
    
    
    Mdouble volumeFraction, particleDiameter, rhop, K1, mu_slid, mu_roll,mu_tor, poly, slidingDissipation; //!material input parameters
    Mdouble tmax, dampingCoefficient = 0.0; //!material/simulation input parameters
    Mdouble dot_strain_xy, gx, gy, gz, Lx, Ly; //!shear strain rate variables
    Mdouble dV, dPressure, dVL, Volume, stressYY, stressYYGoal, slow_factor, velocity_yy, alpha, cap_factor; //!Stress control variables
    Mdouble stressYY_static, stressYY_kinetic, RHOP, Jy, Fy; //!stress components calculation variables
    
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
		
		
		//! particles properties and initial positions
		double N = particleHandler.getNumberOfObjects();
		Mdouble Vp = 0;
		
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
		
        
		
        double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = constants::pi/std::sqrt(K1/mass-mathsFunc::square(slidingDissipation/(2*mass)));
        
        setTimeStep(tc/50);
        //setTimeStep(1e-3);
		setTimeMax(tmax);
		setGravity(Vec3D(0.0,0.0,0.0));
     	double YMax =  getYMax();
     	Lx = getXMax()-getXMin();
     	Ly = getYMax()-getYMin();
        Mdouble velocity_xy = dot_strain_xy*Ly;

        for (int i=0; i < N; i++) {
			gx = (particleHandler.getObject(i)->getPosition().Y - getYMax()/2.0)/Ly;
			//particleHandler.getObject(i)->setVelocity(Vec3D(dot_strain_xy*getTimeStep()*gx, 0.0, 0.0));
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
        //! particleSpecies (d average = 2)
        particleSpecies->setDensity(rhop);
        particleSpecies->setSlidingFrictionCoefficient(mu_slid);
		particleSpecies->setSlidingFrictionCoefficientStatic(mu_slid);
		particleSpecies->setRollingFrictionCoefficient(mu_roll);
		particleSpecies->setRollingFrictionCoefficientStatic(mu_roll);
		particleSpecies->setTorsionFrictionCoefficient(mu_tor);
		particleSpecies->setTorsionFrictionCoefficientStatic(mu_tor);
		particleSpecies->setLiquidBridgeVolumeMax(90e-12);
		particleSpecies->setDistributionCoefficient(1.0);
		particleSpecies->setSurfaceTension(0.010);
		particleSpecies->setContactAngle(20.0*constants::pi/180.0);
		particleSpecies->setSlidingStiffness(120);
		
		boundaryHandler.clear();								//! Delete all exist boundaries
		
         //! Lees Edwards bc in y direction & periodic boundary in x direction
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
            [velocity_xy] (double time) { return time*velocity_xy; },
            [velocity_xy] (double time UNUSED) { return 0; },
            getXMin(),getXMax(),getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);
        
        //! periodic boundary in z direction
		PeriodicBoundary normWall;
        normWall.set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(normWall);		
		
		
		
		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;
		
        std::cout << "Lx = " << getXMax()-getXMin() << ", Ly = " << getYMax()-getYMin() << ", Lz = " << getZMax()-getZMin() << std::endl;
        std::cout << "Vwall = " << velocity_xy << std::endl;
        std::cout << "nu = " << Vp/((getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())) << std::endl;
        std::cout << "k1 = " << K1 << std::endl;
        std::cout << "output = " << getName() << std::endl;
        
        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		 
    }
    
	void actionsBeforeTimeLoop () override {
		// do a first force computation
		hGridActionsBeforeTimeLoop();
		checkAndDuplicatePeriodicParticles();
		hGridActionsBeforeTimeStep();
		computeAllForces();
		removeDuplicatePeriodicParticles();
		unsigned count = 0;
		for (BaseInteraction* j : interactionHandler) 
		{
			if (j->getOverlap()>=0) count++;
		}
		
		Mdouble value = 50e-12 * particleHandler.getNumberOfObjects() / count;
		for (BaseInteraction* j : interactionHandler)
		{
			LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(j);
			if (l==nullptr)
			{
				std::cerr << "Warning ParticleType needs to be LiquidFilmParticle" << std::endl;
				exit(-1);
			}
			if (l->getOverlap()>=0) l->setLiquidBridgeVolume(value);
			
		}
		
	}
    
    
    void actionsAfterTimeStep() override
    {
		
		//! pressure control in yy direction
		static Mdouble meanRadius =1.0;
		static Mdouble crossArea = constants::pi * mathsFunc::square(meanRadius);
		static Mdouble stiffness = particleSpecies->getSlidingStiffness();
		stressYY = 0.0;
		stressYY_static = 0.0;
		stressYY_kinetic = 0.0;
		RHOP = 0.0;
		Jy = 0.0;
		Fy = 0.0;
		Ly = getYMax()-getYMin();
		static Mdouble integrated_velocity_XY = 0.0;
		Mdouble velocity_xy = dot_strain_xy*Ly;
		integrated_velocity_XY +=velocity_xy*getTimeStep();
		Volume = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
			
			//!calculate stress_yy for kinetic part
			double N = particleHandler.getNumberOfObjects();
			for (int i=0; i < N; i++) {
			RHOP += rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0; 
			Jy += particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fy += particleHandler.getObject(i)->getVelocity().Y*particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			}
			stressYY_kinetic = Fy - (Jy*Jy/RHOP);
			
			//!calculate stress_yy for static part
			for (auto i : interactionHandler) {
				stressYY_static += i->getForce().Y * i->getNormal().Y * i->getDistance();		
			}
			
			//! calculate the stress_yy total and average over the volume
			stressYY = stressYY_kinetic + stressYY_static;
			stressYY /= Volume;
			
			//! amount by which the pressure has to be increased
			dPressure = stressYY - stressYYGoal;
			
			
			if (dPressure/stressYYGoal > cap_factor)
			{dPressure = stressYYGoal*cap_factor;}
			else if (dPressure/stressYYGoal < -cap_factor)
			{dPressure = -stressYYGoal*cap_factor;}
			else 
			{dPressure = dPressure;}
			
			if (getTime()> 50 && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 5e-1)
			{alpha = 0.00914, cap_factor = 0.02;} 
			else if(getTime()> 50 && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 1e-1 && dot_strain_xy < 5e-1)
			{alpha = 4.57, cap_factor = 0.02;}
			else if (getTime()> 100 && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy > 1e-2 && dot_strain_xy < 1e-1 )
			{alpha = 4.57, cap_factor = 0.02;}
			else if (getTime()> 200 && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 5e-3 && dot_strain_xy <= 1e-2 )
			{alpha = 4.57, cap_factor = 0.02;}
			else if (getTime()> 400 && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 1e-3 && dot_strain_xy < 5e-3 )
			{alpha = 0.457, cap_factor = 0.02;} 
			else if (getTime()> 800 && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy < 1e-3)
			{alpha = 0.0457, cap_factor = 0.02;} 
			else
			{alpha = alpha, cap_factor = cap_factor;}
			
			//! amount by which position should be changed to achieve the right pressure
			//dV = dPressure * crossArea / stiffness /getTimeStep();
			dV = dPressure * crossArea / stiffness;
			
			//! velocity in yy direction with pressure tune parameters
			//Mdouble velocity_yy = dV*fabs(dPressure/stressYYGoal);
			velocity_yy = dV/getTimeStep();
			//std::cout << "stressYY=" << stressYY << " stressYY_static=" << stressYY_static << " stressYY_kinetic=" << stressYY_kinetic<<  " dPressure=" << dPressure << " alpha=" << alpha <<" velocity_xy=" << velocity_xy <<  "Ly" << Ly <<std::endl;
			//std::cout << "integrated_velocity" << integrated_velocity_XY << std::endl;
			//std::cout << " velocity_xy=" << velocity_xy << std::endl;
			//std::cout << " move_y_old=" << slow_factor*velocity_yy*getTimeStep() << " move_y_new=" << dV*getTimeStep()/alpha<< std::endl;
		//! move all the particles in x and y direction   	
        for (int i=0; i < N; i++) {
			gx = (particleHandler.getObject(i)->getPosition().Y - (getYMax()+getYMin())/2.0)/Ly;
			gy = (particleHandler.getObject(i)->getPosition().Y - getYMin())/Ly;
			
			//!comment next line for volume const simple shear
			particleHandler.getObject(i)->move(Vec3D(velocity_xy*getTimeStep()*gx,dV*getTimeStep()*gy/alpha, 0.0));
			//!uncomment next line for volume const simple shear
			//particleHandler.getObject(i)->move(Vec3D(velocity_xy*getTimeStep()*gx,0.0, 0.0));
		}
		//setYMax(getYMax()+slow_factor*velocity_yy*getTimeStep());
		
		//!comment next line for volume const simple shear
			setYMax(getYMax()+dV*getTimeStep()/alpha);
			
			//! move the lees-edwards boundary in y direction
			LeesEdwardsBoundary* leesEdwardsBoundary = dynamic_cast<LeesEdwardsBoundary*>(boundaryHandler.getObject(0));
			leesEdwardsBoundary->set(
				[&] (double time) { return integrated_velocity_XY; },
				[] (double time UNUSED) { return 0; },
				getXMin(),getXMax(),getYMin(),getYMax());
        
           
	}

    
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	poly_simpleshear problem;
	
	
	//!  --------------------------------------------------
	
	problem.slidingDissipation = 0.0005; //set dissipation coefficient
	problem.particleDiameter = 0.0022;		//set particle diameter
    problem.rhop = 2000.0;				//set particle density
    problem.K1 = 120;				//set loading stiffness
    problem.mu_slid = 0.01;				//set sliding friction coefficient
    problem.mu_roll = 0.0;				//set rolling friction coefficient
    problem.mu_tor = 0.0;				//set torsional friction coefficient
    problem.poly = 2;		
    problem.stressYYGoal = 250;
    
    
    //! ----------------------------------------------------------------
    problem.alpha = 0.00914; //!choose between 5 to 0.01
	problem.cap_factor = 1.0;
    
    
    problem.setName("Plane_LiquidMigrations03_lbridge_new");
    problem.tmax = 50;
    problem.dot_strain_xy = 0.03;
    
    problem.setSaveCount(8000);
    problem.eneFile.setSaveCount(2000);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    problem.eneFile.setFileType(FileType::ONE_FILE);
 
    problem.solve();
}
