
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class poly_simpleshear: public Mercury3D{
public:
/* In 's4-3_SimpleShearMoveBoundary', it will read the restart file for a relaxed
 * configuration from s3_IsoRelax and start doing the shear with stress control in yy direction,
 * shear rate will be defined by strain-rate input which is du_x/dy = dot_strain_xy =gamma_dot
 * The shear is applied by pre-defined shear velocity_xy on the Lees Edwards boundary. 
 * the stress value will be set by stressYYGoal and kept constant after reaching the steady state.
 * */    
    poly_simpleshear(std::string restartName)
    {
        setName(restartName);
        readRestartFile();
        setRestarted(true);
        particleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
        
    }
    
    LinearPlasticViscoelasticFrictionSpecies* particleSpecies;
    
    Mdouble volumeFraction, particleDiameter, rhop, en, K1, K2, Kc, Phic, mu_slid, mu_roll,mu_tor, poly; //!material input parameters
    Mdouble tmax, dampingCoefficient = 0.0; //!material/simulation input parameters
    Mdouble dot_strain_xy, gx, gy, gz, Lx, Ly; //!shear strain rate variables
    Mdouble dV, dPressure, dVL, Volume, stressYY, stressYYGoal, slow_factor, velocity_yy, alpha, cap_factor, beta_stop; //!Stress control variables
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

    void setupInitialConditions() override
    {		
        double Rmin = particleHandler.getObject(0)->getRadius();
		double Rmax = particleHandler.getObject(0)->getRadius();
		//! particles properties and initial positions
		double N = particleHandler.getNumberOfObjects();
		Mdouble Vp = 0;
		Mdouble vx = 0;
		
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
		
        for (int i=0; i < N; i++) {
			vx = 2*(particleHandler.getObject(i)->getPosition().Y - getYMax()/2.0);
			particleHandler.getObject(i)->setVelocity(Vec3D(dot_strain_xy*getTimeStep()*vx, 0.0, 0.0));
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
        double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = std::sqrt( mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));
        
		double YMax = getYMax();
		Lx = getXMax()-getXMin();
     	Ly = getYMax()-getYMin();
		Mdouble velocity = dot_strain_xy*Ly;
        setTimeStep(tc/50);
		setTimeMax(tmax);
        		
        
        
		
        //! particleSpecies    set the species parameters
        particleSpecies->setDensity(rhop);
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(tc,en,mass);
        particleSpecies->setPlasticParameters(K1,K2,Kc,Phic);
		//particleSpecies->setDissipation(dissipation);
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
		setGravity(Vec3D(0.0,0.0,0.0));
		
		boundaryHandler.clear();								//! Delete all exist boundaries
		
         //! Lees Edwards bc in y direction & periodic boundary in x direction
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
            [velocity] (double time) { return time*velocity; },
            [velocity] (double time UNUSED) { return velocity; },
            getXMin(),getXMax(),getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);
        
        //! periodic boundary in z direction
		PeriodicBoundary normWall;
        normWall.set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(normWall);
			
		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;
		
        std::cout << "Lx = " << getXMax() << ", Ly = " << getYMax() << ", Lz = " << getZMax() << std::endl;
        std::cout << "Vwall = " << velocity << std::endl;
        std::cout << "nu = " << Vp/(getXMax()*getYMax()*getZMax()) << std::endl;
		std::cout << "gp = " << velocity/getYMax() << std::endl;
        std::cout << "k1 = " << K1 << "k2 = " << K2 << std::endl;
        std::cout << "output = " << getName() << std::endl;
        
        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		 
    }
     void actionsAfterTimeStep() override
    {
		
		
		static Mdouble meanRadius =1.0;
		static Mdouble crossArea = constants::pi * mathsFunc::square(meanRadius);
		static Mdouble stiffness = particleSpecies->getLoadingStiffness();
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
			
			if (getTime()> 0.005*getTimeMax() && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 5e-1)
			{alpha = 2*getTimeStep(), cap_factor = 0.02;} 
			else if(getTime()> 0.005*getTimeMax() && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 1e-1 && dot_strain_xy < 5e-1)
			{alpha = 1000*getTimeStep(), cap_factor = 0.02;}
			else if (getTime()> 0.01*getTimeMax() && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy > 1e-2 && dot_strain_xy < 1e-1 )
			{alpha = 1000*getTimeStep(), cap_factor = 0.02;}
			else if (getTime()> 0.02*getTimeMax() && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 5e-3 && dot_strain_xy <= 1e-2 )
			{alpha = 1000*getTimeStep(), cap_factor = 0.02;}
			else if (getTime()> 0.04*getTimeMax() && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy >= 1e-3 && dot_strain_xy < 5e-3 )
			{alpha = 100*getTimeStep(), cap_factor = 0.02;} 
			else if (getTime()> 0.08*getTimeMax() && dPressure > -0.001*stressYYGoal && dPressure < 0.001*stressYYGoal && dot_strain_xy < 1e-3)
			{alpha = 2*getTimeStep(), cap_factor = 0.02;} 
			else
			{alpha = alpha, cap_factor = cap_factor;}

			
			// amount by which position should be changed to achieve the right pressure
			dV = dPressure * crossArea / stiffness /getTimeStep();
			
			velocity_yy = dV/getTimeStep();
			//std::cout << "stressYY=" << stressYY << " stressYY_static=" << stressYY_static << " stressYY_kinetic=" << stressYY_kinetic<<  " dPressure=" << dPressure << " alpha=" << alpha <<std::endl;
			
			setYMax(getYMax()+dV*getTimeStep()/alpha);
			double YMax = getYMax();
			Ly = getYMax()-getYMin();

			Mdouble velocity = dot_strain_xy*Ly;
			LeesEdwardsBoundary* leesEdwardsBoundary = dynamic_cast<LeesEdwardsBoundary*>(boundaryHandler.getObject(0));
			leesEdwardsBoundary->set(
				[velocity] (double time) { return integrated_velocity_XY; },
				[velocity] (double time UNUSED) { return velocity; },
				getXMin(),getXMax(),getYMin(),getYMax());
        }			
};
        

int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::string restartName ("p250");
	
	poly_simpleshear problem(restartName);
	//!  --------------------------------------------------
    problem.particleDiameter = 2.0;		//set particle diameter
    problem.rhop = 2000.0;				//set particle density
    problem.en = 0.804;					//set restitution coefficient
    problem.K1 = 100000;				//set loading stiffness
    problem.K2 = 100000;				//set unloading stiffness
    problem.Kc = 0.0;					//set cohesive stiffness
    
    problem.mu_slid = 0.5;				//set sliding friction coefficient
    problem.mu_roll = 0.0;				//set rolling friction coefficient
    problem.mu_tor = 0.0;				//set torsional friction coefficient
    problem.Phic = 0.5;					// penetration DepthMax, the maximum depth of linear plastic-viscoelastic normal force
    problem.poly = 3;		

    
    
    //! ----------------------------------------------------------------
	problem.alpha = 0.00914; //!choose initially as 2*dt
	problem.cap_factor = 1.0;
    
	problem.stressYYGoal = 250;
    problem.dot_strain_xy = 5e-2;
    problem.tmax = 20000;
    problem.setName("p250");
    
    problem.setSaveCount(8000);
    problem.eneFile.setSaveCount(2000);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    //problem.dataFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    //problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);
    problem.solve();
}
