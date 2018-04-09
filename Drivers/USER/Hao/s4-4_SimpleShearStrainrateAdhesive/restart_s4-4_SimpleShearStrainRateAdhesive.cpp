
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"

class restart_simpleshear: public Mercury3D{
public:
    restart_simpleshear()
    {
        setName("p500");
        readRestartFile();
        setRestarted(false);
        speciesHandler.clear();
        particleSpecies = speciesHandler.copyAndAddObject(LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies());
        
        //particleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        //double N = particleHandler.getNumberOfObjects();
        //for (int i=0; i < N; i++) {
			//particleHandler.getObject(i)->setSpecies(particleSpecies);
		//}
		
		 for (BaseParticle* p : particleHandler)
            p->setSpecies(particleSpecies);
        //std::cout << "input = " << getName() << std::endl;
		
    }
    
   LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies* particleSpecies;
    
    
    Mdouble volumeFraction, particleDiameter, rhop, en, K1, K2, Kc, Phic, mu_slid, mu_roll,mu_tor, poly, K_adh, f_adh_max; //!material input parameters
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
        double tc = std::sqrt( mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));
        
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
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(tc,en,mass);
        particleSpecies->setPlasticParameters(K1,K2,Kc,Phic);
        particleSpecies->setAdhesionStiffness(K_adh);
        particleSpecies->setAdhesionForceMax(f_adh_max);
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
	}
	
	void actionsAfterTimeStep() override
    {
		
		//! pressure control in yy direction
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
	restart_simpleshear problem;
	
	problem.particleDiameter = 2.0;		//set particle diameter
    problem.rhop = 2000.0;				//set particle density
    problem.en = 0.804;					//set restitution coefficient
    problem.K1 = 100000;				//set loading stiffness
    problem.K2 = 100000;				//set unloading stiffness
    problem.Kc = 0.0;					//set cohesive stiffness
    problem.K_adh = 100000;
    problem.f_adh_max = 1;
    
    problem.mu_slid = 0.5;				//set sliding friction coefficient
    problem.mu_roll = 0.0;				//set rolling friction coefficient
    problem.mu_tor = 0.0;				//set torsional friction coefficient
    problem.Phic = 0.5;					// penetration DepthMax, the maximum depth of linear plastic-viscoelastic normal force
    problem.poly = 3;		
    problem.stressYYGoal = 500;
    
    
    //! ----------------------------------------------------------------
    problem.alpha = 0.00914; //!choose between 5 to 0.01
	problem.cap_factor = 1.0;
    
    
    problem.setName("p500-restart");
    problem.tmax = 10000;
    problem.dot_strain_xy = 5e-3;
    
    problem.setSaveCount(4000);
    problem.eneFile.setSaveCount(1000);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    //problem.dataFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    //problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);
 
    problem.solve();
}
