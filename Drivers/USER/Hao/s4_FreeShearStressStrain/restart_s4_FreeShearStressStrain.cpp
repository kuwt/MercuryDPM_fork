
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
        double N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
			particleHandler.getObject(i)->setSpecies(particleSpecies);
		}
        //std::cout << "input = " << getName() << std::endl;
    }
    
    LinearPlasticViscoelasticFrictionSpecies* particleSpecies;
    
    Mdouble volumeFraction, particleDiameter, rhop, en, K1, K2, Kc, Phic, mu_slid, mu_roll,mu_tor, poly; //!material input parameters
    Mdouble tmax, dampingCoefficient = 0.0; //!material/simulation input parameters
    Mdouble Lx, Ly, Lz, Cx, Cy, Cz, rpx, rpy, rpz; //!Box length in x-y-z, center position of the box and particle position relative to the center of box
    Mdouble dot_strain_xx, dot_strain_yy, dot_strain_zz, dot_strain_xy; //!strainrate tensor with four components
    Mdouble gx, gy, gz, velocity_xy, velocity_yy; //!shear strain rate variables
    Mdouble gain_xx, gain_yy, gain_zz, gain_xy, Volume; //!Stress control gain factors, which will be multiplied by time step dt
    Mdouble dStrain_xx, dStrain_yy, dStrain_zz, dStrain_xy; //!Stress control strainrate tensor change per time step
    Mdouble stressXXGoal,stressYYGoal,stressZZGoal,stressXYGoal;//!Stress control constants
    Mdouble stressXX,stressXX_static, stressXX_kinetic, stressYY,stressYY_static, stressYY_kinetic;
    Mdouble stressZZ,stressZZ_static, stressZZ_kinetic, stressXY,stressXY_static, stressXY_kinetic; //!stress components calculation variables
    Mdouble RHOP, Jx, Jy, Jz, Fx, Fy, Fz, Fxy; //!stress components calculation variables
    
    //! Add Background damping to the system
    void computeExternalForces (BaseParticle * CI) override	
	{
		DPMBase::computeExternalForces (CI);
		if (!CI->isFixed())
        {

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
		Lx = getXMax()-getXMin();
     	Ly = getYMax()-getYMin();
     	Lz = getZMax()-getZMin();
		
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
		
        for (int i=0; i < N; i++) {
			//vx = 2*(particleHandler.getObject(i)->getPosition().Y - getYMax()/2.0);
			//particleHandler.getObject(i)->setVelocity(Vec3D(dot_strain_xy*getTimeStep()*vx, 0.0, 0.0));
			particleHandler.getObject(i)->setVelocity(Vec3D(0.0, 0.0, 0.0));
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
        double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = std::sqrt( mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));
		double velocity_xy = dot_strain_xy*Ly;
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
		
		
		if ( stressXYGoal!=0 && dot_strain_xy !=0 || stressXXGoal!=0 && dot_strain_xx !=0 || stressYYGoal!=0 && dot_strain_yy !=0 || stressZZGoal!=0 && dot_strain_zz !=0)
		{
		std::cout << "wrong control combination, stressGoal and dot_strain can not be controlled at same time"  <<std::endl;
		return;
		}
		else if (stressXYGoal==0 && dot_strain_xy ==0)
		{
		//! Set up new box periodic boundaries, no simple shear activated
        PeriodicBoundary normWall;
        normWall.set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(normWall);
        
        std::cout << "Shear rate is zero, set up normal walls for periodic box" << std::endl;
        
		}
		else if (stressXYGoal!=0 || dot_strain_xy !=0)
		{
         //! Lees Edwards bc in y direction & periodic boundary in x direction, simple shear boundary activated
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
        
        std::cout << "Shear rate is none zero, set up Lees-Edwards periodic boundary in x-y directions" << std::endl;
		}
			
		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;
		
        std::cout << "Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << std::endl;
        std::cout << "shear velocity = " << velocity_xy << std::endl;
        std::cout << "nu = " << Vp/(getXMax()*getYMax()*getZMax()) << std::endl;
		std::cout << "gp = " << velocity_xy/getYMax() << std::endl;
        std::cout << "k1 = " << K1 << "k2 = " << K2 << std::endl;
        std::cout << "output = " << getName() << std::endl;
        
        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;
		 
    }
    
	void actionsAfterTimeStep() override
    {
		if (stressXYGoal==0 && dot_strain_xy ==0) //!this activate only the strainrate control in x-y-z
		{
		Lx = (getXMax()-getXMin()); Ly = (getYMax()-getYMin()); Lz = (getZMax()-getZMin());
		Cx = (getXMax()+getXMin())/2.0; Cy = (getYMax()+getYMin())/2.0; Cz = (getZMax()+getZMin())/2.0;//! Box length Lx, Lyand Lz, and center point C of the box
		Volume = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin()); //!Get the volume of the box for later averaging
		dStrain_xx = 0.0; dStrain_yy = 0.0; dStrain_zz = 0.0; dStrain_xy = 0.0;
		//std::cout << "Lx = " << Lx << std::endl;
		//std::cout << "Ly = " << Ly << std::endl;
		//std::cout << "Lz = " << Lz << std::endl;

       	if (stressXXGoal!=0 || stressYYGoal!=0 || stressZZGoal!=0)
       	{
		stressXX = 0.0; stressYY = 0.0; stressZZ = 0.0; stressXY = 0.0;
		stressXX_static = 0.0; stressYY_static = 0.0; stressZZ_static = 0.0; stressXY_static = 0.0;
		stressXX_kinetic = 0.0; stressYY_kinetic = 0.0; stressZZ_kinetic = 0.0; stressXY_kinetic = 0.0;
		RHOP = 0.0;
		Jx = 0.0; Jy = 0.0; Jz = 0.0;
		Fx = 0.0; Fy = 0.0; Fz = 0.0, Fxy = 0.0; 
		Lx = (getXMax()-getXMin()); Ly = (getYMax()-getYMin()); Lz = (getZMax()-getZMin());
		Cx = (getXMax()+getXMin())/2.0; Cy = (getYMax()+getYMin())/2.0; Cz = (getZMax()+getZMin())/2.0;//! Box length Lx, Lyand Lz, and center point C of the box
		
		//!calculate stress for kinetic part
			double N = particleHandler.getNumberOfObjects();
			for (int i=0; i < N; i++) {
			RHOP += rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0; 
			Jx += particleHandler.getObject(i)->getVelocity().X*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fx += particleHandler.getObject(i)->getVelocity().X*particleHandler.getObject(i)->getVelocity().X*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			Jy += particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fy += particleHandler.getObject(i)->getVelocity().Y*particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			Jz += particleHandler.getObject(i)->getVelocity().Z*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fz += particleHandler.getObject(i)->getVelocity().Z*particleHandler.getObject(i)->getVelocity().Z*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			Fxy += particleHandler.getObject(i)->getVelocity().X*particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			}
			stressXX_kinetic = Fx - (Jx*Jx/RHOP);
			stressYY_kinetic = Fy - (Jy*Jy/RHOP);
			stressZZ_kinetic = Fz - (Jz*Jz/RHOP);
			stressXY_kinetic = Fxy - (Jx*Jy/RHOP);
			
			//!calculate stress_yy for static part
			for (auto i : interactionHandler) {
				stressXX_static += i->getForce().X * i->getNormal().X * i->getDistance();
				stressYY_static += i->getForce().Y * i->getNormal().Y * i->getDistance();
				stressZZ_static += i->getForce().Z * i->getNormal().Z * i->getDistance();
				stressXY_static += i->getForce().X * i->getNormal().Y * i->getDistance();		
			}
			
			
			//! calculate the stress total and average over the volume
			stressXX = stressXX_kinetic + stressXX_static;
			stressXX /= Volume;
			stressYY = stressYY_kinetic + stressYY_static;
			stressYY /= Volume;
			stressZZ = stressZZ_kinetic + stressZZ_static;
			stressZZ /= Volume;
			stressXY = stressXY_kinetic + stressXY_static;
			stressXY /= Volume;

			
			//! amount by which the strainrate tensor has to be changed
			if (stressXXGoal!=0)
			{
			dStrain_xx = dStrain_xx + gain_xx*getTimeStep()*(stressXX - stressXXGoal);
			std::cout << "StressXX = " << stressXX << std::endl;
			}
			else if (stressYYGoal!=0)
			{
			dStrain_yy = dStrain_yy + gain_yy*getTimeStep()*(stressYY - stressYYGoal);
			std::cout << "StressYY = " << stressYY << std::endl;
			}
			else if (stressZZGoal!=0)
			{
			dStrain_zz = dStrain_zz + gain_zz*getTimeStep()*(stressZZ - stressZZGoal);
			std::cout << "StressZZ = " << stressZZ << std::endl;
			}
			else if (stressXYGoal!=0)
			{
			dStrain_xy = dStrain_xy + gain_xy*getTimeStep()*(stressXY - stressXYGoal);
			std::cout << "StressXY = " << stressXY << std::endl;
			}
			
		}
		//!  Update the strainrate tensor
		dot_strain_xx = dot_strain_xx + dStrain_xx;
		dot_strain_yy = dot_strain_yy + dStrain_yy;
		dot_strain_zz = dot_strain_zz + dStrain_zz;
		dot_strain_xy = dot_strain_xy + dStrain_xy;
		
		
        //!  Change the system size according to next time step
		setXMax(getXMax()+0.5*Lx*dot_strain_xx*getTimeStep());
        setXMin(getXMin()-0.5*Lx*dot_strain_xx*getTimeStep());
        setYMax(getYMax()+0.5*Ly*dot_strain_yy*getTimeStep());
        setYMin(getYMin()-0.5*Ly*dot_strain_yy*getTimeStep());
        setZMax(getZMax()+0.5*Lz*dot_strain_zz*getTimeStep());
        setZMin(getZMin()-0.5*Lz*dot_strain_zz*getTimeStep());
        
        //!  Move the boundary in x direction to next time step
        PeriodicBoundary* normWall;
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(0));
        normWall->set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        
         //!  Move the boundary in y direction to next time step
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(1));
        normWall->set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
  
        //!  Move the boundary in z direction to next time step
        normWall = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(2));
        normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
		
		//!  Give the strain-rate for all particles and move them to next time step before integration
        double N = particleHandler.getNumberOfObjects();
        for (int i=0; i < N; i++) {
			rpx = particleHandler.getObject(i)-> getPosition().X - Cx;
			rpy = particleHandler.getObject(i)-> getPosition().Y - Cy;
			rpz = particleHandler.getObject(i)-> getPosition().Z - Cz;
			particleHandler.getObject(i)->move(Vec3D(dot_strain_xx*getTimeStep()*rpx,dot_strain_yy*getTimeStep()*rpy, dot_strain_zz*getTimeStep()*rpz));				
        
		}	
	} 
	else
	{ 
		Lx = (getXMax()-getXMin()); Ly = (getYMax()-getYMin()); Lz = (getZMax()-getZMin());
		Cx = (getXMax()+getXMin())/2.0; Cy = (getYMax()+getYMin())/2.0; Cz = (getZMax()+getZMin())/2.0;//! Box length Lx, Lyand Lz, and center point C of the box
		
		dStrain_xx = 0.0; dStrain_yy = 0.0; dStrain_zz = 0.0; dStrain_xy = 0.0;
		static Mdouble integrated_velocity_XY = 0.0;
		Volume = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
		
		
			if (stressXXGoal!=0 || stressYYGoal!=0 || stressZZGoal!=0 || stressXYGoal!=0)
			{
			stressXX = 0.0; stressYY = 0.0; stressZZ = 0.0; stressXY = 0.0;
			stressXX_static = 0.0; stressYY_static = 0.0; stressZZ_static = 0.0; stressXY_static = 0.0;
			stressXX_kinetic = 0.0; stressYY_kinetic = 0.0; stressZZ_kinetic = 0.0; stressXY_kinetic = 0.0;
			RHOP = 0.0;
			Jx = 0.0; Jy = 0.0; Jz = 0.0;
			Fx = 0.0; Fy = 0.0; Fz = 0.0, Fxy = 0.0; 
		
			//!calculate stress for kinetic part
			double N = particleHandler.getNumberOfObjects();
			for (int i=0; i < N; i++) {
			RHOP += rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0; 
			Jx += particleHandler.getObject(i)->getVelocity().X*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fx += particleHandler.getObject(i)->getVelocity().X*particleHandler.getObject(i)->getVelocity().X*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			Jy += particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fy += particleHandler.getObject(i)->getVelocity().Y*particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			Jz += particleHandler.getObject(i)->getVelocity().Z*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			Fz += particleHandler.getObject(i)->getVelocity().Z*particleHandler.getObject(i)->getVelocity().Z*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			Fxy += particleHandler.getObject(i)->getVelocity().X*particleHandler.getObject(i)->getVelocity().Y*rhop*constants::pi*mathsFunc::cubic(particleHandler.getObject(i)->getRadius()*2)/6.0;
			
			}
			stressXX_kinetic = Fx - (Jx*Jx/RHOP);
			stressYY_kinetic = Fy - (Jy*Jy/RHOP);
			stressZZ_kinetic = Fz - (Jz*Jz/RHOP);
			stressXY_kinetic = Fxy - (Jx*Jy/RHOP);
			
			//!calculate stress_yy for static part
			for (auto i : interactionHandler) {
				stressXX_static += i->getForce().X * i->getNormal().X * i->getDistance();
				stressYY_static += i->getForce().Y * i->getNormal().Y * i->getDistance();
				stressZZ_static += i->getForce().Z * i->getNormal().Z * i->getDistance();
				stressXY_static += i->getForce().X * i->getNormal().Y * i->getDistance();		
			}
			
			
			//! calculate the stress total and average over the volume
			stressXX = stressXX_kinetic + stressXX_static;
			stressXX /= Volume;
			stressYY = stressYY_kinetic + stressYY_static;
			stressYY /= Volume;
			stressZZ = stressZZ_kinetic + stressZZ_static;
			stressZZ /= Volume;
			stressXY = stressXY_kinetic + stressXY_static;
			stressXY /= Volume;
			
			
			//! amount by which the strainrate tensor has to be changed
				if (stressXXGoal!=0)
				{
				dStrain_xx = dStrain_xx + gain_xx*getTimeStep()*(stressXX - stressXXGoal);
				std::cout << "StressXX = " << stressXX << std::endl;
				}
				else if (stressYYGoal!=0)
				{
				dStrain_yy = dStrain_yy + gain_yy*getTimeStep()*(stressYY - stressYYGoal);
				std::cout << "StressYY = " << stressYY << std::endl;
				}
				else if (stressZZGoal!=0)
				{
				dStrain_zz = dStrain_zz + gain_zz*getTimeStep()*(stressZZ - stressZZGoal);
				std::cout << "StressZZ = " << stressZZ << std::endl;
				}
				else if (stressXYGoal!=0)
				{
				dStrain_xy = dStrain_xy + gain_xy*getTimeStep()*(stressXY - stressXYGoal);
				std::cout << "StressXY = " << stressXY << std::endl;
				}
			}
			
			//!  Update the strainrate tensor
			dot_strain_xx = dot_strain_xx + dStrain_xx;
			dot_strain_yy = dot_strain_yy + dStrain_yy;
			dot_strain_zz = dot_strain_zz + dStrain_zz;
			dot_strain_xy = dot_strain_xy + dStrain_xy;
		
		
			//!  Change the system size according to next time step
			setXMax(getXMax()+0.5*Lx*dot_strain_xx*getTimeStep());
			setXMin(getXMin()-0.5*Lx*dot_strain_xx*getTimeStep());
			setYMax(getYMax()+0.5*Ly*dot_strain_yy*getTimeStep());
			setYMin(getYMin()-0.5*Ly*dot_strain_yy*getTimeStep());
			setZMax(getZMax()+0.5*Lz*dot_strain_zz*getTimeStep());
			setZMin(getZMin()-0.5*Lz*dot_strain_zz*getTimeStep());
			
			Lx = (getXMax()-getXMin());					//! Box length Lx, Lyand Lz for next time step
			Ly = (getYMax()-getYMin());
			Lz = (getZMax()-getZMin());

			double velocity_xy = dot_strain_xy*Ly;
			integrated_velocity_XY +=velocity_xy*getTimeStep();
			LeesEdwardsBoundary* leesEdwardsBoundary = dynamic_cast<LeesEdwardsBoundary*>(boundaryHandler.getObject(0));
			leesEdwardsBoundary->set(
				[velocity_xy] (double time) { return integrated_velocity_XY; },
				[velocity_xy] (double time UNUSED) { return 0; },
				getXMin(),getXMax(),getYMin(),getYMax());
				
			//!  Give the strain-rate for all particles and move them to next time step before integration
			double N = particleHandler.getNumberOfObjects();
			for (int i=0; i < N; i++) {
			rpx = particleHandler.getObject(i)-> getPosition().X - Cx;
			rpy = particleHandler.getObject(i)-> getPosition().Y - Cy;
			rpz = particleHandler.getObject(i)-> getPosition().Z - Cz;
			particleHandler.getObject(i)->move(Vec3D(dot_strain_xx*getTimeStep()*rpx+dot_strain_xy*getTimeStep()*rpy,dot_strain_yy*getTimeStep()*rpy, dot_strain_zz*getTimeStep()*rpz));				
        
		}	
        }
	
     
	}			
};
        

int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::string restartName ("p500");
	
	poly_simpleshear problem(restartName);
	//!  --------------------------------------------------
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
    problem.poly = 3;		

    
    
    //! ----------------------------------------------------------------
    //!set up the stress tensor (activate stress control)
	problem.stressXXGoal = 500;
	problem.stressYYGoal = 0.0;
	problem.stressZZGoal = 0.0;
	problem.stressXYGoal = 0.0; //! the input for shear stress control should be always negative sign
	
	//!set up the strainrate tensor
	problem.dot_strain_xx = 0.0005;
    problem.dot_strain_yy = 0.0;
    problem.dot_strain_zz = 0.0;
    problem.dot_strain_xy = 0.0;
    
    //!set up the gain factors for stress control
    problem.gain_xx = 0.0001; //!these gain factors will be multiplied by dt, they are case dependent and set by user
    problem.gain_yy = 0.0001;
    problem.gain_zz = 0.0001;
    problem.gain_xy = 0.0001;
    
    //!set up the simulation time and name of output
    problem.tmax = 100;
    problem.setName("p500");
    
    problem.setSaveCount(1000);
    problem.eneFile.setSaveCount(1000);
    problem.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    //problem.dataFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    //problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);
    problem.solve();
}
