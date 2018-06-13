
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/StressStrainControlBoundary.h>
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
        setRestarted(false);
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
	Vec3D lengthBox; //!Box length in x-y-z
	Vec3D centerBox; //!center position of the box
	Vec3D relativeToCenter; //!particle position relative to the center of box
	Matrix3D strainRate; //!strainrate tensor
    Matrix3D gainFactor; //!Stress control gain factors, which will be multiplied by timestep dt
    Matrix3D dstrainRate; //!Stress control strainrate tensor change per timestep
    Matrix3D stressGoal; //!Stress control constants
	Matrix3D stressTotal;  //!stress components calculation variables
	Matrix3D stressKinetic;  //!stress components calculation variables
	Matrix3D stressStatic;  //!stress components calculation variables

    
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
		lengthBox.X = getXMax()-getXMin();
     	lengthBox.Y = getYMax()-getYMin();
     	lengthBox.Z = getZMax()-getZMin();
		
		for (int i=1; i < (N-1); i++) {
			Rmin = std::min(Rmin,particleHandler.getObject(i)->getRadius());
			Rmax = std::max(Rmax,particleHandler.getObject(i)->getRadius());
		}
		
        for (int i=0; i < N; i++) {
			//vx = 2*(particleHandler.getObject(i)->getPosition().Y - getYMax()/2.0);
			//particleHandler.getObject(i)->setVelocity(Vec3D(strainRate.XY*getTimeStep()*vx, 0.0, 0.0));
			particleHandler.getObject(i)->setVelocity(Vec3D(0.0, 0.0, 0.0));
			particleHandler.getObject(i)->setSpecies(particleSpecies);
			Vp = Vp + particleHandler.getObject(i)->getVolume();
		}
		
        double particleDiameter = 2.0*Rmin;
        double mass = rhop*constants::pi*mathsFunc::cubic(particleDiameter)/6.0;
        double tc = std::sqrt( mass/2.0/K1*( mathsFunc::square(constants::pi) + mathsFunc::square(log(en) ) ));
        double velocity_xy = strainRate.XY * lengthBox.Y;
        setTimeStep(tc/50);
		setTimeMax(tmax);
        		
        
        
		
        //! particleSpecies    set the species parameters
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
		
		boundaryHandler.clear();								//! Delete all exist boundaries


		if ( stressGoal.XY!=0 && strainRate.XY !=0 || stressGoal.XX!=0 && strainRate.XX !=0 || stressGoal.YY!=0 && strainRate.YY !=0 || stressGoal.ZZ!=0 && strainRate.ZZ !=0)
		{
            logger(ERROR,"Wrong control combination, stressGoal and dot_strain can not be controlled at same time");
		}
		else if (stressGoal.XY==0 && strainRate.XY ==0)
		{
		//! Set up new box periodic boundaries, no simple shear activated
        PeriodicBoundary normWall;
        normWall.set(Vec3D(1.0, 0.0, 0.0), getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 1.0, 0.0), getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(normWall);
        normWall.set(Vec3D(0.0, 0.0, 1.0), getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(normWall);

            logger(INFO,"Shear rate is zero, set up normal walls for periodic box");
		}
		else if (stressGoal.XY!=0 || strainRate.XY !=0)
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

            logger(INFO,"Shear rate is none zero, set up Lees-Edwards periodic boundary in x-y directions");
		}else {
            logger(ERROR,"This should not happen; please check implementation of StressStrainControlBoundary::set");
        }

		std::cout << "N = " << N << std::endl;
		std::cout << "Rmin = " << Rmin << std::endl;
		std::cout << "Rmax = " << Rmax << std::endl;

        std::cout << "lengthBox.X = " << lengthBox.X << ", lengthBox.Y = " << lengthBox.Y << ", lengthBox.Z = " << lengthBox.Z << std::endl;
        std::cout << "shear velocity = " << velocity_xy << std::endl;
        std::cout << "nu = " << Vp/(getXMax()*getYMax()*getZMax()) << std::endl;
		std::cout << "gp = " << velocity_xy/getYMax() << std::endl;
        std::cout << "k1 = " << K1 << ", k2 = " << K2 << std::endl;
        std::cout << "output = " << getName() << std::endl;

        std::cout << "delta t = " << getTimeStep() << std::endl;
        std::cout << "saving time = " <<  dataFile.getSaveCount()*getTimeStep() << std::endl;

//        StressStrainControlBoundary stressStrainControlBoundary;
//        boundaryHandler.copyAndAddObject(stressStrainControlBoundary);
		 
    }
    
	void actionsAfterTimeStep() override {
		//start writing new version of control code

		lengthBox.X = (getXMax() - getXMin());
		lengthBox.Y = (getYMax() - getYMin());
		lengthBox.Z = (getZMax() - getZMin());
		centerBox.X = (getXMax() + getXMin()) / 2.0;
		centerBox.Y = (getYMax() + getYMin()) / 2.0;
		centerBox.Z = (getZMax() + getZMin()) / 2.0;//! Box length Lx, Lyand Lz, and center point C of the box
		static Mdouble integrated_velocity_XY = 0.0; //!Used for updating the shift of Lees-Edwards boundary

		if (stressGoal.XX == 0 && stressGoal.YY == 0 && stressGoal.ZZ == 0 &&
			stressGoal.XY == 0) //!this activate only the strainrate control in x-y-z
		{
			//!  Change the system size according to next time step
			setXMax(getXMax() + 0.5 * lengthBox.X * strainRate.XX * getTimeStep());
			setXMin(getXMin() - 0.5 * lengthBox.X * strainRate.XX * getTimeStep());
			setYMax(getYMax() + 0.5 * lengthBox.Y * strainRate.YY * getTimeStep());
			setYMin(getYMin() - 0.5 * lengthBox.Y * strainRate.YY * getTimeStep());
			setZMax(getZMax() + 0.5 * lengthBox.Z * strainRate.ZZ * getTimeStep());
			setZMin(getZMin() - 0.5 * lengthBox.Z * strainRate.ZZ * getTimeStep());

			//! Box length Lx, Lyand Lz for next time step
			lengthBox.X = (getXMax() - getXMin());
			lengthBox.Y = (getYMax() - getYMin());
			lengthBox.Z = (getZMax() - getZMin());

			//!  Move the boundaries to next time step
			if (strainRate.XY != 0) {
				//!  Move the Lees-Edwards boundary in z direction to next time step
				double velocity_xy = strainRate.XY * lengthBox.Y;
				integrated_velocity_XY += velocity_xy * getTimeStep();
				LeesEdwardsBoundary *leesEdwardsBoundary = dynamic_cast<LeesEdwardsBoundary *>(boundaryHandler.getObject(
						0));
				leesEdwardsBoundary->set(
						[velocity_xy](double time) { return integrated_velocity_XY; },
						[velocity_xy](double time UNUSED) { return 0; },
						getXMin(), getXMax(), getYMin(), getYMax());

				//!  Move the boundary in z direction to next time step
				PeriodicBoundary *normWall;

				normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(1));
				normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
			} else {
				//!  Move the boundary in x direction to next time step
				PeriodicBoundary *normWall;
				normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(0));
				normWall->set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());

				//!  Move the boundary in y direction to next time step
				normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(1));
				normWall->set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());

				//!  Move the boundary in z direction to next time step
				normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(2));
				normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
			}

			//!  Give the strain-rate for all particles and move them to next timestep before integration
			for (auto& p : particleHandler) {
				relativeToCenter.X = p->getPosition().X - centerBox.X;
				relativeToCenter.Y = p->getPosition().Y - centerBox.Y;
				relativeToCenter.Z = p->getPosition().Z - centerBox.Z;
				p->move(Vec3D(strainRate.XX * getTimeStep() * relativeToCenter.X + strainRate.XY * getTimeStep() * relativeToCenter.Y,
							  strainRate.YY * getTimeStep() * relativeToCenter.Y, strainRate.ZZ * getTimeStep() * relativeToCenter.Z));
			}
		} else {
			dstrainRate.setZero();
			stressTotal.setZero();
			stressStatic.setZero();
			stressKinetic.setZero();

			//!calculate stress for kinetic part
			stressKinetic = getKineticStress();

			//!calculate the static stress tensor
			stressStatic = getStaticStress();

			//! calculate the stress total and average over the volume
			stressTotal = getTotalStress();

			//! amount by which the strainrate tensor has to be changed
			if (stressGoal.XX != 0) {
				dstrainRate.XX = dstrainRate.XX + gainFactor.XX * getTimeStep() * (stressTotal.XX - stressGoal.XX);
				std::cout << "StressXX = " << stressTotal.XX << std::endl;
			}
			if (stressGoal.YY != 0) {
				dstrainRate.YY = dstrainRate.YY + gainFactor.YY * getTimeStep() * (stressTotal.YY - stressGoal.YY);
				std::cout << "StressYY = " << stressTotal.YY << std::endl;
			}
			if (stressGoal.ZZ != 0) {
				dstrainRate.ZZ = dstrainRate.ZZ + gainFactor.ZZ * getTimeStep() * (stressTotal.ZZ - stressGoal.ZZ);
				std::cout << "StressZZ = " << stressTotal.ZZ << std::endl;
			}
			if (stressGoal.XY != 0) {
				dstrainRate.XY = dstrainRate.XY + gainFactor.XY * getTimeStep() * (stressTotal.XY - stressGoal.XY);
                std::cout << "dstrainRate.XY = " << dstrainRate.XY << std::endl;
				std::cout << "StressXY = " << stressTotal.XY << std::endl;
				std::cout << "StressYX = " << stressTotal.YX << std::endl;
			}
		}

		//!  Update the strainrate tensor
		strainRate = strainRate + dstrainRate;

		//!  Change the system size according to next time step
		setXMax(getXMax() + 0.5 * lengthBox.X * strainRate.XX * getTimeStep());
		setXMin(getXMin() - 0.5 * lengthBox.X * strainRate.XX * getTimeStep());
		setYMax(getYMax() + 0.5 * lengthBox.Y * strainRate.YY * getTimeStep());
		setYMin(getYMin() - 0.5 * lengthBox.Y * strainRate.YY * getTimeStep());
		setZMax(getZMax() + 0.5 * lengthBox.Z * strainRate.ZZ * getTimeStep());
		setZMin(getZMin() - 0.5 * lengthBox.Z * strainRate.ZZ * getTimeStep());

		//! Box length Lx, Lyand Lz for next time step
		lengthBox.X = (getXMax() - getXMin());
		lengthBox.Y = (getYMax() - getYMin());
		lengthBox.Z = (getZMax() - getZMin());

		//!  Move the boundaries to next time step
		if (strainRate.XY != 0) {
			//!  Move the Lees-Edwards boundary in z direction to next time step
			double velocity_xy = strainRate.XY * lengthBox.Y;
			integrated_velocity_XY += velocity_xy * getTimeStep();
			LeesEdwardsBoundary *leesEdwardsBoundary = dynamic_cast<LeesEdwardsBoundary *>(boundaryHandler.getObject(
					0));
			leesEdwardsBoundary->set(
					[velocity_xy](double time) { return integrated_velocity_XY; },
					[velocity_xy](double time UNUSED) { return 0; },
					getXMin(), getXMax(), getYMin(), getYMax());

			//!  Move the boundary in z direction to next time step
			PeriodicBoundary *normWall;

			normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(1));
			normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());

		} else {
			//  Move the boundary in x direction to next time step
			PeriodicBoundary *normWall;
			normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(0));
			normWall->set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());

			//  Move the boundary in y direction to next time step
			normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(1));
			normWall->set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());

			//  Move the boundary in z direction to next time step
			normWall = dynamic_cast<PeriodicBoundary *>(boundaryHandler.getObject(2));
			normWall->set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
		}
		// Give the strain-rate for all particles and move them to next timestep before integration
		for (auto& p : particleHandler) {
			relativeToCenter.X = p->getPosition().X - centerBox.X;
			relativeToCenter.Y = p->getPosition().Y - centerBox.Y;
			relativeToCenter.Z = p->getPosition().Z - centerBox.Z;
			p->move(Vec3D(strainRate.XX * getTimeStep() * relativeToCenter.X + strainRate.XY * getTimeStep() * relativeToCenter.Y,
						  strainRate.YY * getTimeStep() * relativeToCenter.Y, strainRate.ZZ * getTimeStep() * relativeToCenter.Z));
			}
	}

	Mdouble getTotalVolume() const
	{
		/**\brief Get the total volume of the system.
 	*/

		return (DPMBase::getXMax() - DPMBase::getXMin()) * (DPMBase::getYMax() - DPMBase::getYMin()) *
				 (DPMBase::getZMax() - DPMBase::getZMin()); //!Get the total volume of the system
	}


 	Matrix3D getKineticStress()
	{
	/**\brief Calculate the kinetic stress tensor in the system
 	* \details The function calculate the kinetic stress tensor based on particle fluctuation velocity.
 	*/
		Matrix3D stressKinetic;
		Matrix3D Fij;
		Matrix3D Jij;
		stressKinetic.setZero();
		Jij.setZero();
		Fij.setZero();

		//!calculate stress for kinetic part

		for (auto& p : particleHandler) {
			//RHOP += rhop * constants::pi * mathsFunc::cubic(particleHandler.getObject(i)->getRadius() * 2) / 6.0;
			Jij.XX += p->getVelocity().X * p->getMass();
			Fij.XX += p->getVelocity().X * p->getVelocity().X * p->getMass();

			Jij.YY += p->getVelocity().Y * p->getMass();
			Fij.YY += p->getVelocity().Y * p->getVelocity().Y * p->getMass();

			Jij.ZZ += p->getVelocity().Z * p->getMass();
			Fij.ZZ += p->getVelocity().Z * p->getVelocity().Z * p->getMass();

			Fij.XY += p->getVelocity().X * p->getVelocity().Y * p->getMass();

			Fij.YX += p->getVelocity().Y * p->getVelocity().X * p->getMass();
		}

		stressKinetic.XX = Fij.XX - (Jij.XX * Jij.XX / DPMBase::getTotalMass());
		stressKinetic.YY = Fij.YY - (Jij.YY * Jij.YY / DPMBase::getTotalMass());
		stressKinetic.ZZ = Fij.ZZ - (Jij.ZZ * Jij.ZZ / DPMBase::getTotalMass());
		stressKinetic.XY = Fij.XY - (Jij.XX * Jij.YY / DPMBase::getTotalMass());
		stressKinetic.YX = Fij.YX - (Jij.XX * Jij.YY / DPMBase::getTotalMass());

		stressKinetic = stressKinetic/getTotalVolume();

		return stressKinetic;
	}
	Matrix3D getStaticStress()
	{
		/**\brief Calculate the static stress tensor in the system
         * \details The function calculate the static stress tensor based on particle contact force and
         * contact normal branch vector.
         */
		Matrix3D stressStatic;  //!stress components calculation variables
		stressStatic.setZero();

		for (auto i : interactionHandler) {
			stressStatic.XX += i->getForce().X * i->getNormal().X * i->getDistance();
			stressStatic.YY += i->getForce().Y * i->getNormal().Y * i->getDistance();
			stressStatic.ZZ += i->getForce().Z * i->getNormal().Z * i->getDistance();
			stressStatic.XY += i->getForce().X * i->getNormal().Y * i->getDistance();
			stressStatic.YX += i->getForce().Y * i->getNormal().X * i->getDistance();
		}

		stressStatic = stressStatic/getTotalVolume();

		return stressStatic;
	}

	Matrix3D getTotalStress()
	{
		/**\brief Calculate the Total stress tensor in the system
         * \details The function calculate the total stress tensor which is
         * the sum of kinetic and static stress tensors.
         */
		Matrix3D stressTotal;  //!stress components calculation variables
		stressTotal.setZero();


		stressTotal = getKineticStress() + getStaticStress();

		return stressTotal;
	}

};
        

int main(int argc UNUSED, char *argv[] UNUSED)
{
	std::string restartName ("StressStrainControlBoundarySelfTestInput");
	
	poly_simpleshear problem(restartName);
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
    problem.poly = 3;

    
    
    //! ----------------------------------------------------------------
    //!set up the stress tensor (activate stress control)
	problem.stressGoal.XX = 0.0;
	problem.stressGoal.YY = 0.0;
	problem.stressGoal.ZZ = 0.0;
	problem.stressGoal.XY = -2.0; //! the input for shear stress control should be always negative sign
	
	//!set up the strainrate tensor
	problem.strainRate.XX = 0.0;
    problem.strainRate.YY = 0.0;
    problem.strainRate.ZZ = 0.0;
    problem.strainRate.XY = 0.0;
    
    //!set up the gain factors for stress control
    problem.gainFactor.XX = 0.0001; //!these gain factors will be multiplied by dt, they are case dependent and set by user
    problem.gainFactor.YY = 0.0001;
    problem.gainFactor.ZZ = 0.0001;
    problem.gainFactor.XY = 0.0001;

    
    //!set up the simulation time and name of output
    problem.tmax = 0.1;
    problem.setName("s4_FreeShearStressStrainSelfTest");

    problem.setSaveCount(1000);
    problem.eneFile.setSaveCount(1000);
	problem.dataFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);

    problem.solve();
	//problem.writeRestartFile();
}
