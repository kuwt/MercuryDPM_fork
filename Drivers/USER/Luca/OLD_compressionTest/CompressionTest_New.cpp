#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

/*
 * ToDo:
 * - Add comments throughout the code
 * - Separate into .h, .cc and .cpp
 * - Check the following since it doesn't make sense on execution:
 * 		tC_BB/dt: 51.44
 *		tC_SS/dt: 53.79
 *		tC_BW/dt: 27.02
 *		tC_SW/dt: 27.02
 *		tC_BS/dt: 53.79
 * - Add a verbose version (not urgent)
 * - implement and check the torque test
 * - add evaluation of phi_min and add a check for phi > phi_min
 * - introduce the inertial number and re-evaluate the piston velocity from there
 * - get rid of k2maxcalibratedarray and every reference to it
 * - put an if/else in the main switching petween k1 calibration (compression) and kC (torque)
 * - add a version of maxPressureRelDeviation for the torque
 * - PUT THE TORQUE CALIBRATION IN ANOTHER DRIVER! DOESN-T BELONG HERE!
 * - CLEAR EVERY REFERENCE TO TORQUE, ROTATION AND SO ON
 * - - unify settling and non-settling problem (maybe generate a driver on its own)
 */
 
 // STAGE 10 IS THE "DEAD" STAGE FOR QUITTING AND REBOOTING
 
 // 01.24 - Substituted kNew = kMean to kNew = (kFinal - kInitial)/2 (and included the computation of kMean)
 /* 
  * 01.26 - Changed teh piston velocity and introduced the input parameter setPistonVelocityScalingFactor
  * 01.26 - Now the system waits also after heightTargetMet, not only after pressureTargetMet
  * 01.26 - changed teh whole tuning routine and implemented the new kC calibration 
  *  
  * 02.20 - debugged displacement and deviation checks
  * 02.20 - introduced .cdat output and changed kdata output
  * 02.20 - introduced phi auto-tuning
  * 02.20 - added computeVariations() and updateVariationArrays() functions
  */
 

class CompressionTest_calibrationRoutine : public Mercury3D
{
	
	private:
	
	void setupInitialConditions() override
    {
		stage = 1;
		nBigInserted = 0;
		nSmallInserted = 0;
		calibrationStep = 0;
		t0 = getTime();
		
		waitForParticleSettling = false;
		isHeightTargetMet = false;
		isPressureTargetMet = false;
		
		std::cout << "SETTING SPECIES..." << std::endl;
		setSpecies();
		std::cout << "MAKING CASING..." << std::endl;
		makeCasing();
		std::cout << "COMPUTING NUMBER OF PARTICLES..." << std::endl;
		computeNumberOfParticles();
		std::cout << "N. OF BIG PARTICLES NEEDED: " << nBig << std::endl;
		std::cout << "N. OF SMALL PARTICLES NEEDED: " << nSmall << std::endl;
		std::cout << "TOTAL PARTICLES NEEDED: " << nSmall + nBig << std::endl;
		
		std::cout << std::endl << "PARTICLE INSERTION" << std::endl << std::endl;
		
		stage++;
    }
    
    void actionsOnRestart() override
    {
		stage = 3;
		calibrationStep = 0;
		t0 = getTime();
		setTimeMax(10.0);
		
		//restartSpecies();
		setSpecies();
		
		isHeightTargetMet = false;
		isPressureTargetMet = false;
	}
    
    void actionsAfterTimeStep() override
    {	
		if (stage == 2) insertParticles();
		if (stage == 3) 
		{
			if (isOnlySettlingFlag)
			{
				removeTopParticleLayer();
				
				std::cout << "PARTICLE SETTLING TERMINATED" << std::endl;
				std::cout << "QUITTING..." << std::endl;
				setTimeMax(getTime() + 10.0*getTimeStep());
				stage = 5;	
			}
			else
			{
				makePiston();
				
				meanK1Big = 0.0;
				setTime(0.0);
				
				// TEST
				for (int i=0; i < nCalibrationDataPoints; i++) std::cout << heightSteps[i] << "   " << pressureSteps[i] << std::endl;
				
				updateTargetDataPoint();
				resetPistonVelocity();
				makeCdatFile();
				
				t0 = getTime();
				stage++;
			}
		}
		
		if (stage >= 4) 
		{
			if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) writeToCdatFile();
		}
		
		if (stage == 4)
		{
			movePiston();
			
			makeDataAnalysis();
			runDataPointsCalibration();

			// TEST DYNAMIC OUTPUT
			if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) 
			{
				std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << 
			", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight << 
			", v = " << pistonVelocity << std::endl <<
			"P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << " , P_target = " << targetPressure << 
			", P/P_target = " << (pistonPressure + 0.01)/(targetPressure + 0.01) << ", |(P - P_target)/P_target| = " << fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << std::endl <<
			"cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap << 
			", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap << 
			", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
			"k1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() << 
			", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kC = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl;
			
			if (isPhiDynamicallyTunedFlag) std::cout << "phi_min = " << phiMin << std::endl;
			
			std::cout << std::endl;
			}
		}
	}
    
    void actionsAfterSolve() override
    {
		if (!isOnlySettlingFlag)
		{			
			delete [] k1CalibratedArray;
			delete [] k2MaxCalibratedArray;
			delete [] heightSteps;
			delete [] pressureSteps;
						
			cdatFile.close();
		}
	}
    
    public:
    
    
    // FUNCTIONS CALLED IN MAIN   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
    
    void setCdatOutputTimeInterval(double dt)
    {
        cdatOutputTimeInterval = dt;
    }
    
    void setEnergyRatioTolerance(double eRatio)
    {
        energyRatioTolerance = eRatio;
    }
    
    void setCasingProperties(double radius, double height, double density)
    {
		casingRadius = radius;
		casingHeight = height;
		densityWall = density;
		
		setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        
        setXMax(casingRadius);
        setYMax(casingRadius);
        setZMax(casingHeight);		
	}
	
	void setPowderBedProperties(double height, double pf)
    {
		particleBedHeight = height;
		particleBedPackingFraction = pf;
	}
    
    void setParticleProperties(double rB, double rS, double dB, double dS, double rhoB, double rhoS, double massRatio)
    {
		radiusBig = rB;
		radiusSmall = rS;
		sizeDispersityBig = dB;
		sizeDispersitySmall = dS;
		densityBig = rhoB;
		densitySmall = rhoS;
		smallToBigMassRatio = massRatio;
		
		volumeBig = 4.0*constants::pi*pow(radiusBig,3.0)/3.0;
		volumeSmall = 4.0*constants::pi*pow(radiusSmall,3.0)/3.0;
		massBig = volumeBig*densityBig;
		massSmall = volumeSmall*densitySmall;
	}
	
    void setParticleWallSlidingFrictionCoefficients(double bigWallMu, double smallWallMu)
    {
        bigWallSlidingFrictionCoeff = bigWallMu;
        smallWallSlidingFrictionCoeff = smallWallMu;
    }
    
    void setParticleWallRollingFrictionCoefficients(double bigWallMu, double smallWallMu)
    {
        bigWallRollingFrictionCoeff = bigWallMu;
        smallWallRollingFrictionCoeff = smallWallMu;
    }
    
    void setParticleWallTorsionFrictionCoefficients(double bigWallMu, double smallWallMu)
    {
        bigWallTorsionFrictionCoeff = bigWallMu;
        smallWallTorsionFrictionCoeff = smallWallMu;
    }
    
    void setParticleParticleSlidingFrictionCoefficients(double bigBigMu, double smallSmallMu, double bigSmallMu)
    {
        bigBigSlidingFrictionCoeff = bigBigMu;
        smallSmallSlidingFrictionCoeff = smallSmallMu;
        bigSmallSlidingFrictionCoeff = bigSmallMu;
    }

    void setParticleParticleRollingFrictionCoefficients(double bigBigMu, double smallSmallMu, double bigSmallMu)
    {
        bigBigRollingFrictionCoeff = bigBigMu;
        smallSmallRollingFrictionCoeff = smallSmallMu;
        bigSmallRollingFrictionCoeff = bigSmallMu;
    }
  
    void setParticleParticleTorsionFrictionCoefficients(double bigBigMu, double smallSmallMu, double bigSmallMu)
    {
        bigBigTorsionFrictionCoeff = bigBigMu;
        smallSmallTorsionFrictionCoeff = smallSmallMu;
        bigSmallTorsionFrictionCoeff = bigSmallMu;
    }
    
    void setWallStiffnessAndRestitutionCoefficients(double kW, double eB, double eS)
    {
        k1Wall = kW;
        bigWallRestitutionCoeff = eB;
        smallWallRestitutionCoeff = eS;
    }
    
    void setParticleParticleRestitutionCoefficients(double bigBigE, double smallSmallE, double bigSmallE)
    {
        bigBigRestitutionCoeff = bigBigE;
        smallSmallRestitutionCoeff = smallSmallE;
        bigSmallRestitutionCoeff = bigSmallE;
    }
    
    void setBigParticlePlasticProperties(double k1, double k2max, double kC, double phi)
    {
        k1Big = k1;
        k2MaxBig = k2max;
        kCBig = kC;
        phiBig = phi;
    }
    
    void setSmallParticlePlasticProperties(double k1, double k2max, double kC, double phi)
    {
        k1Small = k1;
        k2MaxSmall = k2max;
        kCSmall = kC;
        phiSmall = phi;
    }    

	void setCalibrationDataPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double pauseDuration)
    {
        nCalibrationDataPoints = n;
        pauseAfterDataPointMet = pauseDuration;
        k1CalibratedArray = new double[n];
        k2MaxCalibratedArray = new double[n];
        heightSteps = new double[n];
        pressureSteps = new double[n];
        
        for (int i = 0; i < n; i++)
        {
            heightSteps[i] = heightArrayPointer[i];
            pressureSteps[i] = pressureArrayPointer[i];
        }
    }

	void setCalibrationMaximumRelativeDeviations(double maxH, double maxP)
	{
		maxHeightRelDeviation = maxH;
		maxPressureRelDeviation = maxP;
	}

	void setCalibrationRoutineParameters(int nRun, int nMax, double *k1Array, double *kCArray, double *phiArray, double *devArray, double *dispArray)
	{
		nMaximumRoutineRuns = nMax;
		nCurrentRoutineRun = nRun;
		
		pointerToInitialK1Array = k1Array;
		pointerToInitialKCArray = kCArray;
		pointerToInitialPhiArray = phiArray;
		
		pointerToDeviationArray = devArray;
		pointerToDisplacementArray = dispArray;
	}
	
	void setNTimeStepsBetweenStiffnesAdjustments(double dt)
	{
		nTimeStepsBetweenStiffnessAdjustments = dt;
	}
	
	void setGlobalFlags(bool settlingFlag, bool phiFlag)
	{
		isOnlySettlingFlag = settlingFlag;
		isPhiDynamicallyTunedFlag = phiFlag;
	}
	
	void setPistonVelocityScalingFactor(double factor)
	{
		pistonVelocityScalingFactor = factor;
	}
	
	void setPistonRotation(double angularVelocity)
	{
		pistonAngularVelocity = angularVelocity;
		isPistonRotatingFlag = true;
	}
	
	void setCConstantAndDeltaHatMax(double C, double Dmax)
	{
		cConstant = C;
		deltaHatMax = Dmax;
	}
	
	
	// FUNCTIONS CALLED IN THE CLASS   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
	
	void setSpecies()
    {
        speciesHandler.clear();

        // BIG-BIG
        speciesBig = new LinearPlasticViscoelasticFrictionSpecies;
        speciesBig -> setDensity(densityBig);
        speciesBig -> setStiffnessAndRestitutionCoefficient(k1Big, bigBigRestitutionCoeff, massBig);
        speciesBig -> setUnloadingStiffnessMax(k2MaxBig);
        speciesBig -> setCohesionStiffness(kCBig);
        speciesBig -> setPenetrationDepthMax(phiBig);

        speciesBig -> setSlidingFrictionCoefficient(bigBigSlidingFrictionCoeff);
        speciesBig -> setSlidingStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
        speciesBig -> setSlidingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setRollingFrictionCoefficient(bigBigRollingFrictionCoeff);
        speciesBig -> setRollingStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
        speciesBig -> setRollingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setTorsionFrictionCoefficient(bigBigTorsionFrictionCoeff);
        speciesBig -> setTorsionStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
        speciesBig -> setTorsionDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesBig);
        
        // SMALL-SMALL
        speciesSmall = new LinearPlasticViscoelasticFrictionSpecies;
        speciesSmall -> setDensity(densitySmall);
        speciesSmall -> setStiffnessAndRestitutionCoefficient(k1Small, smallSmallRestitutionCoeff, massSmall);
        speciesSmall -> setUnloadingStiffnessMax(k2MaxSmall);
        speciesSmall -> setCohesionStiffness(kCSmall);
        speciesSmall -> setPenetrationDepthMax(phiSmall);

        speciesSmall -> setSlidingFrictionCoefficient(smallSmallSlidingFrictionCoeff);
        speciesSmall -> setSlidingStiffness(speciesSmall -> getLoadingStiffness()*2.0/7.0);
        speciesSmall -> setSlidingDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesSmall -> setRollingFrictionCoefficient(smallSmallRollingFrictionCoeff);
        speciesSmall -> setRollingStiffness(speciesSmall -> getLoadingStiffness()*2.0/7.0);
        speciesSmall -> setRollingDissipation(speciesSmall -> getDissipation()*2.0/7.0);        
        speciesSmall -> setTorsionFrictionCoefficient(smallSmallTorsionFrictionCoeff);
        speciesSmall -> setTorsionStiffness(speciesSmall -> getLoadingStiffness()*2.0/7.0);
        speciesSmall -> setTorsionDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesSmall);
        
        // WALL-WALL
        speciesWall = new LinearPlasticViscoelasticFrictionSpecies;
        speciesWall -> setDensity(densityWall);
        speciesWall -> setStiffnessAndRestitutionCoefficient(k1Wall, 1.0, massSmall);
        speciesWall -> setUnloadingStiffnessMax(k1Wall);
        speciesWall -> setCohesionStiffness(0.0);
        speciesWall -> setPenetrationDepthMax(0.0);
        
        speciesWall -> setSlidingFrictionCoefficient(0.0);
        speciesWall -> setSlidingStiffness(speciesWall -> getLoadingStiffness()*2.0/7.0);
        speciesWall -> setSlidingDissipation(speciesWall -> getDissipation()*2.0/7.0);
        speciesWall -> setRollingFrictionCoefficient(0.0);
        speciesWall -> setRollingStiffness(speciesWall -> getLoadingStiffness()*2.0/7.0);
        speciesWall -> setRollingDissipation(speciesWall -> getDissipation()*2.0/7.0);
        speciesWall -> setTorsionFrictionCoefficient(0.0);
        speciesWall -> setTorsionStiffness(speciesWall -> getLoadingStiffness()*2.0/7.0);
        speciesWall -> setTorsionDissipation(speciesWall -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesWall);
        
        // BIG-WALL
        //speciesMixedBigWall = speciesHandler.getMixedObject(speciesBig, speciesWall);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Wall), bigWallRestitutionCoeff, massBig);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setUnloadingStiffnessMax(k2MaxBig);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kCBig);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setPenetrationDepthMax(phiBig);

        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingFrictionCoefficient(bigWallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingFrictionCoefficient(bigWallRollingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionFrictionCoefficient(bigWallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        
        // SMALL-WALL
        //speciesMixedSmallWall = speciesHandler.getMixedObject(speciesSmall, speciesWall);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Small + k1Wall), smallWallRestitutionCoeff, massSmall);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setUnloadingStiffnessMax(k2MaxSmall);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setCohesionStiffness(kCSmall);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setPenetrationDepthMax(phiSmall);

        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingFrictionCoefficient(smallWallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingFrictionCoefficient(smallWallRollingFrictionCoeff);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionFrictionCoefficient(smallWallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);
        
        // BIG-SMALL
        //speciesMixedBigSmall = speciesHandler.getMixedObject(speciesBig, speciesSmall);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Small), bigSmallRestitutionCoeff, 0.5*(massBig + massSmall));
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setUnloadingStiffnessMax(0.5*(k2MaxBig + k2MaxSmall));
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kCBig + kCSmall));
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setPenetrationDepthMax(0.5*(phiBig + phiSmall));
        
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingFrictionCoefficient(bigSmallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingFrictionCoefficient(bigSmallRollingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionFrictionCoefficient(bigSmallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
		
		std::cout << "BIG-BIG stiffness and dissipation coefficients: " << speciesBig -> getLoadingStiffness() << " " << speciesBig -> getDissipation() << "\n";
		std::cout << "BIG-BIG friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
		std::cout << "BIG-BIG tangential stiffnesses: " << speciesBig -> getSlidingStiffness() << " " << speciesBig -> getRollingStiffness() << " " << speciesBig -> getTorsionStiffness() << "\n";
		std::cout << "BIG-BIG tangential dissipation coefficients: " << speciesBig -> getSlidingDissipation() << " " << speciesBig -> getRollingDissipation() << " " << speciesBig -> getTorsionDissipation() << "\n";
		std::cout << "BIG-BIG collision time: " << std::setprecision(4) << speciesBig -> getCollisionTime(massBig) << "\n\n";
		
		std::cout << "SMALL-SMALL stiffness and dissipation coefficients: " << speciesSmall -> getLoadingStiffness() << " " << speciesSmall -> getDissipation() << "\n";
		std::cout << "SMALL-SMALL friction coefficients: " << smallSmallSlidingFrictionCoeff << " " << smallSmallRollingFrictionCoeff << " " << smallSmallTorsionFrictionCoeff << "\n";
		std::cout << "SMALL-SMALL tangential stiffnesses: " << speciesSmall -> getSlidingStiffness() << " " << speciesSmall -> getRollingStiffness() << " " << speciesSmall -> getTorsionStiffness() << "\n";
		std::cout << "SMALL-SMALL tangential dissipation coefficients: " << speciesSmall -> getSlidingDissipation() << " " << speciesSmall -> getRollingDissipation() << " " << speciesSmall -> getTorsionDissipation() << "\n";
		std::cout << "SMALL-SMALL collision time: " << std::setprecision(4) << speciesSmall -> getCollisionTime(massSmall) << "\n\n";
		
		std::cout << "BIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation() << "\n";
		std::cout << "BIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
		std::cout << "BIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionStiffness() << "\n";
		std::cout << "BIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionDissipation() << "\n";
		std::cout << "BIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(massBig) << "\n\n";
		
		std::cout << "SMALL-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation() << "\n";
		std::cout << "SMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
		std::cout << "SMALL-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getTorsionStiffness() << "\n";
		std::cout << "SMALL-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getTorsionDissipation() << "\n";
		std::cout << "SMALL-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getCollisionTime(massSmall) << "\n\n";
		
		std::cout << "BIG-SMALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation() << "\n";
		std::cout << "BIG-SMALL friction coefficients: " << bigSmallSlidingFrictionCoeff << " " << bigSmallRollingFrictionCoeff << " " << bigSmallTorsionFrictionCoeff << "\n";
		std::cout << "BIG-SMALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getTorsionStiffness() << "\n";
		std::cout << "BIG-SMALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getTorsionDissipation() << "\n";
		std::cout << "BIG-SMALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getCollisionTime(0.5*(massBig + massSmall)) << "\n\n";
		
		std::cout << "tC_BB/dt: " << std::setprecision(4) << speciesBig -> getCollisionTime(massBig)/getTimeStep() << "\n";
		std::cout << "tC_SS/dt: " << std::setprecision(4) << speciesSmall -> getCollisionTime(massSmall)/getTimeStep() << "\n";
		std::cout << "tC_BW/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(massBig)/getTimeStep() << "\n";
		std::cout << "tC_SW/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getCollisionTime(massSmall)/getTimeStep() << "\n";
		std::cout << "tC_BS/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getCollisionTime(0.5*(massBig + massSmall))/getTimeStep() << "\n\n";
    }
    
    /*void restartSpecies()
    {
		speciesBig = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
		speciesSmall = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1));
		speciesWall = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(2));
		speciesMixedBigWall = speciesHandler.getMixedObject(speciesBig, speciesWall);
		speciesMixedSmallWall = speciesHandler.getMixedObject(speciesSmall, speciesWall);
		speciesMixedBigSmall = speciesHandler.getMixedObject(speciesBig, speciesSmall);
		
		std::cout << "BIG-BIG stiffness and dissipation coefficients: " << speciesBig -> getLoadingStiffness() << " " << speciesBig -> getDissipation() << "\n";
		std::cout << "BIG-BIG friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
		std::cout << "BIG-BIG tangential stiffnesses: " << speciesBig -> getSlidingStiffness() << " " << speciesBig -> getRollingStiffness() << " " << speciesBig -> getTorsionStiffness() << "\n";
		std::cout << "BIG-BIG tangential dissipation coefficients: " << speciesBig -> getSlidingDissipation() << " " << speciesBig -> getRollingDissipation() << " " << speciesBig -> getTorsionDissipation() << "\n";
		std::cout << "BIG-BIG collision time: " << std::setprecision(4) << speciesBig -> getCollisionTime(massBig) << "\n\n";
		
		std::cout << "SMALL-SMALL stiffness and dissipation coefficients: " << speciesSmall -> getLoadingStiffness() << " " << speciesSmall -> getDissipation() << "\n";
		std::cout << "SMALL-SMALL friction coefficients: " << smallSmallSlidingFrictionCoeff << " " << smallSmallRollingFrictionCoeff << " " << smallSmallTorsionFrictionCoeff << "\n";
		std::cout << "SMALL-SMALL tangential stiffnesses: " << speciesSmall -> getSlidingStiffness() << " " << speciesSmall -> getRollingStiffness() << " " << speciesSmall -> getTorsionStiffness() << "\n";
		std::cout << "SMALL-SMALL tangential dissipation coefficients: " << speciesSmall -> getSlidingDissipation() << " " << speciesSmall -> getRollingDissipation() << " " << speciesSmall -> getTorsionDissipation() << "\n";
		std::cout << "SMALL-SMALL collision time: " << std::setprecision(4) << speciesSmall -> getCollisionTime(massSmall) << "\n\n";
		
		std::cout << "BIG-WALL stiffness and dissipation coefficients: " << speciesMixedBigWall -> getLoadingStiffness() << " " << speciesMixedBigWall -> getDissipation() << "\n";
		std::cout << "BIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
		std::cout << "BIG-WALL tangential stiffnesses: " << speciesMixedBigWall -> getSlidingStiffness() << " " << speciesMixedBigWall -> getRollingStiffness() << " " << speciesMixedBigWall -> getTorsionStiffness() << "\n";
		std::cout << "BIG-WALL tangential dissipation coefficients: " << speciesMixedBigWall -> getSlidingDissipation() << " " << speciesMixedBigWall -> getRollingDissipation() << " " << speciesMixedBigWall -> getTorsionDissipation() << "\n";
		std::cout << "BIG-WALL collision time: " << std::setprecision(4) << speciesMixedBigWall -> getCollisionTime(massBig) << "\n\n";
		
		std::cout << "SMALL-WALL stiffness and dissipation coefficients: " << speciesMixedSmallWall -> getLoadingStiffness() << " " << speciesMixedSmallWall -> getDissipation() << "\n";
		std::cout << "SMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
		std::cout << "SMALL-WALL tangential stiffnesses: " << speciesMixedSmallWall -> getSlidingStiffness() << " " << speciesMixedSmallWall -> getRollingStiffness() << " " << speciesMixedSmallWall -> getTorsionStiffness() << "\n";
		std::cout << "SMALL-WALL tangential dissipation coefficients: " << speciesMixedSmallWall -> getSlidingDissipation() << " " << speciesMixedSmallWall -> getRollingDissipation() << " " << speciesMixedSmallWall -> getTorsionDissipation() << "\n";
		std::cout << "SMALL-WALL collision time: " << std::setprecision(4) << speciesMixedSmallWall -> getCollisionTime(massSmall) << "\n\n";
		
		std::cout << "BIG-SMALL stiffness and dissipation coefficients: " << speciesMixedBigSmall -> getLoadingStiffness() << " " << speciesMixedBigSmall -> getDissipation() << "\n";
		std::cout << "BIG-SMALL friction coefficients: " << bigSmallSlidingFrictionCoeff << " " << bigSmallRollingFrictionCoeff << " " << bigSmallTorsionFrictionCoeff << "\n";
		std::cout << "BIG-SMALL tangential stiffnesses: " << speciesMixedBigSmall -> getSlidingStiffness() << " " << speciesMixedBigSmall -> getRollingStiffness() << " " << speciesMixedBigSmall -> getTorsionStiffness() << "\n";
		std::cout << "BIG-SMALL tangential dissipation coefficients: " << speciesMixedBigSmall -> getSlidingDissipation() << " " << speciesMixedBigSmall -> getRollingDissipation() << " " << speciesMixedBigSmall -> getTorsionDissipation() << "\n";
		std::cout << "BIG-SMALL collision time: " << std::setprecision(4) << speciesMixedBigSmall -> getCollisionTime(0.5*(massBig + massSmall)) << "\n\n";
		
		std::cout << "tC_BB/dt: " << std::setprecision(4) << speciesBig -> getCollisionTime(massBig)/getTimeStep() << "\n";
		std::cout << "tC_SS/dt: " << std::setprecision(4) << speciesSmall -> getCollisionTime(massSmall)/getTimeStep() << "\n";
		std::cout << "tC_BW/dt: " << std::setprecision(4) << speciesMixedBigWall -> getCollisionTime(massBig)/getTimeStep() << "\n";
		std::cout << "tC_SW/dt: " << std::setprecision(4) << speciesMixedSmallWall -> getCollisionTime(massSmall)/getTimeStep() << "\n";
		std::cout << "tC_BS/dt: " << std::setprecision(4) << speciesMixedBigSmall -> getCollisionTime(0.5*(massBig + massSmall))/getTimeStep() << "\n\n";
	}*/
    
	void computeNumberOfParticles()
    {
		double particleBedVolume = 1.1*constants::pi*pow(casingRadius,2.0)*particleBedHeight;
        nSmall = (particleBedPackingFraction*particleBedVolume/volumeSmall)*(smallToBigMassRatio*densityBig/(densitySmall + smallToBigMassRatio*densityBig));
        nBig = (particleBedPackingFraction*particleBedVolume - nSmall*volumeSmall)/volumeBig;
	}
	
	void makeCasing()
    {
        wallHandler.clear();
        
        base.setSpecies(speciesWall);
        base.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
        wallHandler.copyAndAddObject(base);
        
        top.setSpecies(speciesWall);
        top.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
        wallHandler.copyAndAddObject(top);
        
        AxisymmetricIntersectionOfWalls(casing);
        casing.setSpecies(speciesWall);
        casing.setPosition(Vec3D(0.0,0.0,0.0));
        casing.setOrientation(Vec3D(0.0,0.0,1.0));
        casing.addObject(Vec3D(1.0,0.0,0.0),Vec3D(casingRadius,0.0,0.0));
        wallHandler.copyAndAddObject(casing);
    }
	
	bool particleInsertionSuccessful(bool isBig)
	{
		int insertionFailCounter = 0;
		Vec3D particlePosition;
		SphericalParticle p0;
		
		p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
		
		if (isBig) 
		{
			p0.setRadius(radiusBig*(1.0 + sizeDispersityBig*random.getRandomNumber(-1.0,1.0)));
			p0.setSpecies(speciesBig);
		}
		else
		{
			p0.setRadius(radiusSmall*(1.0 + sizeDispersitySmall*random.getRandomNumber(-1.0,1.0)));
			p0.setSpecies(speciesSmall);
		}
		
		do
		{
			particlePosition.X = (casingRadius - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
			particlePosition.Y = sqrt(pow(casingRadius - 1.01*(p0.getRadius()), 2.0) - pow(particlePosition.X, 2.0))*random.getRandomNumber(-1.0,1.0);
			particlePosition.Z = random.getRandomNumber(particleBedHeight + 1.01*(p0.getRadius()),casingHeight - 1.01*(p0.getRadius()));
			
			p0.setPosition(particlePosition);
			
			if (checkParticleForInteraction(p0)) 
			{
				particleHandler.copyAndAddObject(p0);
				return true;
			}
			
			insertionFailCounter++;
		}
		while (insertionFailCounter < 1000);
		
		return false;
	}
	
	void insertParticles()
	{		
		while(!waitForParticleSettling && (nBigInserted < nBig || nSmallInserted < nSmall))
		{
			if (nBigInserted < nBig && nSmallInserted < nSmall)
			{
				if (random.getRandomNumber(0.0,1.0) < nBig/(nBig + nSmall))
				{
					if (particleInsertionSuccessful(true)) 
					{
						nBigInserted++;
					}
					else 
					{
						t0 = getTime();
						waitForParticleSettling = true;
					}
				}
				else
				{
					if (particleInsertionSuccessful(false)) 
					{
						nSmallInserted++;
					}
					else 
					{
						t0 = getTime();
						waitForParticleSettling = true;
					}
				}
			}
			else if (nBigInserted < nBig && nSmallInserted >= nSmall)
			{
				if (particleInsertionSuccessful(true)) 
				{
					nBigInserted++;
				}
				else 
				{
					t0 = getTime();
					waitForParticleSettling = true;
				}
			}
			else if (nBigInserted >= nBig && nSmallInserted < nSmall)
			{
				if (particleInsertionSuccessful(false)) 
				{
					nSmallInserted++;
				}
				else 
				{
					t0 = getTime();
					waitForParticleSettling = true;
				}
			}
			else 
			{
				t0 = getTime();
				waitForParticleSettling = true;
			}
		}
		
		if (waitForParticleSettling && getTime() > t0 + sqrt(2.0*(casingHeight - particleBedHeight)/9.81)) waitForParticleSettling = false;
		if (!waitForParticleSettling && nBigInserted >= nBig && nSmallInserted >= nSmall && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance) 
		{
			std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
			stage++;
		}
	}

	void removeTopParticleLayer()
	{
		std::cout << std::endl << "REMOVING TOP PARTICLE LAYER" << std::endl;
		std::cout << "INITIAL NUMBER OF PARTICLES: " << particleHandler.getNumberOfObjects() << std::endl;
		
		for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            if (particleHandler.getObject(i) -> getPosition().Z + particleHandler.getObject(i) -> getRadius() > particleBedHeight) particleHandler.removeObject(i);
        }
        
        std::cout << "NEW NUMBER OF PARTICLES: " << particleHandler.getNumberOfObjects() << std::endl << std::endl;
	}
	
	void makePiston()
    {
		std::cout << std::endl << "CREATING PISTON" << std::endl << std::endl;
		
        pistonHeight = particleBedHeight + radiusSmall;
        pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

        compressionPiston.setSpecies(speciesWall);
        compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
        pistonPointer = wallHandler.copyAndAddObject(compressionPiston);
    }
	
	void makeDataAnalysis()
    {
        meanCoordinationNumber = 0.0;
        
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();
        
        pistonForce = 0.0;
        pistonPressure = 0.0;
        pistonTorque = 0.0;
        basePressure = 0.0;
        cylinderPressure = 0.0;
        meanTotalRelativeOverlap = 0.0;
        maxTotalRelativeOverlap = 0.0;
        meanPistonRelativeOverlap = 0.0;
        maxPistonRelativeOverlap = 0.0;
        meanBaseRelativeOverlap = 0.0;
        maxBaseRelativeOverlap = 0.0;
        meanCylinderRelativeOverlap = 0.0;
        maxCylinderRelativeOverlap = 0.0;
        
        int totalInteractionCounter = 0;
        int pistonInteractionCounter = 0;
        int baseInteractionCounter = 0;
        int cylinderInteractionCounter = 0;
        
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            
            // piston interactions
            if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
            {
                pistonPressure -= ((*i) -> getForce()).Z;
                
                meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                if (isPistonRotatingFlag) pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;
					
                pistonInteractionCounter++;
            }
            
            // base interactions
            if ((*i) -> getI() -> getIndex() == base.getIndex())
            {
                basePressure += ((*i) -> getForce()).Z;
                
                meanBaseRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxBaseRelativeOverlap) maxBaseRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                baseInteractionCounter++;
            }
            
            // cylinder interactions
            if ((*i) -> getI() -> getIndex() == casing.getIndex())
            {
                cylinderPressure -= (((*i) -> getContactPoint()).X * ((*i) -> getForce()).X + ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).Y)/casingRadius;
                // sqrt(pow(((*i) -> getContactPoint()).X,2.0) + pow(((*i) -> getContactPoint()).Y,2.0))
                
                meanCylinderRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxCylinderRelativeOverlap) maxCylinderRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                cylinderInteractionCounter++;
            }
            
            totalInteractionCounter++;
        }
        
        pistonForce = pistonPressure;
        pistonPressure /= constants::pi*pow(casingRadius,2.0);
        basePressure /= constants::pi*pow(casingRadius,2.0);
        if (pistonHeight > particleBedHeight) cylinderPressure /= 2.0*constants::pi*casingRadius*particleBedHeight;
        else cylinderPressure /= 2.0*constants::pi*casingRadius*pistonHeight;
        meanTotalRelativeOverlap /= totalInteractionCounter;
        meanPistonRelativeOverlap /= pistonInteractionCounter;
        meanBaseRelativeOverlap /= baseInteractionCounter;
        meanCylinderRelativeOverlap /= cylinderInteractionCounter;
    }
	
    void updateTargetDataPoint()
    {
		targetHeight = heightSteps[calibrationStep];
		targetPressure = pressureSteps[calibrationStep];
			
		if (calibrationStep > 0) previousTargetHeight = heightSteps[calibrationStep - 1];
		else previousTargetHeight = pistonHeight;
		
		std::cout << "PREVIOUS PISTON HEIGHT: " << previousTargetHeight << std::endl;
		std::cout << "NEXT PISTON TARGET HEIGHT: " << targetHeight << std::endl;
		std::cout << "NEXT PISTON TARGET PRESSURE: " << targetPressure << std::endl << std::endl;
		
		isHeightTargetMet = false;
        isPressureTargetMet = false;
	}
	
    void resetPistonVelocity()
	{
		if (previousTargetHeight < targetHeight) pistonVelocity = pistonMaximumVelocity;
		else pistonVelocity = -pistonMaximumVelocity;
		
		std::cout << "RESETTING PISTON VELOCITY. NEW VELOCITY: " << pistonVelocity << std::endl << std::endl;
	}
	
	void movePiston()
    {
		pistonHeight += pistonVelocity*getTimeStep();
		pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
		pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));     
    }
    
    void rotatePiston()
    {
		pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,pistonAngularVelocity));
		pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));        
    }
    
    void refreshSpecies()
    {
		// BIG-BIG
        speciesBig -> setStiffnessAndRestitutionCoefficient(k1Big, bigBigRestitutionCoeff, massBig);
        speciesBig -> setUnloadingStiffnessMax(k2MaxBig);
        speciesBig -> setCohesionStiffness(kCBig);
        speciesBig -> setPenetrationDepthMax(phiBig);

        speciesBig -> setSlidingStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
        speciesBig -> setSlidingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setRollingStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
        speciesBig -> setRollingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setTorsionStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
        speciesBig -> setTorsionDissipation(speciesBig -> getDissipation()*2.0/7.0);       
        
        // BIG-WALL
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Wall), bigWallRestitutionCoeff, massBig);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setUnloadingStiffnessMax(k2MaxBig);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kCBig);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setPenetrationDepthMax(phiBig);

        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        
        // BIG-SMALL
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Small), bigSmallRestitutionCoeff, 0.5*(massBig + massSmall));
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setUnloadingStiffnessMax(0.5*(k2MaxBig + k2MaxSmall));
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kCBig + kCSmall));
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setPenetrationDepthMax(0.5*(phiBig + phiSmall));

        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
	}
	
	void tuneK1()
	{
		if (pistonPressure < targetPressure) k1Big *= 1.001;
		else k1Big *= 0.999;
		
		if (isPhiDynamicallyTunedFlag) tunePhi();
		
		refreshSpecies();

		if (k1Big > k2MaxBig) 
		{
			std::cout << "\n\nFATAL ERROR: K1_BIG > K2MAX_BIG. QUITTING.";
			setTimeMax(getTime() + 10.0*getTimeStep());
			stage = 10;
		}
		
		if (k1Big < 0.0) 
		{
			std::cout << "\n\nFATAL ERROR: K1_BIG < 0.0. QUITTING.";
			setTimeMax(getTime() + 10.0*getTimeStep());
			stage = 10;
		}
	}
	
	void tunePhi()
	{
		phiBig = pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax/(cConstant - 1.0);
		phiMin = (1.0 - k1Big/k2MaxBig)*deltaHatMax;
	}
	
   void runDataPointsCalibration()
   {
		// reverse and halve the velocity if the target height was passed
        if (!isHeightTargetMet && ((pistonVelocity < 0.0 && previousTargetHeight < targetHeight) || (pistonVelocity > 0.0 && previousTargetHeight > targetHeight))) pistonVelocity *= -0.5;

		// upon reaching the target height update the bool and sets the velocity to 0
        if (!isHeightTargetMet && !isPressureTargetMet)
        {
			if (fabs((pistonHeight - targetHeight)/targetHeight) > maxHeightRelDeviation) movePiston();
			else
			{
				std::cout << "Height target reached.\n";
				std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
                t0 = getTime();
                    
				pistonVelocity = 0.0;
				isHeightTargetMet = true;
			}
        }

		// when the height is fine deal with the pressure
        if (isHeightTargetMet && !isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
        {
			if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
			{
				// rotates teh piston in torque mode
				if (isPistonRotatingFlag) rotatePiston();
				
				// the stiffness is changed every n time steps
				if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneK1();
			}
			else
			{
				std::cout << "Pressure target reached, now chillin' for " << pauseAfterDataPointMet << " seconds\n";
				t0 = getTime();
				isPressureTargetMet = true;
			}
		}
  
        // switch to next calibration point if the test of the previous is passed
        if (isHeightTargetMet && isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
        {
			// updates the mean value of the stiffnesses at the calibrated data points
			k1CalibratedArray[calibrationStep] = k1Big;
			k2MaxCalibratedArray[calibrationStep] = k2MaxBig;
			calibrationStep++;	
					
            if (calibrationStep < nCalibrationDataPoints)
            {
                updateTargetDataPoint();
				resetPistonVelocity();
            }
            else
            {
				computeVariations();
				updateVariationArrays();
				
                std::cout << "Compression cycles finished. Quitting...\n";
                //meanK1Big = computeMean(k1CalibratedArray, nCalibrationDataPoints);
                //pointerToDeviationArray[nCurrentRoutineRun] = computeDeviation(k1CalibratedArray, nCalibrationDataPoints);
                //pointerToDisplacementArray[nCurrentRoutineRun] = computeDisplacement(pointerToInitialK1Array[nCurrentRoutineRun], k1CalibratedArray, nCalibrationDataPoints);
                
                makeKdataFile();
                if (nCurrentRoutineRun < nMaximumRoutineRuns) updateInitialKArray();
                
                setTimeMax(getTime() + 10.0*getTimeStep());
                stage = 10;
            }
        }
   }  
   
   void computeVariations()
   {
	   meanK1Big = 0.0;
	   meanRelativeDeviationK1Big = 0.0;
	   meanRelativeDisplacementeK1Big = 0.0;
	   
	   for (int i = 0; i < nCalibrationDataPoints; i++) meanK1Big += k1CalibratedArray[i];
	   meanK1Big /= nCalibrationDataPoints;
	   
	   for (int i = 0; i < nCalibrationDataPoints; i++) meanRelativeDeviationK1Big += fabs(k1CalibratedArray[i] - meanK1Big);
	   meanRelativeDeviationK1Big /= nCalibrationDataPoints*meanK1Big;
	   
	   meanRelativeDisplacementeK1Big = fabs(meanK1Big - pointerToInitialK1Array[nCurrentRoutineRun])/pointerToInitialK1Array[nCurrentRoutineRun];
   }
   
   void updateVariationArrays()
   {
	   pointerToDeviationArray[nCurrentRoutineRun] = meanRelativeDeviationK1Big;
       pointerToDisplacementArray[nCurrentRoutineRun] = meanRelativeDisplacementeK1Big;
   }
   
   void updateInitialKArray()
   {
	   tunePhi();
	   pointerToInitialK1Array[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialK1Array[nCurrentRoutineRun] + meanK1Big);
	   pointerToInitialKCArray[nCurrentRoutineRun + 1] = pointerToInitialKCArray[nCurrentRoutineRun]; // this will be changed during the shear test
	   pointerToInitialPhiArray[nCurrentRoutineRun + 1] = phiBig;
   }

	// creates the kdata output file
    void makeKdataFile()
    {
        std::ostringstream kdataName;
        std::cout.unsetf(std::ios::floatfield);
        kdataName << getName() << ".kdata";
        
        kdataFile.open(kdataName.str(), std::ios::out);
        kdataFile << "k1_init = " << pointerToInitialK1Array[nCurrentRoutineRun] << std::endl;
        kdataFile << "k2_max = " << k2MaxBig << std::endl;
        kdataFile << "kC = " << pointerToInitialKCArray[nCurrentRoutineRun] << std::endl;
        kdataFile << "phi = " << pointerToInitialPhiArray[nCurrentRoutineRun] << std::endl;
        kdataFile << "|k1| = " << meanK1Big << std::endl;
        kdataFile << "MD(k1)/|k1| = " << pointerToDeviationArray[nCurrentRoutineRun] << std::endl;
        kdataFile << "||k1| - k1_init|/k1_init = " << pointerToDisplacementArray[nCurrentRoutineRun] << std::endl;
        kdataFile << "height \t pressure  \t k1" << std::endl;
        for (int i = 0; i < nCalibrationDataPoints; i++)
        {
            kdataFile << heightSteps[i] << "\t" << pressureSteps[i] << "\t" << k1CalibratedArray[i] << std::endl;
        }
        
        kdataFile.close();
    }
      
    // creates the data output file and writes the first row
    void makeCdatFile()
    {
        std::ostringstream cdatName;
        std::cout.unsetf(std::ios::floatfield);
        cdatName << getName() << ".cdat";
        
        cdatFile.open(cdatName.str(), std::ios::out);
        cdatFile << "time \t Eel \t Ekin/Eel \t h \t v \t P \t Pbase \t Pmax \t P/Pmax \t |(P - Pmax)/Pmax| \t dTotMean \t dTotMax \t dPistonMean \t dPistonMax \t dBaseMean \t dBaseMax \t k1 \t k2 \t kC \t phi" << std::endl;
    }
    
    // writes the compression data to the output file
    void writeToCdatFile()
    {
        cdatFile <<
        getTime() << "   " <<
        getElasticEnergy() << "   " <<
        getKineticEnergy()/getElasticEnergy() << "   " <<
        pistonHeight << "   " <<
        pistonVelocity << "   " <<
        pistonPressure << "   " <<
        basePressure << "   " <<
        targetPressure << "   " <<
        (pistonPressure + 0.01)/(targetPressure + 0.01) << "   " <<
        fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << "   " <<
        meanTotalRelativeOverlap << "   " <<
        maxTotalRelativeOverlap << "   " <<
        meanPistonRelativeOverlap << "   " <<
        maxPistonRelativeOverlap << "   " <<
        meanBaseRelativeOverlap << "   " <<
        maxBaseRelativeOverlap << "   " <<
        k1Big << "   " <<
        k2MaxBig << "   " <<
        kCBig << "   " <<
        pointerToInitialPhiArray[nCurrentRoutineRun] << "   " <<
        std::endl;
    }
    

    //  ----- GLOBAL FUNCTIONS -----
    void printTime() const override
    {
		if (stage < 4)
		{
			std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << 
			", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
			std::cout.flush();
		}
		else
		{
			std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << 
			", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight << 
			", v = " << pistonVelocity << std::endl <<
			"P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << " , P_target = " << targetPressure << 
			", P/P_target = " << (pistonPressure + 0.01)/(targetPressure + 0.01) << ", |(P - P_target)/P_target| = " << fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << std::endl <<
			"cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap << 
			", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap << 
			", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
			"k1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() << 
			", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kC = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl << std::endl;
			std::cout.flush();
		}
    }
  
  
	// VARIABLES   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
	// particles
	double radiusBig, radiusSmall;
	double sizeDispersityBig, sizeDispersitySmall;
    double densityBig, densitySmall;
    double volumeBig, volumeSmall;
    double massBig, massSmall;
    double smallToBigMassRatio;
    double nBig, nSmall;
    int nBigInserted, nSmallInserted;
    
    // powder bed
    double particleBedHeight;
    double particleBedPackingFraction;
    
    // geometry
    double casingRadius, casingHeight;
    double densityWall;
    InfiniteWall base, top;
    AxisymmetricIntersectionOfWalls casing;
    
    // piston
    double pistonHeight;
    double pistonVelocity;
    double pistonVelocityScalingFactor;
    double pistonMaximumVelocity;
    double pistonForce;
    double pistonPressure;
    double pistonAngularVelocity;
    double pistonTorque;
    InfiniteWall compressionPiston;
    InfiniteWall* pistonPointer;
    
    // friction 
    double bigWallSlidingFrictionCoeff, smallWallSlidingFrictionCoeff;
    double bigWallRollingFrictionCoeff, smallWallRollingFrictionCoeff;
    double bigWallTorsionFrictionCoeff, smallWallTorsionFrictionCoeff;
    double bigBigSlidingFrictionCoeff, smallSmallSlidingFrictionCoeff, bigSmallSlidingFrictionCoeff;
    double bigBigRollingFrictionCoeff, smallSmallRollingFrictionCoeff, bigSmallRollingFrictionCoeff;
    double bigBigTorsionFrictionCoeff, smallSmallTorsionFrictionCoeff, bigSmallTorsionFrictionCoeff;
    
    // collision
    double k1Wall, k1Big, k1Small;
    double k2MaxBig, k2MaxSmall;
    double kCBig, kCSmall;
    double phiBig, phiSmall;
    double phiMin;
    double bigWallRestitutionCoeff, smallWallRestitutionCoeff;
    double bigBigRestitutionCoeff, smallSmallRestitutionCoeff, bigSmallRestitutionCoeff;
    
    // species
    LinearPlasticViscoelasticFrictionSpecies *speciesBig, *speciesSmall, *speciesWall;
    
    // output
    double cdatOutputTimeInterval;
    std::ofstream kdataFile;
    std::ofstream cdatFile;
    
    // data analysis
    double basePressure;
    double cylinderPressure;
    double meanCoordinationNumber;
    double meanTotalRelativeOverlap;
    double maxTotalRelativeOverlap;
    double meanPistonRelativeOverlap;
    double maxPistonRelativeOverlap;
    double meanBaseRelativeOverlap;
    double maxBaseRelativeOverlap;
    double meanCylinderRelativeOverlap;
    double maxCylinderRelativeOverlap;
    
    // calibration
    int nCalibrationDataPoints;
    double pauseAfterDataPointMet;
    double *k1CalibratedArray;
    double *kCCalibratedArray;
    double *k2MaxCalibratedArray;
    double *heightSteps;
    double *pressureSteps;
    double *torqueSteps;
    bool isHeightTargetMet;
    bool isPressureTargetMet;
    bool isTorqueTargetMet;
    int calibrationStep;
    double targetHeight;
    double targetPressure;
    double targetTorque;
    double previousTargetHeight;
    double maxHeightRelDeviation;
    double maxPressureRelDeviation;
    double nTimeStepsBetweenStiffnessAdjustments;
    double *pointerToInitialK1Array;
    double *pointerToInitialKCArray;
    double *pointerToInitialPhiArray;
    double *pointerToDeviationArray;
    double *pointerToDisplacementArray;
    double meanK1Big;
    double meanRelativeDeviationK1Big;
    double meanRelativeDisplacementeK1Big;
    double cConstant;
    double deltaHatMax;
    
    // global
    int stage;
    double energyRatioTolerance;
    double t0;
    bool waitForParticleSettling;
    int nMaximumRoutineRuns;
    int nCurrentRoutineRun;
    bool isPistonRotatingFlag;
    bool isOnlySettlingFlag;
    bool isPhiDynamicallyTunedFlag;
};


int main(int argc, char *argv[])
{
	// GLOBAL FLAGS
	bool settleOnly = false;
	bool rotatePiston = false;
	bool tunePhiDynamically = true;
	
	// TIME STEP
	double timeStep = 1.0e-6;
	
	// SETUP PARAMETERS
	double particleBedHeight = 0.019;
	double particleBedPackingFraction = 0.60;
	double casingRadius = 0.0125;
	double casingHeight = 2.0*particleBedHeight;
	double densityWalls = 5000.0;
	
	double bigToCasingSizeRatio = 40.0;
	double bigToSmallSizeRatio = 1.0;
	double smallToBigMassRatio = 0.0;
	double sizeDispersionBigParticles = 0.1;
	double sizeDispersionSmallParticles = 0.1;
	double densityBigParticles = 1452.7;
	double densitySmallParticles = 1452.7;
	
	// PISTON PARAMETERS
	double pistonVelocityScalingFactor = 1.0;
	double pistonAngularVelocity = 0.1*constants::pi;
	
	// CALIBRATION DATA POINTS
	double pauseDuration = 0.005;
	//const int nDataPoints = 17;
    //double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886, 0.0183917, 0.0184044, 0.0184047, 0.0183935, 0.0183803, 0.0183798, 0.018384, 0.0183902, 0.0183958, 0.0184039, 0.0184891};
    //double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71, 7929.76, 6001.37, 5993.35, 8013.84, 9934.54, 9056.17, 8030.74, 7008.09, 5953.65, 4975.92, 0.0};
	//const int nDataPoints = 16;
    //double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886, 0.0183917, 0.0184044, 0.0184047, 0.0183935, 0.0183803, 0.0183798, 0.018384, 0.0183902, 0.0183958, 0.0184039};
    //double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71, 7929.76, 6001.37, 5993.35, 8013.84, 9934.54, 9056.17, 8030.74, 7008.09, 5953.65, 4975.92};
    //const int nDataPoints = 10;
    //double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886, 0.0183917, 0.0184044, 0.0184047, 0.0183935};
    //double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71, 7929.76, 6001.37, 5993.35, 8013.84};
    const int nDataPoints = 6;
    double heightArray[nDataPoints] = {0.0185328, 0.0184762, 0.0184496, 0.0184252, 0.0184058, 0.0183886};
    double pressureArray[nDataPoints] = {4790.42, 5763.35, 6749.81, 7771.34, 8751.82, 9717.71};
    //const int nDataPoints = 2;
    //double heightArray[nDataPoints] = {0.0185328, 0.0184762};
    //double pressureArray[nDataPoints] = {4790.42, 5763.35};
    //double torqueArray[nDataPoints] = {0.002, 0.003};
    
    double deltaHatMax = (particleBedHeight - 0.0183886)/particleBedHeight;
    double cConstant = (particleBedHeight - 0.0183886)/(0.0184891 - 0.0183886);
    
    // INTERACTION PARAMETERS
	double k1Big = 400.0;
	double k2MaxBig = 5000.0;
	double kCBig = 40.0;
	double phiBig = pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax/(cConstant - 1.0);
	
	double k1Small = 400.0;
	double k2MaxSmall = 5000.0;
	double kCSmall = 40.0;
	double phiSmall = 0.05;
	
	double k1Wall = 5000.0;
	
	double eBigWall = 0.7;
	double eSmallWall = 0.7;
	double eBigBig = 0.5;
	double eSmallSmall = 0.5;
	double eBigSmall = 0.5;
	
	// FRICTION PARAMETERS
	double muSBigBig = 0.16;
	double muRBigBig = 0.05;
	double muSBigWall = 0.10;
	double muRBigWall = 0.01;
	
    // CALIBRATION PARAMETERS
    double maximumStiffnessRelativeDeviation = 0.05;
    double maximumHeightRelativeDeviation = 1.0e-4;
    double maximumPressureRelativeDeviation = 1.0e-4;
    double timeStepsBetweenStiffnessTunings = 10.0;
    
    int runNumber = 0; 
    const int nMaximumRoutineCycles = 10;
    double initialK1Array[nMaximumRoutineCycles];
    double initialKCArray[nMaximumRoutineCycles];
    double initialPhiArray[nMaximumRoutineCycles];
    double deviationArray[nMaximumRoutineCycles];
    double displacementArray[nMaximumRoutineCycles];
  
	for (int i = 0; i < nMaximumRoutineCycles; i++) 
	{
		initialK1Array[i] = 0.0;
		initialKCArray[i] = 0.0;
		initialPhiArray[i] = 0.0;
		deviationArray[i] = 1.0;
		displacementArray[i] = 1.0;
	}
	
	initialK1Array[0] = k1Big;
	initialKCArray[0] = kCBig;
	initialPhiArray[0] = phiBig;
	
	// STRING VARIABLES
    std::ostringstream name;
    name.str("");
    name.clear();
    
    if (argc > 1) // THE CALIBRATION ROUTINE
    {
		settleOnly = false;
		
		do 
		{
			// --- TESTING FUNCTIONS ---
			std::cout << std::endl << std::endl;
			std::cout << "C constant is " << cConstant << std::endl;
			std::cout << "DeltaHat_Max is " << deltaHatMax << std::endl << std::endl;
			
			std::cout << "BEFORE EXECUTION" << std::endl;
			std::cout << "K1_init \t KC_init \t phi_init \t devK1 \t dispK1" << std::endl;
			for (int i = 0; i < nMaximumRoutineCycles; i++) std::cout << initialK1Array[i] << "\t\t" << initialKCArray[i] << "\t\t" << initialPhiArray[i] << "\t\t" << deviationArray[i] << "\t\t" << displacementArray[i] << std::endl;
			std::cout << std::endl << std::endl;
			// --- --- --- --- --- ---
    
			CompressionTest_calibrationRoutine problem;
			
			problem.readArguments(argc, argv);
			
			// INITIALIZATION
			problem.setTimeStep(timeStep);
			problem.setTimeMax(10.0);
			problem.setGravity(Vec3D(0.00,0.00,-9.81));
			problem.setSystemDimensions(3);
			//problem.setVerbose(true);
			problem.setCdatOutputTimeInterval(0.001);
			problem.setEnergyRatioTolerance(1.0e-4);
			problem.setSaveCount(0.01/problem.getTimeStep());
			
			problem.setCasingProperties(casingRadius, casingHeight, densityWalls);
			problem.setPowderBedProperties(particleBedHeight, particleBedPackingFraction);
			problem.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);
			
			problem.setParticleWallSlidingFrictionCoefficients(muSBigWall, muSBigWall);
			problem.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
			problem.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);
			
			problem.setParticleParticleSlidingFrictionCoefficients(muSBigBig, muSBigBig, muSBigBig);
			problem.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
			problem.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);
			
			problem.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
			problem.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
			problem.setBigParticlePlasticProperties(initialK1Array[runNumber], k2MaxBig, initialKCArray[runNumber], initialPhiArray[runNumber]);
			problem.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);
			
			problem.setCalibrationDataPoints(nDataPoints, heightArray, pressureArray, pauseDuration);
			problem.setCConstantAndDeltaHatMax(cConstant, deltaHatMax);
			problem.setCalibrationMaximumRelativeDeviations(maximumHeightRelativeDeviation, maximumPressureRelativeDeviation);
			problem.setCalibrationRoutineParameters(runNumber, nMaximumRoutineCycles, initialK1Array, initialKCArray, initialPhiArray, deviationArray, displacementArray);
			problem.setNTimeStepsBetweenStiffnesAdjustments(timeStepsBetweenStiffnessTunings);
			
			problem.setGlobalFlags(settleOnly, tunePhiDynamically);
			problem.setPistonVelocityScalingFactor(pistonVelocityScalingFactor);
			
			// NAME SETTING
			name.str("");
			name.clear();
			std::cout.unsetf(std::ios::floatfield);
			name << "CompressionTest_NEWROUTINE_CARRIER_Nrun_" << runNumber << "_rRatio_" << std::fixed << std::setprecision(0) <<  bigToCasingSizeRatio << 
			"_k1_" << initialK1Array[runNumber] << "_k2max_" << k2MaxBig << "_kc_" << initialKCArray[runNumber] << "_phi_" << std::fixed << std::setprecision(4) << phiBig << 
			"_muSBB_" << std::fixed << std::setprecision(2) << muSBigBig << "_muRBB_" << muRBigBig << "_muSBW_" << muSBigWall << "_muRBW_" << muRBigWall << "___TESTING_phiAutoTuning";
			problem.setName(name.str());
			
			std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl << std::endl;
			
			problem.solve();
			runNumber++;
			
			// --- TESTING FUNCTIONS ---
			std::cout << std::endl << std::endl;
			std::cout << "AFTER EXECUTION" << std::endl;
			std::cout << "K1_init \t KC_init \t phi_init \t devK1 \t dispK1" << std::endl;
			for (int i = 0; i < nMaximumRoutineCycles; i++) std::cout << initialK1Array[i] << "\t" << initialKCArray[i] << "\t" << initialPhiArray[i] << "\t" << deviationArray[i] << "\t" << displacementArray[i] << std::endl;
			std::cout << std::endl << std::endl;
			// --- --- --- --- --- ---
		}
		while (runNumber < nMaximumRoutineCycles && !(displacementArray[runNumber - 1] < maximumStiffnessRelativeDeviation && deviationArray[runNumber - 1] < maximumStiffnessRelativeDeviation));
    }
    else // THE SETTLING ONLY
    {
		settleOnly = true;
		
		CompressionTest_calibrationRoutine problem;
		
		// INITIALIZATION
		problem.setTimeStep(timeStep);
		problem.setTimeMax(10.0);
		problem.setGravity(Vec3D(0.00,0.00,-9.81));
		problem.setSystemDimensions(3);
		//problem.setVerbose(true);
		problem.setCdatOutputTimeInterval(0.005);
		problem.setEnergyRatioTolerance(1.0e-4);
		problem.setSaveCount(0.01/problem.getTimeStep());
		
		problem.setCasingProperties(casingRadius, casingHeight, densityWalls);
		problem.setPowderBedProperties(particleBedHeight, particleBedPackingFraction);
		problem.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);
		
		problem.setParticleWallSlidingFrictionCoefficients(0.3, 0.3);
		problem.setParticleWallRollingFrictionCoefficients(0.01, 0.01);
		problem.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);
		
		problem.setParticleParticleSlidingFrictionCoefficients(muSBigBig, 0.0, 0.0);
		problem.setParticleParticleRollingFrictionCoefficients(muRBigBig, 0.0, 0.0);
		problem.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);
		
		problem.setWallStiffnessAndRestitutionCoefficients(k1Wall, 0.7, 0.7);
		problem.setParticleParticleRestitutionCoefficients(0.5, 0.3, 0.4);
		problem.setBigParticlePlasticProperties(initialK1Array[runNumber], k2MaxBig, initialKCArray[runNumber], initialPhiArray[runNumber]);
		problem.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);
		
		problem.setGlobalFlags(settleOnly, tunePhiDynamically);
		
		// NAME SETTING
		std::cout.unsetf(std::ios::floatfield);
        name << "CompressionTest_SETTLEDPARTICLEBED_rRatio_" << std::fixed << std::setprecision(0) <<  bigToCasingSizeRatio << 
        "_k1_" << initialK1Array[runNumber] << "_k2max_" << k2MaxBig << "_kc_" << initialKCArray[runNumber] << "_phi_" << std::fixed << std::setprecision(2) << phiBig <<
        "_muSBB_" << muSBigBig << "_muRBB_" << muRBigBig << "___TESTING";
        problem.setName(name.str());
        
		std::cout << "SETTLING OF THE PARTICLE BED" << std::endl;

		problem.solve();
	}
	return 0;
}


// this is the calibration routine where the calibration of k2max was included
/*
void runDataPointsCalibration()
	{

		  
		// reverse and halve the velocity if the target height was passed
        if (!isHeightTargetMet && ((pistonVelocity < 0.0 && previousTargetHeight < targetHeight) || (pistonVelocity > 0.0 && previousTargetHeight > targetHeight))) pistonVelocity *= -0.5;

		// upon reaching the target height update the bool and sets the velocity to 0
        if (!isHeightTargetMet && !isPressureTargetMet)
        {
			if (fabs((pistonHeight - targetHeight)/targetHeight) > maxHeightRelDeviation) movePiston();
			else
			{
				std::cout << "Height target reached.\n";
				std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
                t0 = getTime();
                    
				pistonVelocity = 0.0;
				isHeightTargetMet = true;
			}
        }

		// when the height is fine deal with the pressure
        if (isHeightTargetMet && !isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
        {
			if (previousTargetHeight > targetHeight) // calibration of k1
			{
				if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviationK1)
				{
					// the stiffness is changed every n time steps
					if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneK1();
				}
				else
				{
					std::cout << "Pressure target reached, now chillin' for " << pauseAfterDataPointMet << " seconds\n";
					t0 = getTime();
					isPressureTargetMet = true;
				}
			}
			else // calibration of k2Max
			{
				if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviationK1)
				{
					// the stiffness is changed every n time steps
					if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneK2Max();
				}
				else
				{
					std::cout << "Pressure target reached, now chillin' for " << pauseAfterDataPointMet << " seconds\n";
					t0 = getTime();
					isPressureTargetMet = true;
				}
			}
		}
  
        // switch to next calibration point if the test of the previous is passed
        if (isHeightTargetMet && isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
        {
			// updates the mean value of the stiffnesses at the calibrated data points
			k1CalibratedArray[calibrationStep] = k1Big;
			k2MaxCalibratedArray[calibrationStep] = k2MaxBig;
			calibrationStep++;	
					
            if (calibrationStep < nCalibrationDataPoints)
            {
                updateTargetDataPoint();
				resetPistonVelocity();
            }
            else
            {
                std::cout << "Compression cycles finished. Quitting...\n";
                meanK1Big = computeMean(k1CalibratedArray, nCalibrationDataPoints);
                pointerToDeviationArray[nCurrentRoutineRun] = computeDeviation(k1CalibratedArray, nCalibrationDataPoints, meanK1Big);
                pointerToDisplacementArray[nCurrentRoutineRun] = computeDisplacement(pointerToInitialK1Array[nCurrentRoutineRun], meanK1Big);
                
                //std::cout << std::endl << nCalibrationDataPoints << "   " << nCurrentRoutineRun << std::endl;
                //std::cout << "DEBUG - mean: " << meanK1Big << ", deviation: " << pointerToDeviationArray[nCurrentRoutineRun] << ", displacement: " << pointerToDisplacementArray[nCurrentRoutineRun] << std::endl;
                
                for (int i = 0; i < nCalibrationDataPoints; i++) std::cout << k1CalibratedArray[i] << "   " << k2MaxCalibratedArray[i] << std::endl;
                
                if (nCurrentRoutineRun < nMaximumRoutineRuns) updateInitialKArray();
                
                setTimeMax(getTime() + 10.0*getTimeStep());
                stage = 10;
            }
        }
   }
   */

/*
 * // when the height is fine deal with the pressure
        if (isHeightTargetMet && !isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
        {
			if (previousTargetHeight > targetHeight) // calibration of k1
			{
				if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviationK1)
				{
					// the stiffness is changed every n time steps
					if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneK1();
				}
				else
				{
					std::cout << "Pressure target reached.\n";
					std::cout << pistonPressure << "\t" << targetPressure << "\t" << maxPressureRelDeviationK1 << "\t" << fabs((pistonPressure - targetPressure)/targetPressure) << std::endl;
					
					// updates the mean value of the stiffnesses at the calibrated data points
					meanK1Big += k1Big;
					meanKCBig += kCBig;
					
					std::cout << std::endl << "MEAN K1: " << meanK1Big/(nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun] + 1);
					std::cout << "MEAN KC: " << meanKCBig/(nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun] + 1) << std::endl;
					
					if (isStiffnessAdjusted) // if the stiffness was adjusted and we have runs left reboot the routine and update the initial stiffness
					{
						if (nCurrentRoutineRun < nMaximumRoutineRuns - 1)
						{
							pointerToInitialK1Array[nCurrentRoutineRun + 1] = 0.5*(k1Big + pointerToInitialK1Array[nCurrentRoutineRun]);
							pointerToInitialKCArray[nCurrentRoutineRun + 1] = pointerToInitialKCArray[nCurrentRoutineRun];
							
							std::cout << "Resetting the routine...\n";
							setTimeMax(getTime() + 10.0*getTimeStep());
							stage = 10;
						}
						else
						{
							std::cout << "Stiffness was adjusted but there are no runs left.\n";
							setTimeMax(getTime() + 10.0*getTimeStep());
							stage = 10;
						}
					}
					else // if the initial stiffness was fine go on with the actual run and update the pressure flag
					{
						std::cout << "The initial stiffness were fine.\n";
						std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
						t0 = getTime();
						
						nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun]++;
						isPressureTargetMet = true;
					}
				}
			}
			else // calibration of kC
			{
				if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviationKC)
				{
					// the stiffness is changed every n time steps
					if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) 
					{
						tuneKC();
						
						if (nCurrentRoutineRun < nMaximumRoutineRuns - 1)
						{
							pointerToInitialK1Array[nCurrentRoutineRun + 1] = pointerToInitialK1Array[nCurrentRoutineRun];
							pointerToInitialKCArray[nCurrentRoutineRun + 1] = kCBig;
							
							std::cout << "Resetting the routine...\n";
							setTimeMax(getTime() + 10.0*getTimeStep());
							stage = 10;
						}
						else
						{
							std::cout << "Stiffness was adjusted but there are no runs left.\n";
							setTimeMax(getTime() + 10.0*getTimeStep());
							stage = 10;
						}
					}
				}
				else // if the initial stiffness was fine go on with the actual run and update the pressure flag
				{
					std::cout << "The initial stiffness were fine.\n";
					std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
					t0 = getTime();
					
					nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun]++;
					isPressureTargetMet = true;
				}
			}
		}
 * */

// old part of runDataPointsCalibration() without the kC calibration
/*
 // when the height is fine deal with the pressure
        if (isHeightTargetMet && !isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
        {
			// if the pressure is not met, calibrate the spring stiffness (either k1 or kC)
            if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
            {
                // the stiffness is changed every n time steps
				if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep())
				{
					// on compression changes k1, otherwise changes kC
					if (previousTargetHeight > targetHeight)
					{
						tuneK1();
					}
					else
					{
						if (pistonPressure < targetPressure) kCBig *= 0.95;
						else kCBig *= 1.05;
						
						isStiffnessAdjusted = true;
						
						if (kCBig < 0.0) 
						{
							std::cout << "\n\nFATAL ERROR: KC_BIG < 0. QUITTING.";
							exit(0);
						}
					}
				}
            }
            else // when the target pressure is reached check what to do
            {
                std::cout << "Pressure target reached.\n";
                std::cout << pistonPressure << "\t" << targetPressure << "\t" << maxPressureRelDeviation << "\t" << fabs((pistonPressure - targetPressure)/targetPressure) << std::endl;
                
                // updates the mean value of the stiffnesses at the calibrated data points
                meanK1Big += k1Big;
                meanKCBig += kCBig;
                
                std::cout << std::endl << "MEAN K1: " << meanK1Big/(nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun] + 1);
                std::cout << "MEAN KC: " << meanKCBig/(nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun] + 1) << std::endl;
                
                if (isStiffnessAdjusted) // if the stiffness was adjusted and we have runs left reboot the routine and update the initial stiffness
                {
					if (nCurrentRoutineRun < nMaximumRoutineRuns - 1)
					{
						pointerToInitialK1Array[nCurrentRoutineRun + 1] = 0.5*(k1Big + pointerToInitialK1Array[nCurrentRoutineRun]);
						pointerToInitialKCArray[nCurrentRoutineRun + 1] = 0.5*(kCBig + pointerToInitialKCArray[nCurrentRoutineRun]);
						
						std::cout << "Resetting the routine...\n";
                        setTimeMax(getTime() + 10.0*getTimeStep());
                        stage = 10;
					}
					else
					{
						std::cout << "Stiffness was adjusted but there are no runs left.\n";
                        setTimeMax(getTime() + 10.0*getTimeStep());
                        stage = 10;
					}
				}
				else // if the initial stiffness was fine go on with the actual run and update the pressure flag
				{
					std::cout << "The initial stiffness were fine.\n";
                    std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
                    t0 = getTime();
                    
                    nCalibratedDataPointsPerRunPointer[nCurrentRoutineRun]++;
                    isPressureTargetMet = true;
				}                
            }
        }
  */


















