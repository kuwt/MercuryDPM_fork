#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

/*
 * ToDo:
 * - fix the restart file (lambda functions? check other codes!)
 * - implement data analysis
 * - fix teh drum filling ratio function 
 */
 

class TaperedDrum : public Mercury3D
{
	
	private:
	
	void setupInitialConditions() override
    {
		stage = 1;
		nBigInserted = 0;
		nSmallInserted = 0;
		t0 = getTime();
		
		waitForParticleSettling = false;
		
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
		insertParticles();
		stage++;
    }
    
    void actionsOnRestart() override
    {
		std::cout << "FILE HAS BEEN RESTARTED" << std::endl << std::endl;
		stage = 3;
	}
    
    void actionsAfterTimeStep() override
    {	
		if (stage == 2) insertParticles();
		if (stage == 3) 
		{
			rotateDrum();
		}
		
		if (stage >= 4) 
		{
			//if (fmod(getTime(),0.01) < getTimeStep()) writeToKdataFile();
		}
		
		if (stage == 4)
		{
			//movePiston();
			
			//makeDataAnalysis();
			//runDataPointsCalibration();
			

		}
	}
    
    void actionsAfterSolve() override
    {
		//if (!isOnlySettlingFlag)
		//{
			//delete [] heightSteps;
			//delete [] pressureSteps;
			
			//kdataFile.close();
		//}
	}
    
    public:
    
    
    // FUNCTIONS CALLED IN MAIN ----------------------------------------
    
    void setCdatOutputTimeInterval(double dt)
    {
        cdatOutputTimeInterval = dt;
    }
    
    void setEnergyRatioTolerance(double eRatio)
    {
        energyRatioTolerance = eRatio;
    }
    
    void setCasingProperties(double r, double R, double h, double density)
    {
		radiusMinor = r;
		radiusMajor = R;
		height = h;
		
		setXMin(-radiusMajor);
        setYMin(-radiusMajor);
        setZMin(0.0);
        
        setXMax(radiusMajor);
        setYMax(radiusMajor);
        setZMax(height);		
	}
	
	void setOperatingParameters(double fr, double ratio, double pf)
    {
		froudeNumber = fr;
		drumAngularVelocity = sqrt(2.0*9.81*froudeNumber/(radiusMinor + radiusMajor));
		drumFillingRatio = ratio;
		particleBedPackingFraction = pf;
	}
	
    void setParticleProperties(double rS, double rB, double dB, double dS, double rhoB, double rhoS, double ratio)
    {
		radiusSmall = rS;
		radiusBig = rB;
		sizeDispersityBig = dB;
		sizeDispersitySmall = dS;
		densityBig = rhoB;
		densitySmall = rhoS;
		bigToSmallVolumeRatio = ratio;
		
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

	void setGlobalFlags(bool settlingFlag)
	{
		isOnlySettlingFlag = settlingFlag;
	}
	
	
	
	// FUNCTIONS CALLED IN THE CLASS -----------------------------------
	
	void setSpecies()
    {
        speciesHandler.clear();

        // BIG-BIG
        speciesBig = new LinearViscoelasticFrictionSpecies;
        speciesBig -> setDensity(densityBig);
        speciesBig -> setStiffnessAndRestitutionCoefficient(k1Big, bigBigRestitutionCoeff, massBig);
        //speciesBig -> setUnloadingStiffnessMax(k2MaxBig);
        //speciesBig -> setCohesionStiffness(kCBig);
        //speciesBig -> setPenetrationDepthMax(phiBig);

        speciesBig -> setSlidingFrictionCoefficient(bigBigSlidingFrictionCoeff);
        speciesBig -> setSlidingStiffness(speciesBig -> getStiffness()*2.0/7.0);
        speciesBig -> setSlidingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setRollingFrictionCoefficient(bigBigRollingFrictionCoeff);
        speciesBig -> setRollingStiffness(speciesBig -> getStiffness()*2.0/7.0);
        speciesBig -> setRollingDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesBig -> setTorsionFrictionCoefficient(bigBigTorsionFrictionCoeff);
        speciesBig -> setTorsionStiffness(speciesBig -> getStiffness()*2.0/7.0);
        speciesBig -> setTorsionDissipation(speciesBig -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesBig);
        
        // SMALL-SMALL
        speciesSmall = new LinearViscoelasticFrictionSpecies;
        speciesSmall -> setDensity(densitySmall);
        speciesSmall -> setStiffnessAndRestitutionCoefficient(k1Small, smallSmallRestitutionCoeff, massSmall);
        //speciesSmall -> setUnloadingStiffnessMax(k2MaxSmall);
        //speciesSmall -> setCohesionStiffness(kCSmall);
        //speciesSmall -> setPenetrationDepthMax(phiSmall);

        speciesSmall -> setSlidingFrictionCoefficient(smallSmallSlidingFrictionCoeff);
        speciesSmall -> setSlidingStiffness(speciesSmall -> getStiffness()*2.0/7.0);
        speciesSmall -> setSlidingDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesSmall -> setRollingFrictionCoefficient(smallSmallRollingFrictionCoeff);
        speciesSmall -> setRollingStiffness(speciesSmall -> getStiffness()*2.0/7.0);
        speciesSmall -> setRollingDissipation(speciesSmall -> getDissipation()*2.0/7.0);        
        speciesSmall -> setTorsionFrictionCoefficient(smallSmallTorsionFrictionCoeff);
        speciesSmall -> setTorsionStiffness(speciesSmall -> getStiffness()*2.0/7.0);
        speciesSmall -> setTorsionDissipation(speciesSmall -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesSmall);
        
        // WALL-WALL
        speciesWall = new LinearViscoelasticFrictionSpecies;
        speciesWall -> setDensity(densityWall);
        speciesWall -> setStiffnessAndRestitutionCoefficient(k1Wall, 1.0, massSmall);
        //speciesWall -> setUnloadingStiffnessMax(k1Wall);
        //speciesWall -> setCohesionStiffness(0.0);
        //speciesWall -> setPenetrationDepthMax(0.0);
        
        speciesWall -> setSlidingFrictionCoefficient(0.0);
        speciesWall -> setSlidingStiffness(speciesWall -> getStiffness()*2.0/7.0);
        speciesWall -> setSlidingDissipation(speciesWall -> getDissipation()*2.0/7.0);
        speciesWall -> setRollingFrictionCoefficient(0.0);
        speciesWall -> setRollingStiffness(speciesWall -> getStiffness()*2.0/7.0);
        speciesWall -> setRollingDissipation(speciesWall -> getDissipation()*2.0/7.0);
        speciesWall -> setTorsionFrictionCoefficient(0.0);
        speciesWall -> setTorsionStiffness(speciesWall -> getStiffness()*2.0/7.0);
        speciesWall -> setTorsionDissipation(speciesWall -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesWall);
        
        // BIG-WALL
        //speciesMixedBigWall = speciesHandler.getMixedObject(speciesBig, speciesWall);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Wall), bigWallRestitutionCoeff, massBig);
        //speciesHandler.getMixedObject(speciesBig, speciesWall) -> setUnloadingStiffnessMax(k2MaxBig);
        //speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kCBig);
        //speciesHandler.getMixedObject(speciesBig, speciesWall) -> setPenetrationDepthMax(phiBig);

        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingFrictionCoefficient(bigWallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingFrictionCoefficient(bigWallRollingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionFrictionCoefficient(bigWallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);
        
        // SMALL-WALL
        //speciesMixedSmallWall = speciesHandler.getMixedObject(speciesSmall, speciesWall);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Small + k1Wall), smallWallRestitutionCoeff, massSmall);
        //speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setUnloadingStiffnessMax(k2MaxSmall);
        //speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setCohesionStiffness(kCSmall);
        //speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setPenetrationDepthMax(phiSmall);

        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingFrictionCoefficient(smallWallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingFrictionCoefficient(smallWallRollingFrictionCoeff);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionFrictionCoefficient(smallWallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);
        
        // BIG-SMALL
        //speciesMixedBigSmall = speciesHandler.getMixedObject(speciesBig, speciesSmall);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Small), bigSmallRestitutionCoeff, 0.5*(massBig + massSmall));
        //speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setUnloadingStiffnessMax(0.5*(k2MaxBig + k2MaxSmall));
        //speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kCBig + kCSmall));
        //speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setPenetrationDepthMax(0.5*(phiBig + phiSmall));
        
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingFrictionCoefficient(bigSmallSlidingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingFrictionCoefficient(bigSmallRollingFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionFrictionCoefficient(bigSmallTorsionFrictionCoeff);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getStiffness()*2.0/7.0);
        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
		
		//std::cout << "BIG-BIG stiffness and dissipation coefficients: " << speciesBig -> getLoadingStiffness() << " " << speciesBig -> getDissipation() << "\n";
		//std::cout << "BIG-BIG friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
		//std::cout << "BIG-BIG tangential stiffnesses: " << speciesBig -> getSlidingStiffness() << " " << speciesBig -> getRollingStiffness() << " " << speciesBig -> getTorsionStiffness() << "\n";
		//std::cout << "BIG-BIG tangential dissipation coefficients: " << speciesBig -> getSlidingDissipation() << " " << speciesBig -> getRollingDissipation() << " " << speciesBig -> getTorsionDissipation() << "\n";
		//std::cout << "BIG-BIG collision time: " << std::setprecision(4) << speciesBig -> getCollisionTime(massBig) << "\n\n";
		
		//std::cout << "SMALL-SMALL stiffness and dissipation coefficients: " << speciesSmall -> getLoadingStiffness() << " " << speciesSmall -> getDissipation() << "\n";
		//std::cout << "SMALL-SMALL friction coefficients: " << smallSmallSlidingFrictionCoeff << " " << smallSmallRollingFrictionCoeff << " " << smallSmallTorsionFrictionCoeff << "\n";
		//std::cout << "SMALL-SMALL tangential stiffnesses: " << speciesSmall -> getSlidingStiffness() << " " << speciesSmall -> getRollingStiffness() << " " << speciesSmall -> getTorsionStiffness() << "\n";
		//std::cout << "SMALL-SMALL tangential dissipation coefficients: " << speciesSmall -> getSlidingDissipation() << " " << speciesSmall -> getRollingDissipation() << " " << speciesSmall -> getTorsionDissipation() << "\n";
		//std::cout << "SMALL-SMALL collision time: " << std::setprecision(4) << speciesSmall -> getCollisionTime(massSmall) << "\n\n";
		
		//std::cout << "BIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation() << "\n";
		//std::cout << "BIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
		//std::cout << "BIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionStiffness() << "\n";
		//std::cout << "BIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionDissipation() << "\n";
		//std::cout << "BIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(massBig) << "\n\n";
		
		//std::cout << "SMALL-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation() << "\n";
		//std::cout << "SMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
		//std::cout << "SMALL-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getTorsionStiffness() << "\n";
		//std::cout << "SMALL-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getTorsionDissipation() << "\n";
		//std::cout << "SMALL-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getCollisionTime(massSmall) << "\n\n";
		
		//std::cout << "BIG-SMALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation() << "\n";
		//std::cout << "BIG-SMALL friction coefficients: " << bigSmallSlidingFrictionCoeff << " " << bigSmallRollingFrictionCoeff << " " << bigSmallTorsionFrictionCoeff << "\n";
		//std::cout << "BIG-SMALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getTorsionStiffness() << "\n";
		//std::cout << "BIG-SMALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getTorsionDissipation() << "\n";
		//std::cout << "BIG-SMALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getCollisionTime(0.5*(massBig + massSmall)) << "\n\n";
		
		std::cout << "tC_BB/dt: " << std::setprecision(4) << speciesBig -> getCollisionTime(massBig)/getTimeStep() << "\n";
		std::cout << "tC_SS/dt: " << std::setprecision(4) << speciesSmall -> getCollisionTime(massSmall)/getTimeStep() << "\n";
		std::cout << "tC_BW/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(massBig)/getTimeStep() << "\n";
		std::cout << "tC_SW/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getCollisionTime(massSmall)/getTimeStep() << "\n";
		std::cout << "tC_BS/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getCollisionTime(0.5*(massBig + massSmall))/getTimeStep() << "\n\n";
    }
    
	void computeNumberOfParticles()
    {
		double particleBedVolume = constants::pi*pow(radiusMinor,2.0)*height*drumFillingRatio;
		nBig = 0.5*particleBedVolume*0.66/volumeBig;
		nSmall = 8.0*nBig;
        //nSmall = (particleBedPackingFraction*particleBedVolume/volumeSmall)*(smallToBigMassRatio*densityBig/(densitySmall + smallToBigMassRatio*densityBig));
        //nBig = (particleBedPackingFraction*particleBedVolume - nSmall*volumeSmall)/volumeBig;
	}
	
	void makeCasing()
    {
        wallHandler.clear();
        
        base.setSpecies(speciesWall);
        base.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
        basePointer = wallHandler.copyAndAddObject(base);
        //wallHandler.copyAndAddObject(base);
        
        top.setSpecies(speciesWall);
        top.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
        topPointer = wallHandler.copyAndAddObject(top);
        //wallHandler.copyAndAddObject(top);
        
        AxisymmetricIntersectionOfWalls(sideWall);
        sideWall.setSpecies(speciesWall);
        sideWall.setPosition(Vec3D(0.0,0.0,0.0));
        sideWall.setOrientation(Vec3D(0.0,0.0,1.0));
        sideWall.addObject(Vec3D(height,0.0,radiusMinor-radiusMajor),Vec3D(radiusMinor,0.0,0.0));
        sideWallPointer = wallHandler.copyAndAddObject(sideWall);
        //wallHandler.copyAndAddObject(sideWall);
    }
	
	bool particleInsertionSuccessful(bool isBig)
	{
		int insertionFailCounter = 0;
		Vec3D particlePosition;
		BaseParticle p0;
		
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
			particlePosition.Z = 0.5*height + (0.5*height - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
			particlePosition.X = (radiusMinor + particlePosition.Z*(radiusMajor - radiusMinor)/height - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
			particlePosition.Y = sqrt(pow(radiusMinor + particlePosition.Z*(radiusMajor - radiusMinor)/height - 1.01*(p0.getRadius()),2.0) - pow(particlePosition.X,2.0))*random.getRandomNumber(-1.0,1.0);
			
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
		
		if (waitForParticleSettling && getKineticEnergy()/getElasticEnergy() < 10.0*energyRatioTolerance) waitForParticleSettling = false;
		if (!waitForParticleSettling && nBigInserted >= nBig && nSmallInserted >= nSmall && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance) 
		{
			std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
			stage++;
		}
	}
	
	void rotateDrum()
    {
		sideWallPointer -> setAngularVelocity(Vec3D(0.0,0.0,-drumAngularVelocity));
        sideWallPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
		basePointer -> setAngularVelocity(Vec3D(0.0,0.0,-drumAngularVelocity));
        basePointer -> setOrientation(Vec3D(0.0,0.0,1.0));
        topPointer -> setAngularVelocity(Vec3D(0.0,0.0,-drumAngularVelocity));
        topPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
    }

	
/*	
	//void makeDataAnalysis()
    //{
        //meanCoordinationNumber = 0.0;
        
        //for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        //{
            //meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        //}
        
        //meanCoordinationNumber /= particleHandler.getNumberOfObjects();
        
        //pistonForce = 0.0;
        //pistonPressure = 0.0;
        //pistonTorque = 0.0;
        //basePressure = 0.0;
        //meanTotalRelativeOverlap = 0.0;
        //maxTotalRelativeOverlap = 0.0;
        //meanPistonRelativeOverlap = 0.0;
        //maxPistonRelativeOverlap = 0.0;
        //meanBaseRelativeOverlap = 0.0;
        //maxBaseRelativeOverlap = 0.0;
        
        //int totalInteractionCounter = 0;
        //int pistonInteractionCounter = 0;
        //int baseInteractionCounter = 0;
        
        //for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        //{
            //meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            //if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            
            //if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
            //{
                //pistonPressure -= ((*i) -> getForce()).Z;
                
                //meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                //if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                //if (isPistonRotatingFlag) pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;
					
                //pistonInteractionCounter++;
            //}
            
            //// 0 is the index of the base of the casing
            //if ((*i) -> getI() -> getIndex() == base.getIndex())
            //{
                //basePressure += ((*i) -> getForce()).Z;
                
                //meanBaseRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                //if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxBaseRelativeOverlap) maxBaseRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
                
                //baseInteractionCounter++;
            //}
            
            //totalInteractionCounter++;
        //}
        
        //pistonForce = pistonPressure;
        //pistonPressure /= constants::pi*pow(casingRadius,2.0);
        //basePressure /= constants::pi*pow(casingRadius,2.0);
        //meanTotalRelativeOverlap /= totalInteractionCounter;
        //meanPistonRelativeOverlap /= pistonInteractionCounter;
        //meanBaseRelativeOverlap /= baseInteractionCounter;
    //}
    


	//// creates the data output file and writes the first row
    //void makeKdataFile()
    //{
        //std::ostringstream kdataName;
        //std::cout.unsetf(std::ios::floatfield);
        //kdataName << getName() << ".kdata";
        
        //kdataFile.open(kdataName.str(), std::ios::out);
        //kdataFile << "time \t k1_init \t k1 \t dk1_rel \t dk1_rel_MAX \t k1_flag \t P \t P_target" << std::endl;
    //}

	//void writeToKdataFile()
    //{
        //kdataFile <<
        //getTime() << "\t" <<
        //initialK1Big << "\t" <<
        //k1Big << "\t" <<
        //relativeStiffnessVariation << "\t" <<
        //maxStiffnessRelDeviation << "\t" <<
        //isStiffnessAdjusted << "\t" <<
        //pistonPressure << "\t" <<
        //targetPressure << std::endl;
    //}
    

    //  ----- GLOBAL FUNCTIONS -----
    //void printTime() const override
    //{
		//if (stage < 4)
		//{
			//std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << 
			//", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
			//std::cout.flush();
		//}
		//else
		//{
			//std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << 
			//", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight << 
			//", v = " << pistonVelocity << "\nP = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << " , P_target = " << targetPressure << 
			//", P/P_target = " << (pistonPressure + 0.01)/(targetPressure + 0.01) << ", |(P - P_target)/P_target| = " << fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << 
			//"\ncN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", meanTotOverlap = " << meanTotalRelativeOverlap << ", maxTotOverlap = " << maxTotalRelativeOverlap << 
			//", meanPistonOverlap = " << meanPistonRelativeOverlap << ", maxPistonOverlap = " << maxPistonRelativeOverlap << ", meanBaseOverlap = " << meanBaseRelativeOverlap << 
			//", maxBaseOverlap = " << maxBaseRelativeOverlap << "\nk1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() << 
			//", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kC = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl << std::endl;
			//std::cout.flush();
		//}
    //}
*/
  
  
	// VARIABLES -------------------------------------------------------
	// particles
	double radiusBig, radiusSmall;
	double sizeDispersityBig, sizeDispersitySmall;
    double densityBig, densitySmall;
    double volumeBig, volumeSmall;
    double massBig, massSmall;
    double bigToSmallVolumeRatio;
    double nBig, nSmall;
    int nBigInserted, nSmallInserted;
    
    // operating parameters
    double froudeNumber;
    double drumAngularVelocity;
	double drumFillingRatio;
	double particleBedPackingFraction;
    
    // geometry
    double radiusMinor, radiusMajor, height;
    double densityWall;
    InfiniteWall base, top;
    InfiniteWall *basePointer, *topPointer;
    AxisymmetricIntersectionOfWalls sideWall;
    AxisymmetricIntersectionOfWalls *sideWallPointer;
    
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
    double bigWallRestitutionCoeff, smallWallRestitutionCoeff;
    double bigBigRestitutionCoeff, smallSmallRestitutionCoeff, bigSmallRestitutionCoeff;
    
    // species
    LinearViscoelasticFrictionSpecies *speciesBig, *speciesSmall, *speciesWall;
    
    // output
    double cdatOutputTimeInterval;
    std::ofstream kdataFile;
    
    //// data analysis
    //double basePressure;
    //double meanCoordinationNumber;
    //double meanTotalRelativeOverlap;
    //double maxTotalRelativeOverlap;
    //double meanPistonRelativeOverlap;
    //double maxPistonRelativeOverlap;
    //double meanBaseRelativeOverlap;
    //double maxBaseRelativeOverlap;
    
    // global
    int stage;
    double energyRatioTolerance;
    double t0;
    bool waitForParticleSettling;
    bool isOnlySettlingFlag;
};


// ---------------------------------------------------------------------
int main(int argc, char *argv[])
{
	// GLOBAL FLAGS
	bool settleOnly = false;
	
	// TIME STEP
	double timeStep = 1.0e-5;
	
	// PROCESS PARAMETERS
	double froudeNumber = 0.1;
	double drumFillingRatio = 0.5;
	double particleBedPackingFraction = 0.66;
	
	// DRUM PARAMETERS
	double radiusMinor = 0.04;
	double radiusMajor = 0.06;
	double height = 0.1;
	double densityWalls = 2000.0;
	
	// PARTICLE PARAMETERS
	double particleRadiusSmall = 0.00125;
	double particleRadiusBig = 2.0*particleRadiusSmall;
	double sizeDispersionBigParticles = 0.1;
	double sizeDispersionSmallParticles = 0.1;
	double bigToSmallVolumeRatio = 0.5;
	double densityBigParticles = 2000.0;
	double densitySmallParticles = 2000.0;
	
	// INTERACTION PARAMETERS
	double k1Big = 1000.0;	// SWITCH TO RESTITUTION COEFFICIENT AND COLLISION TIME
	double k2MaxBig = 1000.0;
	double kCBig = 0.0;
	double phiBig = 0.05;
	
	double k1Small = 1000.0;
	double k2MaxSmall = 1000.0;
	double kCSmall = 0.0;
	double phiSmall = 0.05;
	
	double k1Wall = 1000.0;
	
	// FRICTION PARAMETERS
	double muSBigBig = 0.10;
	double muRBigBig = 0.05;
	double muSBigWall = 0.40;
	double muRBigWall = 0.05;

	// STRING VARIABLES
    std::ostringstream name;
    name.str("");
    name.clear();
    
   
	// -----------------------------------------------------------------
	TaperedDrum problem;
	
	problem.readArguments(argc, argv);
	
	// INITIALIZATION
	problem.setTimeStep(timeStep);
	problem.setTimeMax(20.0);
	problem.setGravity(Vec3D(0.00,-9.81,0.00));
	problem.setSystemDimensions(3);
	//problem.setVerbose(true);
	problem.setCdatOutputTimeInterval(0.005);
	problem.setEnergyRatioTolerance(1.0e-4);
	problem.setSaveCount(0.01/problem.getTimeStep());
	
	problem.setCasingProperties(radiusMinor, radiusMajor, height, densityWalls);
	problem.setOperatingParameters(froudeNumber, drumFillingRatio, particleBedPackingFraction);
	problem.setParticleProperties(particleRadiusSmall, particleRadiusBig, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, bigToSmallVolumeRatio);
	
	problem.setParticleWallSlidingFrictionCoefficients(muSBigWall, muSBigWall);
	problem.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
	problem.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);
	
	problem.setParticleParticleSlidingFrictionCoefficients(muSBigBig, muSBigBig, muSBigBig);
	problem.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
	problem.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);
	
	problem.setWallStiffnessAndRestitutionCoefficients(k1Wall, 0.7, 0.7); // SETUP VARIABLES
	problem.setParticleParticleRestitutionCoefficients(0.5, 0.3, 0.4); // SETUP VARIABLES
	problem.setBigParticlePlasticProperties(k1Big, k2MaxBig, kCBig, phiBig);
	problem.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);
	
	problem.setGlobalFlags(settleOnly);
	
	// NAME SETTING
	name.str("");
	name.clear();
	std::cout.unsetf(std::ios::floatfield);
	name << "TaperedDrum_LONG_TESTING";
	problem.setName(name.str());
	
	//std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl << std::endl;
	
	problem.solve();

	return 0;
}



















