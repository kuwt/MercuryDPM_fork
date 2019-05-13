#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/IntersectionOfWalls.h>
#include <Boundaries/PeriodicBoundary.h>
#include <math.h>
#include <fstream>
#include <chrono>
#include <ctime>

/*
NOTES:
12.07:
- reworked everything: used the CompressionTest template
- added automatic convergence cycle based on bisection method
*/

class CalibrationRoutine_AngleOfReposeTest : public Mercury3D
{
private:

   void setupInitialConditions() override
   {
      stage = 0;
		nBigInserted = 0;
		nSmallInserted = 0;
		t0 = getTime();

		waitForParticleSettling = false;

		std::cout << "SETTING SPECIES..." << std::endl;
		setSpecies();
		std::cout << "MAKING GEOMETRY..." << std::endl;
		makeGeometry();
		std::cout << "COMPUTING NUMBER OF PARTICLES..." << std::endl << std::endl;
		computeNumberOfParticles();
		std::cout << "N. OF BIG PARTICLES NEEDED: " << nBig << std::endl;
		std::cout << "N. OF SMALL PARTICLES NEEDED: " << nSmall << std::endl;
		std::cout << "TOTAL PARTICLES NEEDED: " << nSmall + nBig << std::endl << std::endl;

		std::cout  << "PARTICLE INSERTION" << std::endl << std::endl;

		stage++;
   }

   void actionsAfterTimeStep() override
   {
      if (stage == 1) insertParticles();

		if (stage == 2)
		{
			std::cout << "OPENING OUTLETS" << std::endl << std::endl;
			openOutlets();
			t0 = getTime();

			stage++;
		}

		if (stage == 3)
		{
			deleteParticles();
			if (getKineticEnergy()/getElasticEnergy() < maximumEnergyRatio && getTime() > t0 + 10.0*getTimeStep())
			{
				std::cout << std::endl << "PARTICLE DISCHARGE TERMINATED" << std::endl << std::endl;
				stage++;
			}
		}

		if (stage == 4)
		{
         std::cout << "Computing angle of repose..." << std::endl;
			computeLocalBedHeights();
			computeAngleOfRepose();

         std::cout << "Computing variation and updating AoR and variation arrays..." << std::endl;
         computeAndUpdateAoRAndVariationArrays();

         std::cout << "Updating initial mu array..." << std::endl;
         updateInitialMuArray();

			std::cout << "Writing data to file..." << std::endl;
			writeDataToOutputFile();

			std::cout << "Quitting..." << std::endl;
			setTimeMax(getTime() + getTimeStep());

			stage++;
		}
   }

   void actionAfterSolve()
   {
      delete [] localAngleOfRepose;
   }

public:

   // FUNCTIONS CALLED IN MAIN   --------------------------------------
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

   void setCasingProperties(double heightRatio, double lengthRatio, double widthRatio, double outletRatio, double bedPackingFraction, double density)
   {
      particleBedHeight = radiusBig*heightRatio;
      casingHeight = 2.0*particleBedHeight;
      plateLength = radiusBig*lengthRatio;
      plateWidth = radiusBig*widthRatio;
      outletLength = radiusBig*outletRatio;

      particleBedPackingFraction = bedPackingFraction;
      densityWall = density;

      setXMin(-0.5*plateLength - outletLength);
      setYMin(-0.5*plateWidth);
      setZMin(-outletLength);

      setXMax(0.5*plateLength + outletLength);
      setYMax(0.5*plateWidth);
      setZMax(casingHeight);
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

   void setTargetAngleOfRepose(double target)
   {
      targetAoR = target;
   }

   void setMaximumEnergyRatio(double maxRatio)
   {
      maximumEnergyRatio = maxRatio;
   }

   void setCalibrationRoutineParameters(int nRun, int nMax, double *muSArray, double *devArray, double *aorArray)
   {
      nMaximumRoutineRuns = nMax;
      nCurrentRoutineRun = nRun;

      pointerToInitialMuArray = muSArray;
      pointerToDeviationArray = devArray;
      pointerToAoRArray = aorArray;
   }

   void setGridRelativeSize(double xRelSize, double yRelSize)
   {
      xGridRelativeSize = xRelSize;
      yGridRelativeSize = yRelSize;

      xGridSize = xGridRelativeSize*radiusBig;
      yGridSize = yGridRelativeSize*radiusBig;

      nXCells = (int)(plateLength/xGridSize);
      nYCells = (int)(plateWidth/yGridSize);
   }


   // FUNCTIONS CALLED IN THE CLASS   ---------------------------------
   // sets up the species with their properties
   void setSpecies()
   {
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

   // computes the number of particles
   void computeNumberOfParticles()
   {
      double particleBedVolume = 1.1*(plateLength + 2.0*outletLength)*plateWidth*particleBedHeight;
      nSmall = (particleBedPackingFraction*particleBedVolume/volumeSmall)*(smallToBigMassRatio*densityBig/(densitySmall + smallToBigMassRatio*densityBig));
      nBig = (particleBedPackingFraction*particleBedVolume - nSmall*volumeSmall)/volumeBig;
   }

   // generates teh geometry of the system
   void makeGeometry()
   {
      wallHandler.clear();

      shutter.setSpecies(speciesWall);
      shutter.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
      wallHandler.copyAndAddObject(shutter);

      floor.setSpecies(speciesWall);
      floor.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
      wallHandler.copyAndAddObject(floor);

      ceiling.setSpecies(speciesWall);
      ceiling.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
      wallHandler.copyAndAddObject(ceiling);

      xMinWall.setSpecies(speciesWall);
      xMinWall.set(Vec3D(-1.0,0.0,0.0),Vec3D(getXMin(),0.0,0.0));
      wallHandler.copyAndAddObject(xMinWall);

      xMaxWall.setSpecies(speciesWall);
      xMaxWall.set(Vec3D(1.0,0.0,0.0),Vec3D(getXMax(),0.0,0.0));
      wallHandler.copyAndAddObject(xMaxWall);

      yBoundary.set(Vec3D(0.0,1.0,0.0), getYMin(), getYMax());
      boundaryHandler.copyAndAddObject(yBoundary);

      plate.setSpecies(speciesWall);
      plate.addObject(Vec3D(1.0,0.0,0.0),Vec3D(-0.5*plateLength,0.0,0.0));
      plate.addObject(Vec3D(-1.0,0.0,0.0),Vec3D(0.5*plateLength,0.0,0.0));
      plate.addObject(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
      plate.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,-0.5*outletLength));
      wallHandler.copyAndAddObject(plate);
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
         particlePosition.X = (getXMax() - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Y = (getYMax() - 1.01*p0.getRadius())*random.getRandomNumber(-1.0,1.0);
         particlePosition.Z = random.getRandomNumber(particleBedHeight + 1.01*(p0.getRadius()),getZMax() - 1.01*(p0.getRadius()));

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
      if (!waitForParticleSettling && nBigInserted >= nBig && nSmallInserted >= nSmall && getKineticEnergy()/getElasticEnergy() < maximumEnergyRatio)
      {
         std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
         stage++;
      }
   }

   void openOutlets()
   {
      wallHandler.removeObject(shutter.getId());
   }

   void deleteParticles()
   {
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         if (particleHandler.getObject(i) -> getPosition().Z < -0.5*outletLength) particleHandler.removeObject(i);
         if (fabs(particleHandler.getObject(i) -> getPosition().X) > getXMax() - 1.1*particleHandler.getObject(i) -> getRadius()) particleHandler.removeObject(i);
      }
   }

   void computeLocalBedHeights()
   {
      offsetX = 0.5*(plateLength - nXCells*xGridSize);
      offsetY = 0.5*(plateWidth - nYCells*yGridSize);

      double localMaxParticleHeight;
      localBedHeight = new Vec3D[nXCells*nYCells];

      std::cout << "N X CELLS " << nXCells << std::endl;
      std::cout << "N Y CELLS " << nYCells << std::endl;

      for (int i=0; i<nXCells; i++)
      {
         for (int j=0; j<nYCells; j++)
         {
            localBedHeight[i + nXCells*j].Z = 0.0;

            for (int n=particleHandler.getNumberOfObjects()-1; n>=0; n--)
            {
               if (fabs(particleHandler.getObject(n) -> getPosition().X - (-0.5*plateLength + xGridSize*(0.5 + i) + offsetX)) <= 0.5*xGridSize &&
               fabs(particleHandler.getObject(n) -> getPosition().Y - (-0.5*plateWidth + yGridSize*(0.5 + j) + offsetY)) <= 0.5*yGridSize &&
               particleHandler.getObject(n) -> getPosition().Z + particleHandler.getObject(n) -> getRadius() > localBedHeight[i + nXCells*j].Z)
               {
                  localBedHeight[i + nXCells*j].Z = particleHandler.getObject(n) -> getPosition().Z + particleHandler.getObject(n) -> getRadius();
                  localBedHeight[i + nXCells*j].X = particleHandler.getObject(n) -> getPosition().X;
                  localBedHeight[i + nXCells*j].Y = particleHandler.getObject(n) -> getPosition().Y;
               }
            }
         }
      }
   }

   // void computeAngleOfRepose()
   // {
   //    localAngleOfRepose = new double[(nXCells-1)*nYCells];
   //
   //    for (int i=0; i<nXCells-1; i++)
   //    {
   //       for (int j=0; j<nYCells; j++)
   //       {
   //          localAngleOfRepose[i + (nXCells - 1)*j] = atan(fabs(localBedHeight[i + nXCells*j].Z - localBedHeight[i + 1 + nXCells*j].Z)/fabs(localBedHeight[i + nXCells*j].X - localBedHeight[i + 1 + nXCells*j].X));
   //       }
   //    }
   //
   //    meanAngleOfRepose = 0.0;
   //    int counter = 0;
   //    for (int i=0; i<nXCells-1; i++)
   //    {
   //       for (int j=0; j<nYCells; j++)
   //       {
   //          if (!isnan(localAngleOfRepose[i + (nXCells - 1)*j]))
   //          {
   //             meanAngleOfRepose += localAngleOfRepose[i + (nXCells - 1)*j];
   //             counter++;
   //          }
   //       }
   //    }
   //
   //    meanAngleOfRepose /= counter;
   //    std::cout << "MEAN ANGLE OF REPOSE IS " << meanAngleOfRepose << std::endl;
   // }

   void computeAngleOfRepose()
   {
      localAngleOfRepose = new double[2*nYCells];
      double dZ, dX;

      // for every Y slice finds P_xMax, P_xMin and P_zMax
      // then computes alpha_1 = ArcTan[(z(PzMax)-z(PxMax))/Abs[x(PzMax)-x(PxMax)]]
      // and alpha_2 = ArcTan[(z(PzMax)-z(PxMin))/Abs[x(PzMax)-x(PxMin)]]

      meanAngleOfRepose = 0.0;
      for (int j=0; j<nYCells; j++)
      {
         // AoR between highest particle and xMax
         dZ = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z - (particleHandler.getHighestPositionComponentParticle(0) -> getPosition()).Z;
         dX = fabs((particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).X - (particleHandler.getHighestPositionComponentParticle(0) -> getPosition()).X);
         localAngleOfRepose[2*j] = atan(dZ/dX);
         meanAngleOfRepose += atan(dZ/dX);

         // AoR between highest particle and xMin
         dZ = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z - (particleHandler.getLowestPositionComponentParticle(0) -> getPosition()).Z;
         dX = fabs((particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).X - (particleHandler.getLowestPositionComponentParticle(0) -> getPosition()).X);
         localAngleOfRepose[2*j + 1] = atan(dZ/dX);
         meanAngleOfRepose += atan(dZ/dX);
      }

      meanAngleOfRepose /= 2*nYCells;
      std::cout << "MEAN ANGLE OF REPOSE IS " << meanAngleOfRepose*180.0/constants::pi << std::endl;
   }

   // updates the external array of the deviations
   void computeAndUpdateAoRAndVariationArrays()
   {
      pointerToAoRArray[nCurrentRoutineRun] = meanAngleOfRepose/constants::pi*180;
      relativeDeviationAoR = fabs(meanAngleOfRepose/constants::pi*180 - targetAoR)/targetAoR;
      pointerToDeviationArray[nCurrentRoutineRun] = relativeDeviationAoR;
   }

   // updates the initial values of the next calibration run
   void updateInitialMuArray()
   {
      if (nCurrentRoutineRun == 1)
      {
         pointerToInitialMuArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialMuArray[0] + pointerToInitialMuArray[1]);
      }

      if (nCurrentRoutineRun >= 2)
      {
         if (pointerToAoRArray[nCurrentRoutineRun] > targetAoR)
         {
            if (pointerToAoRArray[nCurrentRoutineRun - 2] > pointerToAoRArray[nCurrentRoutineRun - 1])
            {
               pointerToInitialMuArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialMuArray[nCurrentRoutineRun] + pointerToInitialMuArray[nCurrentRoutineRun - 1]);
            }
            else
            {
               pointerToInitialMuArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialMuArray[nCurrentRoutineRun] + pointerToInitialMuArray[nCurrentRoutineRun - 2]);
            }
         }
         else
         {
            if (pointerToAoRArray[nCurrentRoutineRun - 2] > pointerToAoRArray[nCurrentRoutineRun - 1])
            {
               pointerToInitialMuArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialMuArray[nCurrentRoutineRun] + pointerToInitialMuArray[nCurrentRoutineRun - 2]);
            }
            else
            {
               pointerToInitialMuArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialMuArray[nCurrentRoutineRun] + pointerToInitialMuArray[nCurrentRoutineRun - 1]);
            }
         }
      }
   }

   void writeDataToOutputFile()
   {
      std::ostringstream aorName;
      std::cout.unsetf(std::ios::floatfield);
      aorName << getName() << ".aordata";

      std::ofstream outputFile;
      outputFile.open(aorName.str(), std::ios::out);
      outputFile << "Target angle of repose: " << targetAoR << std::endl;
      outputFile << "N. X cells: " << nXCells << std::endl;
      outputFile << "N. Y cells: " << nYCells << std::endl;
      outputFile << "X grid size: " << xGridSize << std::endl;
      outputFile << "Y grid size: " << yGridSize << std::endl;
      outputFile << "Plate length: " << plateLength << std::endl;
      outputFile << "Plate width: " << plateWidth << std::endl;
      outputFile << "Particle bed height: " << particleBedHeight << std::endl << std::endl;

      for (int i=0; i<nXCells; i++)
      {
         for (int j=0; j<nYCells; j++)
         {
            outputFile <<
            i << "   " << j << "   " <<
            -0.5*plateLength + xGridSize*(0.5 + i) + offsetX << "   " <<
            -0.5*plateWidth + yGridSize*(0.5 + j) + offsetY << "   " <<
            localBedHeight[i + nXCells*j].X << "   " <<
            localBedHeight[i + nXCells*j].Y << "   " <<
            localBedHeight[i + nXCells*j].Z << "   " <<
            std::endl;
         }
      }

      // outputFile << std::endl << "Local angle of repose (radians and degrees):" << std::endl;;
      //
      // for (int i=0; i<nXCells-1; i++)
      // {
      //    for (int j=0; j<nYCells; j++)
      //    {
      //       outputFile << i << "   " << j << "   " << localAngleOfRepose[i + (nXCells - 1)*j] << "   (" << localAngleOfRepose[i + (nXCells - 1)*j]/constants::pi*180 << ")" << std::endl;
      //    }
      // }

      outputFile << std::endl << "---------- THIS IS THE NEW AOR COMPUTATION ----------" << std::endl;
      outputFile << "MEAN ANGLE OF REPOSE IS: " << meanAngleOfRepose << "   (" << meanAngleOfRepose/constants::pi*180 << ")" << std::endl;
      outputFile << "|<alpha> - alpha_target|/alpha_target = " << relativeDeviationAoR << std::endl;

      outputFile.close();
   }


   //  ----- GLOBAL FUNCTIONS -----
    void printTime() const override
    {
        std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
		", N. of particles = " << particleHandler.getNumberOfObjects() << ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
        std::cout.flush();
    }


   // VARIABLES  ------------------------------------------------------
   // particles
   double radiusBig, radiusSmall;
   double sizeDispersityBig, sizeDispersitySmall;
   double densityBig, densitySmall;
   double volumeBig, volumeSmall;
   double massBig, massSmall;
   double smallToBigMassRatio;
   double particleTotalVolume;

   // powder bed
   double particleBedHeight;
   double particleBedPackingFraction;

   // geometry
   double casingHeight, plateLength, plateWidth, outletLength;
   double densityWall;
   InfiniteWall floor, shutter, ceiling;
   InfiniteWall xMinWall, xMaxWall;
   InfiniteWall yMinWall, yMaxWall;
   IntersectionOfWalls plate;
   PeriodicBoundary yBoundary;

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

   // loading variables
   double nBig, nSmall;
   int nBigInserted, nSmallInserted;
   double maximumEnergyRatio;

   // calibration
   double targetAoR;
   double relativeDeviationAoR;
   double *pointerToInitialMuArray;
   double *pointerToAoRArray;
   double *pointerToDeviationArray;

   // angle of repose
   double xGridRelativeSize, yGridRelativeSize;
   double xGridSize, yGridSize;
   double offsetX, offsetY;
   int nXCells, nYCells;
   Vec3D *localBedHeight;
   double *localAngleOfRepose;
   double meanAngleOfRepose;

   // global
   int stage;
   double t0;
   bool waitForParticleSettling;
   int nMaximumRoutineRuns;
   int nCurrentRoutineRun;
};

// returns the coefficient to find the cohesiveness based on the target angle of repose
double returnInitialCohesionStiffnessMultiplier(double targetAoR)
{
   if (targetAoR < 25.0) return 1000.0;
   else if (targetAoR < 30.0) return 500.0;
   else if (targetAoR < 35.0) return 200.0;
   else if (targetAoR < 40.0) return 100.0;
   else if (targetAoR < 45.0) return 50.0;
   else if (targetAoR < 50.0) return 20.0;
   else return 10.0;
}

// initializes the calibration result arrays
void initializeCalibrationArrays(int nRunsMax, double *initMuS, double *aorArray, double *devArray, double muSBB)
{
   for (int i = 0; i < nRunsMax; i++)
   {
      initMuS[i] = 0.0;
      aorArray[i] = 0.0;
      devArray[i] = 1.0;
   }

   initMuS[0] = muSBB;
   initMuS[1] = 0.1*muSBB;
}

int main(int argc, char *argv[])
{
   // DATA FILENAME  SPECIFICATION
   std::string material = "CAPSULAC";

   // TARGET AoR
   double targetAoR = 28.81;

   // TIME STEP
   double timeStep = 5.0e-6;

   // GEOMETRY PARAMETERS
   double casingRadius = 0.0125;
   double particleBedHeightToParticleRadiusSizeRatio = 80.0; //80
	double plateLengthToParticleRadiusSizeRatio = 100.0; //100
	double plateWidthToParticleRadiusSizeRatio = 20.0; //20
	double outletLengthToParticleradiusSizeRatio = 20.0; //20
	double particleBedPackingFraction = 0.60;
   double densityWalls = 5000.0;

   // PARTICLE STATIC PARAMETERS
   double bigToCasingSizeRatio = 40.0;
   double bigToSmallSizeRatio = 1.0;
   double smallToBigMassRatio = 0.0;
   double sizeDispersionBigParticles = 0.1;
   double sizeDispersionSmallParticles = 0.1;
   double densityBigParticles = 1500.0;
   double densitySmallParticles = 1500.0;

   // INTERACTION PARAMETERS
   double k2MaxBig = 3000.0;
   double k1Big = k2MaxBig/10.0;
   double kCBig = 12.0;
   double phiBig = 0.05;

   double k2MaxSmall = 3000.0;
   double k1Small = k2MaxSmall/10.0;
   double kCSmall = k2MaxSmall/returnInitialCohesionStiffnessMultiplier(targetAoR);
   double phiSmall = 0.05;

   double k1Wall = k2MaxBig;

   double eBigWall = 0.7;
   double eSmallWall = 0.7;
   double eBigBig = 0.5;
   double eSmallSmall = 0.5;
   double eBigSmall = 0.5;

   // FRICTION PARAMETERS
   double muSBigBig = 0.50;
   double muRBigBig = 0.05;
   double muSBigWall = muSBigBig;
   double muRBigWall = muRBigBig;

   // GRID PARAMETERS
	double xGridRelativeSize = 2.0;
	double yGridRelativeSize = 2.0;

   // CALIBRATION PARAMETERS
   double maximumAoRRelativeDeviation = 0.05;
   double maximumEnergyRatio = 1.0e-4;

   // CALIBRATION ARRAYS
   int runNumber = 0;
   int nMaximumRoutineCycles = 10;
   double initialMuArray[nMaximumRoutineCycles];
   double angleOfReposeArray[nMaximumRoutineCycles];
   double deviationArray[nMaximumRoutineCycles];

   // CALIBRATION DATA POINTS AND ARRAYS INITIALIZATION
   initializeCalibrationArrays(nMaximumRoutineCycles, initialMuArray, angleOfReposeArray, deviationArray, muSBigBig);

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

   // THE CALIBRATION LOOP
   do
   {
      CalibrationRoutine_AngleOfReposeTest crAoR;

      // NAME SETTING
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      name << "CalibrationTest_AoRTEST_" << material << "_Nrun_" << runNumber << "_kc_" << std::fixed << std::setprecision(0) << kCBig <<
      "_muS_" << std::fixed << std::setprecision(2) << initialMuArray[runNumber] << "_muR_" << muRBigBig << "___corrected";
      crAoR.setName(name.str());

      std::cout << "The name of the associated output of this simulation will be:" << std::endl << name.str() << std::endl << std::endl;

      // INITIALIZATION
      crAoR.setTimeStep(timeStep);
      crAoR.setTimeMax(100.0);
      crAoR.setGravity(Vec3D(0.00,0.00,-9.81));
      crAoR.setSystemDimensions(3);

      crAoR.setSaveCount(0.01/crAoR.getTimeStep());

      crAoR.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);
      crAoR.setCasingProperties(particleBedHeightToParticleRadiusSizeRatio, plateLengthToParticleRadiusSizeRatio, plateWidthToParticleRadiusSizeRatio, outletLengthToParticleradiusSizeRatio, particleBedPackingFraction, densityWalls);

      crAoR.setParticleWallSlidingFrictionCoefficients(initialMuArray[runNumber], initialMuArray[runNumber]);
      crAoR.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
      crAoR.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);

      crAoR.setParticleParticleSlidingFrictionCoefficients(initialMuArray[runNumber], initialMuArray[runNumber], initialMuArray[runNumber]);
      crAoR.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
      crAoR.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);

      crAoR.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
      crAoR.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
      crAoR.setBigParticlePlasticProperties(k1Big, k2MaxBig, kCBig, phiBig);
      crAoR.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);

      crAoR.setGridRelativeSize(xGridRelativeSize, yGridRelativeSize);

      crAoR.setTargetAngleOfRepose(targetAoR);
      crAoR.setMaximumEnergyRatio(maximumEnergyRatio);
      crAoR.setCalibrationRoutineParameters(runNumber, nMaximumRoutineCycles, initialMuArray, deviationArray, angleOfReposeArray);

      crAoR.solve();

      // BOUNDARY TESTS
      if (runNumber == 0 && angleOfReposeArray[0] < targetAoR)
      {
         std::cout << "alpha_target > alpha(mu_S = 0.50)" << std::endl << "kC needs to be increased. Quitting..." << std::endl;
         exit(1);
      }
      if (runNumber == 1 && angleOfReposeArray[1] > targetAoR)
      {
         std::cout << "alpha_target < alpha(mu_S = 0.05)" << std::endl << "kC needs to be decreased. Quitting..." << std::endl;
         exit(1);
      }

      runNumber++;
   }
   while (runNumber < nMaximumRoutineCycles && deviationArray[runNumber - 1] > maximumAoRRelativeDeviation);

   return 0;
}
