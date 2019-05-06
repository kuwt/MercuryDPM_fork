#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

// VERSION: PHI DYNAMICALLY TUNED

/*
NOTES:
09.07:
- the phi change in refreshSpeciesBig() was commented out, now it is in
- corrected the updateInitialKArray() that would save the wrong phi value
*/

class CompressionTest_calibrationRoutine : public Mercury3D
{
private:
   void setupInitialConditions() override
   {
      std::cout << "THIS MESSAGE SHOULD NOT BE DISPLAYED." << std::endl;
      std::cout << "SINCE YOU ARE READING IT SOMETHING WENT WRONG." << std::endl;
      std::cout << "EXITING..." << std::endl;

      exit(1);
   }

   void actionsOnRestart() override
   {
      stage = 0;
      setTime(0.0);
      setTimeMax(10.0);

      t0 = getTime();
      calibrationStep = 0;
      meanK1Big = 0.0;
      phiMin = (1.0 - k1Big/k2MaxBig)*deltaHatMax;
      isHeightTargetMet = false;
      isPressureTargetMet = false;

      // species set and reset
      setSpecies();
      resetSpecies();

      // setup of all the rtest
      makePiston();
      updateTargetDataPoint();
      resetPistonVelocity();
      makeCdatFile();
   }

   void actionsAfterTimeStep() override
   {
      makeDataAnalysis();
      if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) writeToCdatFile();

      if (stage == 0) runPressureCalibration();
      if (stage == 1) unloadPiston();

      // TEST DYNAMIC OUTPUT
      if (fmod(getTime(), 2.0*cdatOutputTimeInterval) < getTimeStep())
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight <<
         ", v = " << pistonVelocity << std::endl <<
         "piston pointer h = " << pistonPointer -> getPosition().Z << ", piston pointer v = " << pistonPointer -> getVelocity().Z << ", hMaxParticle = " << (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z << std::endl <<
         "P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << " , P_target = " << targetPressure <<
         ", P/P_target = " << (pistonPressure + 0.01)/(targetPressure + 0.01) << ", |(P - P_target)/P_target| = " << fabs((pistonPressure - targetPressure)/(targetPressure + 0.01)) << std::endl <<
         "cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
         ", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap <<
         ", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
         "k1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() <<
         ", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kC = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() <<
         ", phi_min = " << phiMin << ", phi_min/phi = " << phiMin/(speciesBig -> getPenetrationDepthMax()) <<std::endl;

         std::cout << std::endl;
      }
   }

   void actionsAfterSolve() override
   {
      delete [] k1CalibratedArray;
      delete [] heightSteps;
      delete [] pressureSteps;

      cdatFile.close();
   }

public:
   // FUNCTIONS CALLED IN MAIN  ---------------------------------------
   void setCdatOutputTimeInterval(double dt)
   {
      cdatOutputTimeInterval = dt;
   }

   void setCasingProperties(double radius, double height, double density)
   {
      casingRadius = radius;
      particleBedHeight = height;
      densityWall = density;
      //
      // setXMin(0.0);
      // setYMin(0.0);
      // setZMin(0.0);
      //
      // setXMax(casingRadius);
      // setYMax(casingRadius);
      // setZMax(casingHeight);
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

   void setStiffnessCalibrationDataPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double pauseDuration)
   {
      nCalibrationDataPoints = n;
      pauseAfterDataPointMet = pauseDuration;
      k1CalibratedArray = new double[n];
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

   void setCalibrationRoutineParameters(int nRun, int nMax, double *k1Array, double *phiArray, double *devArray, double *dispArray)
   {
      nMaximumRoutineRuns = nMax;
      nCurrentRoutineRun = nRun;

      pointerToInitialK1Array = k1Array;
      pointerToInitialPhiArray = phiArray;

      pointerToDeviationArray = devArray;
      pointerToDisplacementArray = dispArray;
   }

   void setNTimeStepsBetweenStiffnesAdjustments(double dt)
   {
      nTimeStepsBetweenStiffnessAdjustments = dt;
   }

   void setPistonVelocityScalingFactor(double factor)
   {
      pistonVelocityScalingFactor = factor;
   }

   void setCConstantAndDeltaHatMax(double C, double Dmax)
   {
      cConstant = C;
      deltaHatMax = Dmax;
   }


   // FUNCTIONS CALLED IN THE CLASS  ----------------------------------
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

   // sets the species of the loaded object to the new ones declared in setSpecies()
   void resetSpecies()
   {
      // wall interactables
      for (int i = 0; i < wallHandler.getNumberOfObjects(); i++) wallHandler.getObject(i) -> setSpecies(speciesWall);

      // particle interactables
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++) particleHandler.getObject(i) -> setSpecies(speciesBig);
   }

   // creates the flat compression piston
   void makePiston()
   {
      std::cout << "CREATING PISTON" << std::endl << std::endl;

      pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius());
      pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

      // compressionPiston.setSpecies(speciesWall);
      compressionPiston.setSpecies(speciesBig);
      compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(compressionPiston);
   }

   // analyzes simulations data of interest
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

            pistonInteractionCounter++;
         }

         // base interactions
         if ((*i) -> getI() -> getIndex() == 0)
         {
            basePressure += ((*i) -> getForce()).Z;

            meanBaseRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxBaseRelativeOverlap) maxBaseRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            baseInteractionCounter++;
         }

         // cylinder interactions
         if ((*i) -> getI() -> getIndex() == 2)
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

   // sets the next data point to calibrate, and resets the height/pressure flags
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

   // resets the piston velocity to the starting value
   void resetPistonVelocity()
   {
      if (previousTargetHeight < targetHeight) pistonVelocity = pistonMaximumVelocity;
      else pistonVelocity = -pistonMaximumVelocity;

      std::cout << "RESETTING PISTON VELOCITY. NEW VELOCITY: " << pistonVelocity << std::endl << std::endl;
   }

   // displaces the piston position according to its velocity
   void movePiston()
   {
      pistonHeight += pistonVelocity*getTimeStep();
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));
   }

   // stops the piston by setting to 0 its velocity
   void stopPiston()
   {
      pistonVelocity = 0.0;
      pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));
   }

   // resets the species for the big particles, according to the new loading stiffness
   void refreshSpeciesBig()
   {
      // BIG-BIG
      speciesBig -> setStiffnessAndRestitutionCoefficient(k1Big, bigBigRestitutionCoeff, massBig);
      speciesBig -> setUnloadingStiffnessMax(k2MaxBig);
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
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setPenetrationDepthMax(0.5*(phiBig + phiSmall));

      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
   }

   // adjusts the loading stiffness of big particles according to the piston pressure
   void tuneK1Big()
   {
      if (pistonPressure < targetPressure) k1Big *= 1.001;
      else k1Big *= 0.999;

      phiBig = pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax/(cConstant - 1.0);
      phiMin = (1.0 - k1Big/k2MaxBig)*deltaHatMax;

      refreshSpeciesBig();

      if (k1Big > k2MaxBig)
      {
         std::cout << "\n\nFATAL ERROR: K1_BIG > K2MAX_BIG. QUITTING.";
         exit(1);
      }

      if (k1Big < 0.0)
      {
         std::cout << "\n\nFATAL ERROR: K1_BIG < 0.0. QUITTING.";
         exit(1);
      }
   }

   // the pressure calibration loop
   void runPressureCalibration()
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

            stopPiston();
            isHeightTargetMet = true;
         }
      }

      // when the height is fine deal with the pressure
      if (isHeightTargetMet && !isPressureTargetMet && getTime() > t0 + pauseAfterDataPointMet)
      {
         if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
         {
            // the stiffness is changed every n time steps
            if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneK1Big();
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
         calibrationStep++;

         if (calibrationStep < nCalibrationDataPoints)
         {
            updateTargetDataPoint();
            resetPistonVelocity();
         }
         else
         {
            targetHeight = particleBedHeight;
            resetPistonVelocity();

            std::cout << "Compression cycles finished. Unloading piston...\n";
            //meanK1Big = computeMean(k1CalibratedArray, nCalibrationDataPoints);
            //pointerToDeviationArray[nCurrentRoutineRun] = computeDeviation(k1CalibratedArray, nCalibrationDataPoints);
            //pointerToDisplacementArray[nCurrentRoutineRun] = computeDisplacement(pointerToInitialK1Array[nCurrentRoutineRun], k1CalibratedArray, nCalibrationDataPoints);

            t0 = getTime();
            stage++;
         }
      }
   }

   // computes averages, deviations and displacements of the calibrated values of K1
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

   // updates the external array of the deviations
   void updateVariationArrays()
   {
      pointerToDeviationArray[nCurrentRoutineRun] = meanRelativeDeviationK1Big;
      pointerToDisplacementArray[nCurrentRoutineRun] = meanRelativeDisplacementeK1Big;
   }

   // unloads the piston by moving it to the original powder bed height
   void unloadPiston()
   {
      if (fabs((pistonHeight - targetHeight)/targetHeight) > maxHeightRelDeviation) movePiston();
      else
      {
         computeVariations();
         updateVariationArrays();

         std::cout << "Piston unloaded. Quitting...\n";
         //meanK1Big = computeMean(k1CalibratedArray, nCalibrationDataPoints);
         //pointerToDeviationArray[nCurrentRoutineRun] = computeDeviation(k1CalibratedArray, nCalibrationDataPoints);
         //pointerToDisplacementArray[nCurrentRoutineRun] = computeDisplacement(pointerToInitialK1Array[nCurrentRoutineRun], k1CalibratedArray, nCalibrationDataPoints);

         makeKdataFile();
         if (nCurrentRoutineRun < nMaximumRoutineRuns) updateInitialKArray();

         setTimeMax(getTime() + getTimeStep());
      }
   }

   // updates the initial values of the next calibration run
   void updateInitialKArray()
   {
      pointerToInitialK1Array[nCurrentRoutineRun + 1] = meanK1Big;
      pointerToInitialPhiArray[nCurrentRoutineRun + 1] = pow(1.0 - pointerToInitialK1Array[nCurrentRoutineRun + 1]/k2MaxBig, 2.0)/(pointerToInitialK1Array[nCurrentRoutineRun + 1]/k2MaxBig)*deltaHatMax/(cConstant - 1.0);
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
      kdataFile << "kC = " << kCBig << std::endl;
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
      phiBig << "   " <<
      std::endl;
   }


   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
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


   // VARIABLES  ------------------------------------------------------
   // particles
   double radiusBig, radiusSmall;
   double sizeDispersityBig, sizeDispersitySmall;
   double densityBig, densitySmall;
   double volumeBig, volumeSmall;
   double massBig, massSmall;
   double smallToBigMassRatio;

   // powder bed
   double particleBedHeight;

   // geometry
   double casingRadius;
   double densityWall;

   // piston
   double pistonHeight;
   double pistonVelocity;
   double pistonVelocityScalingFactor;
   double pistonMaximumVelocity;
   double pistonForce;
   double pistonPressure;
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
   double *heightSteps;
   double *pressureSteps;
   bool isHeightTargetMet;
   bool isPressureTargetMet;
   int calibrationStep;
   double targetHeight;
   double targetPressure;
   double previousTargetHeight;
   double maxHeightRelDeviation;
   double maxPressureRelDeviation;
   double nTimeStepsBetweenStiffnessAdjustments;
   double *pointerToInitialK1Array;
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
   double t0;
   int nMaximumRoutineRuns;
   int nCurrentRoutineRun;
};


// computes the delta intersection
double deltaIntersection(double d1, double d2, double p1, double p2)
{
   return d1 - p1*(d2 - d1)/(p2 - p1);
}

// reads the experiment data from file, initializes the target values for the calibration and computes the deltaHatMax and the cConstant
void readExperimentDataAndComputeDeltaHatMaxAndCConstant(int nTot, std::string inputFileName, double *hArray, double *pArray, double h0, double *dMax, double *cC)
{
   int indexOfHighestPressurePoint;
   double highestPressure = 0.0;

   std::ifstream inputFile(inputFileName);
   if (inputFile.is_open())
   {
      for (int i = 0; i < nTot; i++)
      {
         inputFile >> hArray[i];
         inputFile >> pArray[i];
         hArray[i] = h0 - hArray[i];
         if (pArray[i] > highestPressure)
         {
            highestPressure = pArray[i];
            indexOfHighestPressurePoint = i;
         }
      }

      inputFile.close();
   }
   else
   {
      std::cout << "DATA FILE NOT FOUND." << std::endl << "QUITTING..." << std::endl;
      exit(1);
   }

   // computes the deltaHatMax
   dMax[0] = (h0 - hArray[indexOfHighestPressurePoint])/h0;

   // computes the cConstant
   cC[0] = (h0 - hArray[indexOfHighestPressurePoint])/
   (h0 - hArray[indexOfHighestPressurePoint] - deltaIntersection(h0 - hArray[nTot - 1], h0 - hArray[indexOfHighestPressurePoint], pArray[nTot - 1], pArray[indexOfHighestPressurePoint]));

   // cC[0] = (h0 - hArray[indexOfHighestPressurePoint] - deltaIntersection(h0 - hArray[0], h0 - hArray[indexOfHighestPressurePoint], pArray[0], pArray[indexOfHighestPressurePoint]))/
   // (h0 - hArray[indexOfHighestPressurePoint] - deltaIntersection(h0 - hArray[nTot - 1], h0 - hArray[indexOfHighestPressurePoint], pArray[nTot - 1], pArray[indexOfHighestPressurePoint]));

   std::cout << "dMaxHat " << dMax[0] << std::endl;
   std::cout << "cC " << cC[0] << std::endl;
}

// initializes the calibration result arrays
void initializeCalibrationArrays(int nRunsMax, double *initK1, double *initPhi, double *devArray, double *dispArray, double k1BigInit, double phiBigInit)
{
   for (int i = 0; i < nRunsMax; i++)
   {
      initK1[i] = 0.0;
      initPhi[i] = 0.0;
      devArray[i] = 1.0;
      dispArray[i] = 1.0;
   }

   initK1[0] = k1BigInit;
   initPhi[0] = phiBigInit;
}

int main(int argc, char *argv[])
{
   // immediately aborts if not enough arguments were given in command line
   if (argc != 3)
   {
      std::cout << "Not enought input arguments given." << std::endl;
      std::cout << "QUITTING..." << std::endl;

      exit(1);
   }

   // DATA FILENAME  SPECIFICATION
   std::string material = "MANNITOL";
   std::string dataPointsFileNameCOMPRESSION = "Calibration_dataPoints/Calibration_dataPoints_COMPRESSION_MANNITOL.dat";

   // DATA POINT SELECTION FOR ANALYSIS
   int nDataPointsToAnalyze = 15; // ALWAYS SKIP THE FINAL UNLOADING DATA POINT
   int nTotalDataPoints = 16;
   double pauseDuration = 0.01;
   double heightArray[nTotalDataPoints];
   double pressureArray[nTotalDataPoints];

   // TIME STEP
   double timeStep = 1.0e-6;

   // GEOMETRY PARAMETERS
   double particleBedHeight = 0.019;
   double casingRadius = 0.0125;
   double densityWalls = 3000.0;

   // COMPRESSION PARAMETERS
   double pistonVelocityScalingFactor = 1.0;

   // PARTICLE STATIC PARAMETERS
   double bigToCasingSizeRatio = 40.0;
   double bigToSmallSizeRatio = 1.0;
   double smallToBigMassRatio = 0.0;
   double sizeDispersionBigParticles = 0.1;
   double sizeDispersionSmallParticles = 0.1;
   double densityBigParticles = 1500.0;
   double densitySmallParticles = 1500.0;

   // INITIALIZATION OF DATA POINTS ARRAY AND COMPUTATION OF dHatMax AND cConstant
   double deltaHatMax[1];
   double cConstant[1];
   readExperimentDataAndComputeDeltaHatMaxAndCConstant(nTotalDataPoints, dataPointsFileNameCOMPRESSION, heightArray, pressureArray, particleBedHeight, deltaHatMax, cConstant);

   // INTERACTION PARAMETERS
   double k2MaxBig = 5000.0;
   double k1Big = k2MaxBig/10.0;
   double kCBig = 0.0;
   double phiBig = pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax[0]/(cConstant[0] - 1.0);

   double k2MaxSmall = k2MaxBig;
   double k1Small = k1Big;
   double kCSmall = kCBig;
   double phiSmall = phiBig;

   double k1Wall = k2MaxBig;

   double eBigWall = 0.7;
   double eSmallWall = 0.7;
   double eBigBig = 0.5;
   double eSmallSmall = 0.5;
   double eBigSmall = 0.5;

   // FRICTION PARAMETERS
   double muSBigBig = 0.20;
   double muRBigBig = 0.10;
   double muSBigWall = muSBigBig;
   double muRBigWall = muRBigBig;

   // CALIBRATION PARAMETERS
   double maximumStiffnessRelativeDeviation = 0.05;
   double maximumHeightRelativeDeviation = 1.0e-4;
   double maximumPressureRelativeDeviation = 0.01;
   double timeStepsBetweenStiffnessTunings = 5.0;

   // CALIBRATION ARRAYS
   int runNumber = 0;
   int nMaximumRoutineCycles = 10;
   double initialK1Array[nMaximumRoutineCycles];
   double initialPhiArray[nMaximumRoutineCycles];
   double deviationArray[nMaximumRoutineCycles];
   double displacementArray[nMaximumRoutineCycles];

   // CALIBRATION ARRAYS INITIALIZATION
   initializeCalibrationArrays(nMaximumRoutineCycles, initialK1Array, initialPhiArray, deviationArray, displacementArray, k1Big, phiBig);

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

   std::string phase;
   if (nDataPointsToAnalyze == 6) phase = "LOADING";
   else if (nDataPointsToAnalyze == 10) phase = "RELOADING";
   else if (nDataPointsToAnalyze == nTotalDataPoints - 1) phase = "FULL";
   else phase = "";

   // THE CALIBRATION LOOP
   do
   {
      CompressionTest_calibrationRoutine crCompression;
      crCompression.readArguments(argc, argv);

      // NAME SETTING
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      name << "CompressionTest_" << material << "_" << phase << "_Nrun_" << runNumber <<
      "_k1_" << std::fixed << std::setprecision(0) << initialK1Array[runNumber] << "_k2max_" << k2MaxBig << "_kc_" << kCBig << "_phi_" << std::fixed << std::setprecision(4) << initialPhiArray[runNumber] <<
      "_muS_" << std::fixed << std::setprecision(2) << muSBigBig << "_muR_" << muRBigBig <<
      "_cConstant_" << std::fixed << std::setprecision(4) << cConstant[0] << "_DhatMax_" << deltaHatMax[0];
      crCompression.setName(name.str());

      std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl << std::endl;


      // INITIALIZATION
      crCompression.setTimeStep(timeStep);
      crCompression.setTimeMax(10.0);
      crCompression.setGravity(Vec3D(0.00,0.00,-9.81));
      crCompression.setSystemDimensions(3);
      //crCompression.setVerbose(true);
      crCompression.setCdatOutputTimeInterval(0.001);
      crCompression.setSaveCount(0.01/crCompression.getTimeStep());

      crCompression.setCasingProperties(casingRadius, particleBedHeight, densityWalls);
      crCompression.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);

      crCompression.setParticleWallSlidingFrictionCoefficients(muSBigWall, muSBigWall);
      crCompression.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
      crCompression.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);

      crCompression.setParticleParticleSlidingFrictionCoefficients(muSBigBig, muSBigBig, muSBigBig);
      crCompression.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
      crCompression.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);

      crCompression.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
      crCompression.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
      crCompression.setBigParticlePlasticProperties(initialK1Array[runNumber], k2MaxBig, kCBig, initialPhiArray[runNumber]);
      crCompression.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);

      crCompression.setStiffnessCalibrationDataPoints(nDataPointsToAnalyze, heightArray, pressureArray, pauseDuration);
      crCompression.setCConstantAndDeltaHatMax(cConstant[0], deltaHatMax[0]);
      crCompression.setCalibrationMaximumRelativeDeviations(maximumHeightRelativeDeviation, maximumPressureRelativeDeviation);
      crCompression.setCalibrationRoutineParameters(runNumber, nMaximumRoutineCycles, initialK1Array, initialPhiArray, deviationArray, displacementArray);
      crCompression.setNTimeStepsBetweenStiffnesAdjustments(timeStepsBetweenStiffnessTunings);

      crCompression.setPistonVelocityScalingFactor(pistonVelocityScalingFactor);

      crCompression.solve();
      runNumber++;
   }
   while (runNumber < nMaximumRoutineCycles && !(displacementArray[runNumber - 1] < maximumStiffnessRelativeDeviation && deviationArray[runNumber - 1] < maximumStiffnessRelativeDeviation));

   return 0;
}
