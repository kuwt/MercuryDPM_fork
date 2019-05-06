#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Walls/IntersectionOfWalls.h>
#include "RotatingIntersectionOfWalls.h"
#include <math.h>
#include <fstream>

// VERSION: PRE-SHEAR + PRESSURE CONSTRAINED

/*
NOTES:
22.06:
- removed every reference to height and target height, now the program runs only on target pressure and torque
- there is still teh height array data import, need to get rid of that

26.06:
- add a constant pressure condition in the torque calibration routine and compare with the unconstraint case

09.07:
- added a check: if, while calibrating, kC drops under 1, then switch to next calibration point
- added a check: if, while calibrating, kC gets over k1, then switch to next calibration point

26.07:
- added the flat/non_flat piston options
*/

class ShearTest_calibrationRoutine : public Mercury3D
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
      meanKCBig = 0.0;
      nOfFreeParticles = particleHandler.getNumberOfObjects();
      isPressureTargetMet = false;
      isTorqueTargetMet = false;

      // species set and reset
      setSpecies();
      resetSpecies();

      // setup of all the rest
      // makePiston();
      makeComplexPiston();
      updateTargetDataPoint();
      resetPistonVelocity();
      makeCdatFile();
   }

   void actionsAfterTimeStep() override
   {
      makeDataAnalysis();
      if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) writeToCdatFile();

      // if (stage == 0) runPreShearing();
      if (stage == 0) runTorqueCalibration();

      // TEST DYNAMIC OUTPUT
      if (fmod(getTime(), 5.0*cdatOutputTimeInterval) < getTimeStep())
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight <<
         ", hEff = " << effectivePistonHeight << ", v = " << pistonVelocity << ", omega = " << pistonAngularVelocity << std::endl <<
         "P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << std::endl <<
         "torque = " << pistonTorque << " , torque_target = " << targetTorque << std::endl <<
         "cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
         ", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap <<
         ", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
         "k1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() <<
         ", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kC = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl;
         std::cout << std::endl;
      }
   }

   void actionsAfterSolve() override
   {
      delete [] kCCalibratedArray;
      delete [] torqueSteps;
      cdatFile.close();
   }

public:
   // FUNCTIONS CALLED IN MAIN   --------------------------------------
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

   void setTorqueCalibrationPoints(int n, double *pressureArrayPointer, double *torqueArrayPointer, double pauseDuration, double pauseDuration2)
   {
      nCalibrationDataPoints = n;
      pauseAfterDataPointMet = pauseDuration;
      pauseBeforeTorqueCalibration = pauseDuration2;
      kCCalibratedArray = new double[n];
      pressureSteps = new double[n];
      torqueSteps = new double[n];

      for (int i = 0; i < n; i++)
      {
         pressureSteps[i] = pressureArrayPointer[i];
         torqueSteps[i] = torqueArrayPointer[i];
      }
   }

   void setCalibrationMaximumRelativeDeviations(double maxP, double maxT)
   {
      maxPressureRelDeviation = maxP;
      maxTorqueRelDeviation = maxT;
   }

   void setCalibrationRoutineParameters(int nRun, int nMax, double *kCArray, double *devArray, double *dispArray)
   {
      nMaximumRoutineRuns = nMax;
      nCurrentRoutineRun = nRun;

      pointerToInitialKCArray = kCArray;

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

   void setPistonRotation(double angularVelocity)
   {
      pistonAngularVelocity = angularVelocity;
   }

   void setPreShearingParameters(double t, double a, double f)
   {
      preShearingDuration = t;
      preShearingOscillationAmplitude = a;
      preShearingOscillationFrequency = f;
   }

   void setPistonType(bool isFlat)
   {
      isPistonFlat = isFlat;
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

   // sets the species of the loaded object to the new ones declared in setSpecies()
   void resetSpecies()
   {
      // wall interactables
      for (int i = 0; i < wallHandler.getNumberOfObjects(); i++) wallHandler.getObject(i) -> setSpecies(speciesWall);

      // particle interactables
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++) particleHandler.getObject(i) -> setSpecies(speciesBig);
   }

   // creates the piston made out of particles
   void makePiston()
   {
      if (isPistonFlat)
      {
         std::cout << std::endl << "CREATING SIMPLE FLAT PISTON" << std::endl << std::endl;

         pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius());
         pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

         // compressionPiston.setSpecies(speciesWall);
         compressionPiston.setSpecies(speciesBig);
         compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
         pistonPointer = wallHandler.copyAndAddObject(compressionPiston);

         std::cout << "PARTICLE BED HEIGHT " << (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z << std::endl;
         std::cout << "PISTON HEIGHT " << pistonHeight << std::endl;
      }
      else
      {
         std::cout << std::endl << "CREATING PARTICLE PISTON FROM SAMPLED PARTICLES" << std::endl << std::endl;

         double samplingHeight = 0.5*(particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z;
         double topParticleDisplacement = 0.0;
         double totalParticlesSphericalCapVolume = 0.0;
         effectivePistonHeight = 0.0;
         pistonLatticeParticleCounter = 0;

         // counts the sampled particles and checks for the lowest edge of the former
         for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
         {
            if (fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight) < particleHandler.getObject(i) -> getRadius())
            {
               if ((particleHandler.getObject(i) -> getRadius() - fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight)) > topParticleDisplacement)
               {
                  topParticleDisplacement = particleHandler.getObject(i) -> getRadius() - fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight);
               }

               totalParticlesSphericalCapVolume += (pow(particleHandler.getObject(i) -> getRadius() - fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight),2.0)*(2.0*(particleHandler.getObject(i) -> getRadius()) + fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight)));
               pistonLatticeParticleCounter++;
            }
         }

         // the piston height is function of the particle bed height and of the max distance between the particles bottom edge and the sampling height
         pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius()) + topParticleDisplacement;
         effectivePistonHeight = pistonHeight - totalParticlesSphericalCapVolume/(3.0*pow(casingRadius,2.0));
         pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

         BaseParticle p0;
         p0.setVelocity(Vec3D(0.0,0.0,pistonMaximumVelocity));
         p0.setSpecies(speciesBig);

         Vec3D sampledParticleOriginalPosition;

         for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
         {
            if (fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight) < particleHandler.getObject(i) -> getRadius())
            {
               sampledParticleOriginalPosition = particleHandler.getObject(i) -> getPosition();

               p0.setRadius(particleHandler.getObject(i) -> getRadius());
               p0.setPosition(Vec3D(sampledParticleOriginalPosition.X, sampledParticleOriginalPosition.Y, pistonHeight + fabs(sampledParticleOriginalPosition.Z - samplingHeight)));
               p0.fixParticle();
               particleHandler.copyAndAddObject(p0);
            }
         }

         std::cout << "PARTICLE BED HEIGHT " << (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z << std::endl;
         std::cout << "SAMPLING HEIGHT " << samplingHeight << std::endl;
         std::cout << "EFFECTIVE PISTON HEIGHT " << effectivePistonHeight << std::endl;
         std::cout << "PISTON HEIGHT " << pistonHeight << std::endl;

         // this is the flat wall
         compressionPiston.setSpecies(speciesWall);
         compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
         pistonPointer = wallHandler.copyAndAddObject(compressionPiston);

         std::cout << "PISTON CREATED; COMPOSED OF " << pistonLatticeParticleCounter << " PARTICLES." << std::endl;
         std::cout << "NUMBER OF FREE PARTICLES " << nOfFreeParticles << std::endl;
         std::cout << "REALITY CHECK: " << particleHandler.getNumberOfObjects() << " = " << nOfFreeParticles + pistonLatticeParticleCounter << " ?" << std::endl	<< std::endl;
      }
   }

   // creates the bladed piston
   void makeComplexPiston()
   {
      std::cout << std::endl << "CREATING DENTED PISTON" << std::endl << std::endl;

      dentsHeight = 8.0e-4;
      dentsThickness = radiusBig/1000.0;
      dentLengthOpeningRatio = 0.25;

      pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius()) + dentsHeight;
      pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

      compressionPiston.setSpecies(speciesBig);
      compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(compressionPiston);

      // wall x=0 [+]
      w1.setSpecies(speciesBig);
      w1.addObject(Vec3D(1.0,0.0,0.0),Vec3D(-0.5*dentsThickness,0.0,0.0));
      w1.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w1.addObject(Vec3D(-1.0,0.0,0.0),Vec3D(0.5*dentsThickness,0.0,0.0));
      w1.addObject(Vec3D(0.0,1.0,0.0),Vec3D(0.0,dentLengthOpeningRatio*casingRadius,0.0));
      w1pointer = wallHandler.copyAndAddObject(w1);

      // wall x=0 [-]
      w1B.setSpecies(speciesBig);
      w1B.addObject(Vec3D(1.0,0.0,0.0),Vec3D(-0.5*dentsThickness,0.0,0.0));
      w1B.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w1B.addObject(Vec3D(-1.0,0.0,0.0),Vec3D(0.5*dentsThickness,0.0,0.0));
      w1B.addObject(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,-dentLengthOpeningRatio*casingRadius,0.0));
      w1Bpointer = wallHandler.copyAndAddObject(w1B);

      // wall y=0 [+]
      w2.setSpecies(speciesBig);
      w2.addObject(Vec3D(0.0,1.0,0.0),Vec3D(0.0,-0.5*dentsThickness,0.0));
      w2.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w2.addObject(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,0.5*dentsThickness,0.0));
      w2.addObject(Vec3D(1.0,0.0,0.0),Vec3D(dentLengthOpeningRatio*casingRadius,0.0,0.0));
      w2pointer = wallHandler.copyAndAddObject(w2);

      // wall y=0 [-]
      w2B.setSpecies(speciesBig);
      w2B.addObject(Vec3D(0.0,1.0,0.0),Vec3D(0.0,-0.5*dentsThickness,0.0));
      w2B.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w2B.addObject(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,0.5*dentsThickness,0.0));
      w2B.addObject(Vec3D(-1.0,0.0,0.0),Vec3D(-dentLengthOpeningRatio*casingRadius,0.0,0.0));
      w2Bpointer = wallHandler.copyAndAddObject(w2B);

      // wall y=-x [+]
      w3.setSpecies(speciesBig);
      w3.addObject(Vec3D(1.0,1.0,0.0),Vec3D(-0.5*dentsThickness*constants::sqrt_2,0.0,0.0));
      w3.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w3.addObject(Vec3D(-1.0,-1.0,0.0),Vec3D(0.0,0.5*dentsThickness*constants::sqrt_2,0.0));
      w3.addObject(Vec3D(-1.0,1.0,0.0),Vec3D(-dentLengthOpeningRatio*casingRadius/constants::sqrt_2,dentLengthOpeningRatio*casingRadius/constants::sqrt_2,0.0));
      w3pointer = wallHandler.copyAndAddObject(w3);

      // wall y=-x [-]
      w3B.setSpecies(speciesBig);
      w3B.addObject(Vec3D(1.0,1.0,0.0),Vec3D(-0.5*dentsThickness*constants::sqrt_2,0.0,0.0));
      w3B.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w3B.addObject(Vec3D(-1.0,-1.0,0.0),Vec3D(0.0,0.5*dentsThickness*constants::sqrt_2,0.0));
      w3B.addObject(Vec3D(1.0,-1.0,0.0),Vec3D(dentLengthOpeningRatio*casingRadius/constants::sqrt_2,-dentLengthOpeningRatio*casingRadius/constants::sqrt_2,0.0));
      w3Bpointer = wallHandler.copyAndAddObject(w3B);

      // wall y=x [+]
      w4.setSpecies(speciesBig);
      w4.addObject(Vec3D(-1.0,1.0,0.0),Vec3D(0.5*dentsThickness*constants::sqrt_2,0.0,0.0));
      w4.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w4.addObject(Vec3D(1.0,-1.0,0.0),Vec3D(0.0,0.5*dentsThickness*constants::sqrt_2,0.0));
      w4.addObject(Vec3D(1.0,1.0,0.0),Vec3D(dentLengthOpeningRatio*casingRadius/constants::sqrt_2,dentLengthOpeningRatio*casingRadius/constants::sqrt_2,0.0));
      w4pointer = wallHandler.copyAndAddObject(w4);

      // wall y=x [-]
      w4B.setSpecies(speciesBig);
      w4B.addObject(Vec3D(-1.0,1.0,0.0),Vec3D(0.5*dentsThickness*constants::sqrt_2,0.0,0.0));
      w4B.addObject(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight - dentsHeight));
      w4B.addObject(Vec3D(1.0,-1.0,0.0),Vec3D(0.0,0.5*dentsThickness*constants::sqrt_2,0.0));
      w4B.addObject(Vec3D(1.0,1.0,0.0),Vec3D(-dentLengthOpeningRatio*casingRadius/constants::sqrt_2,-dentLengthOpeningRatio*casingRadius/constants::sqrt_2,0.0));
      w4Bpointer = wallHandler.copyAndAddObject(w4B);
   }

   // displaces the dented piston
   void moveComplexPiston()
   {
      pistonHeight += pistonVelocity*getTimeStep();

      w1pointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w1Bpointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w2pointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w2Bpointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w3pointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w3Bpointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w4pointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));
      w4Bpointer -> move(Vec3D(0.0,0.0,pistonVelocity*getTimeStep()));

      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }

   // analyzes simulations data of interest
   void makeDataAnalysis()
   {
      meanCoordinationNumber = 0.0;

      for (int i=nOfFreeParticles-1; i>=0; i--)
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

         // piston interactions -- fixed particles
         if ((*i) -> getI() -> getIndex() > nOfFreeParticles - 1)
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;
            //pistonTorque += ((*i) -> getTorque()).Z;

            //if (fabs(((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X - ((*i) -> getTorque()).Z)/(((*i) -> getTorque()).Z) > 1.0e-3) std::cout <<  (*i) -> getI() -> getIndex() << "\t";

            pistonInteractionCounter++;
         }

         // piston interactions -- flat wall
         if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         // ----------
         // piston interactions -- dents
         if ((*i) -> getI() -> getIndex() == w1pointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w1Bpointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w2pointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w2Bpointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w3pointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w3Bpointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w4pointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }

         if ((*i) -> getI() -> getIndex() == w4Bpointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;

            meanPistonRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxPistonRelativeOverlap) maxPistonRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

            pistonInteractionCounter++;
         }
         // ----------

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

   // sets the next data point to calibrate, and resets the height/pressure/torque flags
   void updateTargetDataPoint()
   {
      targetPressure = pressureSteps[calibrationStep];
      targetTorque = torqueSteps[calibrationStep];

      std::cout << "NEXT PISTON TARGET PRESSURE: " << targetPressure << std::endl;
      std::cout << "NEXT PISTON TARGET TORQUE: " << targetTorque << std::endl << std::endl;

      isPressureTargetMet = false;
      isTorqueTargetMet = false;
   }

   // resets the piston velocity to the starting value
   void resetPistonVelocity()
   {
      if (pistonPressure > targetPressure) pistonVelocity = pistonMaximumVelocity;
      else pistonVelocity = -pistonMaximumVelocity;

      std::cout << "RESETTING PISTON VELOCITY. NEW VELOCITY: " << pistonVelocity << std::endl << std::endl;
   }

   // displaces the piston position according to its velocity
   void movePiston()
   {
      pistonHeight += pistonVelocity*getTimeStep();

      if (!isPistonFlat)
      {
         Vec3D newLatticeParticlePosition;
         effectivePistonHeight += pistonVelocity*getTimeStep();

         for (int i=particleHandler.getNumberOfObjects()-1; i>particleHandler.getNumberOfObjects()-pistonLatticeParticleCounter-1; i--)
         {
            newLatticeParticlePosition = particleHandler.getObject(i) -> getPosition();
            newLatticeParticlePosition.Z += pistonVelocity*getTimeStep();
            particleHandler.getObject(i) -> setPosition(newLatticeParticlePosition);
            particleHandler.getObject(i) -> setVelocity(Vec3D(0.0,0.0,0.0));
         }
      }

      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      // pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));
   }

   // rotates the piston
   void rotatePiston() // check
   {
      if (!isPistonFlat)
      {
         Vec3D newLatticeParticlePosition;
         Vec3D oldLatticeParticlePosition;

         double angularDisplacement = pistonAngularVelocity*getTimeStep();

         for (int i=particleHandler.getNumberOfObjects()-1; i>particleHandler.getNumberOfObjects()-pistonLatticeParticleCounter-1; i--)
         {
            oldLatticeParticlePosition = particleHandler.getObject(i) -> getPosition();

            newLatticeParticlePosition.X = oldLatticeParticlePosition.X*cos(angularDisplacement) - oldLatticeParticlePosition.Y*sin(angularDisplacement);
            newLatticeParticlePosition.Y = oldLatticeParticlePosition.X*sin(angularDisplacement) + oldLatticeParticlePosition.Y*cos(angularDisplacement);
            newLatticeParticlePosition.Z = oldLatticeParticlePosition.Z;

            particleHandler.getObject(i) -> setPosition(newLatticeParticlePosition);
         }
      }

      pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,pistonAngularVelocity));
      pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));

      w1pointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w1Bpointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w2pointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w2Bpointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w3pointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w3Bpointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w4pointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
      w4Bpointer -> rotateAroundZ(pistonAngularVelocity*(getTimeStep()-t0));
   }

   // wiggles the piston to pre-shear the particle bed
   void wigglePiston(double t0) // check
   {
      double wigglingAngularVelocity = preShearingOscillationFrequency*preShearingOscillationAmplitude*cos(2.0*constants::pi*preShearingOscillationFrequency*(getTime()-t0));

      if (!isPistonFlat)
      {
         double angularDisplacement = wigglingAngularVelocity*getTimeStep();

         // t0 is the time the piston starts wiggling
         Vec3D newLatticeParticlePosition;
         Vec3D oldLatticeParticlePosition;

         for (int i=particleHandler.getNumberOfObjects()-1; i>particleHandler.getNumberOfObjects()-pistonLatticeParticleCounter-1; i--)
         {
            oldLatticeParticlePosition = particleHandler.getObject(i) -> getPosition();

            newLatticeParticlePosition.X = oldLatticeParticlePosition.X*cos(angularDisplacement) - oldLatticeParticlePosition.Y*sin(angularDisplacement);
            newLatticeParticlePosition.Y = oldLatticeParticlePosition.X*sin(angularDisplacement) + oldLatticeParticlePosition.Y*cos(angularDisplacement);
            newLatticeParticlePosition.Z = oldLatticeParticlePosition.Z;

            particleHandler.getObject(i) -> setPosition(newLatticeParticlePosition);
         }
      }

      pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,wigglingAngularVelocity));
      pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));

      w1pointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w1Bpointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w2pointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w2Bpointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w3pointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w3Bpointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w4pointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
      w4Bpointer -> rotateAroundZ(wigglingAngularVelocity*(getTimeStep()-t0));
   }

   // resets the species for the big particles, according to the new cohesion stiffness
   void refreshSpeciesBig()
   {
      // BIG-BIG
      speciesBig -> setCohesionStiffness(kCBig);

      // BIG-WALL
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kCBig);

      // BIG-SMALL
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kCBig + kCSmall));
   }

   // adjusts the cohesion stiffness of big particles according to the piston torque
   void tuneKCBig()
   {
      if (pistonTorque < targetTorque) kCBig *= 1.001;
      else kCBig *= 0.999;

      refreshSpeciesBig();

      if (kCBig < 0.0)
      {
         std::cout << "\n\nFATAL ERROR: KC_BIG < 0.0. QUITTING.";
         setTimeMax(getTime() + 10.0*getTimeStep());
         stage = 10;
      }
   }

   // the pre-shearing loop
   void runPreShearing()
   {
      // upon reaching the target height update the bool and sets the velocity to 0
      if (!isPressureTargetMet)
      {
         if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
         {
            if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -0.5;
            // movePiston();
            moveComplexPiston();
         }
         else
         {
            std::cout << "Pressure target reached.\n";
            std::cout << "Starting wiggling...\n";
            t0 = getTime();

            pistonVelocity = 0.0;
            resetPistonVelocity();
            isPressureTargetMet = true;
         }
      }

      // when in position wiggle and keep the pressure constant
      if (isPressureTargetMet && getTime() < t0 + preShearingDuration)
      {
         wigglePiston(t0);
      }

      // after the pre-shearing phase go on with teh calibration
      if (getTime() > t0 + preShearingDuration)
      {
         std::cout << "Pre-shearing complete...\n";
         std::cout << "Starting the torque calibration...\n";
         isPressureTargetMet = false;
         resetPistonVelocity();
         stage++;
      }
   }

   // the torque calibration loop
   void runTorqueCalibration() // check
   {
      if (!isPressureTargetMet && !isTorqueTargetMet)
      {
         if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
         {
            if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -0.5;
            movePiston();
         }
         else
         {
            std::cout << "Pressure target reached, now rotating the piston\n";
            t0 = getTime();
            isPressureTargetMet = true;
         }
      }

      // when the pressure is met start rotating
      if (isPressureTargetMet && !isTorqueTargetMet)
      {
         if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
         {
            if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -1.0;
            movePiston();
         }
         rotatePiston();
      }

      // after the piston rotated for a while deal with the torque
      if (isPressureTargetMet && !isTorqueTargetMet && getTime() > t0 + pauseBeforeTorqueCalibration)
      {
         if (fabs((pistonTorque - targetTorque)/targetTorque) > maxTorqueRelDeviation)
         {
            // the stiffness is changed every n time steps
            if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneKCBig();
         }
         else
         {
            std::cout << "Torque target reached, now chillin' for " << pauseAfterDataPointMet << " seconds\n";
            t0 = getTime();
            isTorqueTargetMet = true;
         }

         // sanity check to keep kC > 1
         if (kCBig < 1.0 && targetTorque < pistonTorque)
         {
            std::cout << "COHESION STIFFNESS FALLING TOO LOW, SWITCHING TO NEXT CALIBRATION POINT." << std::endl;
            std::cout << "Chillin' for " << pauseAfterDataPointMet << " seconds" << std::endl;
            t0 = getTime();
            isTorqueTargetMet = true;
         }

         // sanity check to keep kC < kP
         if (kCBig > k1Big && targetTorque > pistonTorque)
         {
            std::cout << "COHESION STIFFNESS INCREASING TOO MUCH, SWITCHING TO NEXT CALIBRATION POINT." << std::endl;
            std::cout << "Chillin' for " << pauseAfterDataPointMet << " seconds" << std::endl;
            t0 = getTime();
            isTorqueTargetMet = true;
         }
      }

      // switch to next calibration point if the test of the previous is passed
      if (isPressureTargetMet && isTorqueTargetMet && getTime() > t0 + pauseAfterDataPointMet)
      {
         // updates the mean value of the stiffnesses at the calibrated data points
         kCCalibratedArray[calibrationStep] = kCBig;
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

            makeKdataFile();
            if (nCurrentRoutineRun < nMaximumRoutineRuns) updateInitialKArray();

            setTimeMax(getTime() + getTimeStep());
            stage = 10;
         }
      }
   }

   // computes averages, deviations and displacements of the calibrated values of KC
   void computeVariations()
   {
      meanKCBig = 0.0;
      meanRelativeDeviationKCBig = 0.0;
      meanRelativeDisplacementeKCBig = 0.0;

      for (int i = 0; i < nCalibrationDataPoints; i++) meanKCBig += kCCalibratedArray[i];
      meanKCBig /= nCalibrationDataPoints;

      for (int i = 0; i < nCalibrationDataPoints; i++) meanRelativeDeviationKCBig += fabs(kCCalibratedArray[i] - meanKCBig);
      meanRelativeDeviationKCBig /= nCalibrationDataPoints*meanKCBig;

      meanRelativeDisplacementeKCBig = fabs(meanKCBig - pointerToInitialKCArray[nCurrentRoutineRun])/pointerToInitialKCArray[nCurrentRoutineRun];
   }

   // updates the external array of the deviations
   void updateVariationArrays()
   {
      pointerToDeviationArray[nCurrentRoutineRun] = meanRelativeDeviationKCBig;
      pointerToDisplacementArray[nCurrentRoutineRun] = meanRelativeDisplacementeKCBig;
   }

   // updates the initial values of the next calibration run
   void updateInitialKArray()
   {
      pointerToInitialKCArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialKCArray[nCurrentRoutineRun] + meanKCBig);
   }

   // creates the kdata output file
   void makeKdataFile()
   {
      std::ostringstream kdataName;
      std::cout.unsetf(std::ios::floatfield);
      kdataName << getName() << ".kdata";

      kdataFile.open(kdataName.str(), std::ios::out);
      kdataFile << "k1 = " << k1Big << std::endl;
      kdataFile << "k2_max = " << k2MaxBig << std::endl;
      kdataFile << "kC_init = " << pointerToInitialKCArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "phi = " << phiBig << std::endl;
      kdataFile << "|kC| = " << meanKCBig << std::endl;
      kdataFile << "MD(kC)/|kC| = " << pointerToDeviationArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "||kC| - kC_init|/kC_init = " << pointerToDisplacementArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "height \t torque  \t kC" << std::endl;
      for (int i = 0; i < nCalibrationDataPoints; i++)
      {
         kdataFile << pressureSteps[i] << "\t" << torqueSteps[i] << "\t" << kCCalibratedArray[i] << std::endl;
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
      cdatFile << "time \t Eel \t Ekin/Eel \t h \t v \t omega \t P \t Pbase \t Pcyl \t torque \t torque_target \t |(T - Ttarg)/Ttarg| \t dTotMean \t dTotMax \t dPistonMean \t dPistonMax \t dBaseMean \t dBaseMax \t k1 \t k2 \t kC \t phi" << std::endl;
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
      pistonAngularVelocity << "   " <<
      pistonPressure << "   " <<
      basePressure << "   " <<
      cylinderPressure << "   " <<
      pistonTorque << "   " <<
      targetTorque << "   " <<
      fabs((pistonTorque - targetTorque)/targetTorque) << "   " <<
      meanTotalRelativeOverlap << "   " <<
      maxTotalRelativeOverlap << "   " <<
      meanPistonRelativeOverlap << "   " <<
      maxPistonRelativeOverlap << "   " <<
      meanBaseRelativeOverlap << "   " <<
      maxBaseRelativeOverlap << "   " <<
      k1Big << "   " <<
      k2MaxBig << "   " <<
      kCBig << "   " <<
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
         ", v = " << pistonVelocity << ", omega = " << pistonAngularVelocity << std::endl <<
         "P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << ", torque = " << pistonTorque << " , torque_target = " << targetTorque << std::endl <<
         "cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
         ", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap <<
         ", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
         "k1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() <<
         ", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kC = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl <<
         "PRESSURE MET: " << isPressureTargetMet << ", TORQUE MET: " << isTorqueTargetMet << std::endl << std::endl;
         std::cout.flush();
      }
   }


   // VARIABLES   -----------------------------------------------------
   // particles
   double radiusBig, radiusSmall;
   double sizeDispersityBig, sizeDispersitySmall;
   double densityBig, densitySmall;
   double volumeBig, volumeSmall;
   double massBig, massSmall;
   double smallToBigMassRatio;
   int nOfFreeParticles;

   // powder bed
   double particleBedHeight;

   // geometry
   double casingRadius;
   double densityWall;

   // piston
   double pistonHeight;
   double effectivePistonHeight;
   double pistonVelocity;
   double pistonVelocityScalingFactor;
   double pistonMaximumVelocity;
   double pistonForce;
   double pistonPressure;
   double pistonAngularVelocity;
   double pistonTorque;
   InfiniteWall compressionPiston;
   InfiniteWall* pistonPointer;
   int pistonLatticeParticleCounter;
   bool isPistonFlat;

   // complex piston
   RotatingIntersectionOfWalls w1;
   RotatingIntersectionOfWalls* w1pointer;
   RotatingIntersectionOfWalls w2;
   RotatingIntersectionOfWalls* w2pointer;
   RotatingIntersectionOfWalls w1B;
   RotatingIntersectionOfWalls* w1Bpointer;
   RotatingIntersectionOfWalls w2B;
   RotatingIntersectionOfWalls* w2Bpointer;
   RotatingIntersectionOfWalls w3;
   RotatingIntersectionOfWalls* w3pointer;
   RotatingIntersectionOfWalls w4;
   RotatingIntersectionOfWalls* w4pointer;
   RotatingIntersectionOfWalls w3B;
   RotatingIntersectionOfWalls* w3Bpointer;
   RotatingIntersectionOfWalls w4B;
   RotatingIntersectionOfWalls* w4Bpointer;
   double dentsHeight;
   double dentsThickness;
   double dentLengthOpeningRatio;

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
   double pauseBeforeTorqueCalibration;
   double *kCCalibratedArray;
   double *pressureSteps;
   double *torqueSteps;
   bool isPressureTargetMet;
   bool isTorqueTargetMet;
   int calibrationStep;
   double targetPressure;
   double targetTorque;
   double maxTorqueRelDeviation;
   double maxPressureRelDeviation;
   double nTimeStepsBetweenStiffnessAdjustments;
   double *pointerToInitialKCArray;
   double *pointerToDeviationArray;
   double *pointerToDisplacementArray;
   double meanKCBig;
   double meanRelativeDeviationKCBig;
   double meanRelativeDisplacementeKCBig;

   // pre-shearing
   double preShearingDuration;
   double preShearingOscillationAmplitude;
   double preShearingOscillationFrequency;

   // global
   int stage;
   double t0;
   int nMaximumRoutineRuns;
   int nCurrentRoutineRun;
};

// reads the experiment data from file and initializes the target values for the calibration
void readExperimentData(int nTot, std::string inputFileName, double *hArray, double *pArray, double *tArray, double h0)
{
   std::ifstream inputFile(inputFileName);
   if (inputFile.is_open())
   {
      for (int i = 0; i < nTot; i++)
      {
         inputFile >> hArray[i];
         inputFile >> pArray[i];
         inputFile >> tArray[i];
         hArray[i] = h0 - hArray[i];
      }

      inputFile.close();
   }
   else
   {
      std::cout << "Unable to open SHEAR data file" << std::endl;
      std::cout << "QUITTING..." << std::endl;
      exit(1);
   }
}

// initializes the calibration result arrays
void initializeCalibrationArrays(int nRunsMax, double *initKC, double *devArray, double *dispArray, double kCBigInit)
{
   for (int i = 0; i < nRunsMax; i++)
   {
      initKC[i] = 0.0;
      devArray[i] = 1.0;
      dispArray[i] = 1.0;
   }

   initKC[0] = kCBigInit;
}

// computes the delta intersection
double deltaIntersection(double d1, double d2, double p1, double p2)
{
   return d1 - p1*(d2 - d1)/(p2 - p1);
}

// computes the deltaHatMax and the cConstant from the COMPRESSION calibration data
void computeDeltaHatMaxAndCConstant(int nTot, std::string inputFileName, double h0, double *dMax, double *cC)
{
   int indexOfHighestPressurePoint;
   double highestPressure = 0.0;
   double *dTempArray;
   double *pTempArray;
   dTempArray = new double[nTot];
   pTempArray = new double[nTot];

   // allocates the temporary pressure and height array from the compression test
   std::ifstream inputFile(inputFileName);
   if (inputFile.is_open())
   {
      for (int i = 0; i < nTot; i++)
      {
         inputFile >> dTempArray[i];
         inputFile >> pTempArray[i];
      }

      inputFile.close();
   }
   else
   {
      std::cout << "Unable to open COMPRESSION data file" << std::endl;
      std::cout << "QUITTING..." << std::endl;
      exit(1);
   }

   // finds the highest pressure point
   for (int i = 0; i < nTot; i++)
   {
      if (pTempArray[i] > highestPressure)
      {
         highestPressure = pTempArray[i];
         indexOfHighestPressurePoint = i;
      }
   }

   // computes the deltaHatMax
   dMax[0] = (dTempArray[indexOfHighestPressurePoint])/h0;

   // computes the cConstant
   cC[0] = (dTempArray[indexOfHighestPressurePoint])/
   (dTempArray[indexOfHighestPressurePoint] - deltaIntersection(dTempArray[nTot - 1], dTempArray[indexOfHighestPressurePoint], pTempArray[nTot - 1], pTempArray[indexOfHighestPressurePoint]));

   // cC[0] = (dTempArray[indexOfHighestPressurePoint] - deltaIntersection(dTempArray[0], dTempArray[indexOfHighestPressurePoint], pTempArray[0], pTempArray[indexOfHighestPressurePoint]))/
   // (dTempArray[indexOfHighestPressurePoint] - deltaIntersection(dTempArray[nTot - 1], dTempArray[indexOfHighestPressurePoint], pTempArray[nTot - 1], pTempArray[indexOfHighestPressurePoint]));

   std::cout << "dMaxHat " << dMax[0] << std::endl;
   std::cout << "cC " << cC[0] << std::endl;

   delete [] dTempArray;
   delete [] pTempArray;
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
   std::string material = "CAPSULAC";
   std::string dataPointsFileNameCOMPRESSION = "Calibration_NEWDATASETS/Calibration_dataPoints_COMPRESSION_CAPSULAC60___RESCALED.dat";
   std::string dataPointsFileNameSHEAR = "Calibration_NEWDATASETS/Calibration_dataPoints_SHEAR_CAPSULAC60.dat";

   // PISTON TEXTURE SELECTION
   bool isPistonFlat = true;

   // DATA POINT SELECTION FOR ANALYSIS
   int nDataPointsToAnalyze = 6; // ALWAYS SKIP THE FINAL UNLOADING DATA POINT
   int nTotalDataPoints = 16;
   double pauseDuration = 0.01;
   double pauseBeforeTorqueCalibration = 0.08;
   double heightArray[nTotalDataPoints];
   double pressureArray[nTotalDataPoints];
   double torqueArray[nTotalDataPoints];

   // TIME STEP
   double timeStep = 1.0e-6;

   // GEOMETRY PARAMETERS
   double particleBedHeight = 0.019;
   double casingRadius = 0.0125;
   double densityWalls = 5000.0;

   // SHEARING PARAMETERS
   double pistonVelocityScalingFactor = 2.0;
   double pistonRealVelocity = 0.3*constants::pi/180.0;
   double pistonAngularVelocity = 20.0*pistonRealVelocity;

   double preShearingDuration = 0.5;
   double preShearingOscillationAmplitude = constants::pi/12.0; // this has to be translated in radians
   double preShearingOscillationFrequency = 5.0/preShearingDuration;

   // PARTICLE STATIC PARAMETERS
   double bigToCasingSizeRatio = 40.0;
   double bigToSmallSizeRatio = 1.0;
   double smallToBigMassRatio = 0.0;
   double sizeDispersionBigParticles = 0.1;
   double sizeDispersionSmallParticles = 0.1;
   double densityBigParticles = 1500.0;
   double densitySmallParticles = 1500.0;

   // COMPUTATION OF dHatMax AND cConstant
   double deltaHatMax[1];
   double cConstant[1];
   computeDeltaHatMaxAndCConstant(nTotalDataPoints, dataPointsFileNameCOMPRESSION, particleBedHeight, deltaHatMax, cConstant);

   // INTERACTION PARAMETERS
   double k2MaxBig = 3000.0;
   double k1Big = 1086.0;
   double kCBig = k2MaxBig/1000.0;
   double phiBig = pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax[0]/(cConstant[0] - 1.0);

   double k2MaxSmall = 3000.0;
   double k1Small = k2MaxSmall/10.0;
   double kCSmall = k2MaxSmall/1000.0;
   double phiSmall = pow(1.0 - k1Small/k2MaxSmall, 2.0)/(k1Small/k2MaxSmall)*deltaHatMax[0]/(cConstant[0] - 1.0);

   double k1Wall = k2MaxBig;

   double eBigWall = 0.7;
   double eSmallWall = 0.7;
   double eBigBig = 0.5;
   double eSmallSmall = 0.5;
   double eBigSmall = 0.5;

   // FRICTION PARAMETERS
   double muSBigBig = 1.00;
   double muRBigBig = 0.10;
   double muSBigWall = muSBigBig;
   double muRBigWall = muRBigBig;

   // CALIBRATION PARAMETERS
   double maximumStiffnessRelativeDeviation = 0.05;
   double maximumPressureRelativeDeviation = 1.0e-4;
   double maximumTorqueRelativeDeviation = 0.05;
   double timeStepsBetweenStiffnessTunings = 10.0;

   // CALIBRATION ARRAYS
   int runNumber = 0;
   int nMaximumRoutineCycles = 1;
   double initialKCArray[nMaximumRoutineCycles];
   double deviationArray[nMaximumRoutineCycles];
   double displacementArray[nMaximumRoutineCycles];

   // CALIBRATION DATA POINTS AND ARRAYS INITIALIZATION
   readExperimentData(nTotalDataPoints, dataPointsFileNameSHEAR, heightArray, pressureArray, torqueArray, particleBedHeight);
   initializeCalibrationArrays(nMaximumRoutineCycles, initialKCArray, deviationArray, displacementArray, kCBig);

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

   std::string phase;
   if (nDataPointsToAnalyze == 6) phase = "LOADING";
   else if (nDataPointsToAnalyze == 10) phase = "RELOADING";
   else if (nDataPointsToAnalyze == nTotalDataPoints - 1) phase = "FULL";
   else phase = "";

   std::string pistonType;
   if (isPistonFlat) pistonType = "flat";
   else pistonType = "rough";

   // TIME WARNINGS
   std::cout << "Minimum pause duration before torque calibration: " << 3.0/(2.0*bigToCasingSizeRatio*pistonAngularVelocity) << std::endl;
   std::cout << "Pause before torque calibration selected: " << pauseBeforeTorqueCalibration << std::endl;

   // THE CALIBRATION LOOP
   do
   {
      ShearTest_calibrationRoutine crShear;
      crShear.readArguments(argc, argv);

      // NAME SETTING
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      name << "ShearTest_" << material << "_dentedPistonTest_8walls_muS_1.00";
      crShear.setName(name.str());

      std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl << std::endl;

      // INITIALIZATION
      crShear.setTimeStep(timeStep);
      crShear.setTimeMax(1.0);
      crShear.setGravity(Vec3D(0.00,0.00,-9.81));
      crShear.setSystemDimensions(3);
      //crShear.setVerbose(true);
      crShear.setCdatOutputTimeInterval(0.001);
      crShear.setSaveCount(0.01/crShear.getTimeStep());

      crShear.setCasingProperties(casingRadius, particleBedHeight, densityWalls);
      crShear.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);

      crShear.setParticleWallSlidingFrictionCoefficients(muSBigWall, muSBigWall);
      crShear.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
      crShear.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);

      crShear.setParticleParticleSlidingFrictionCoefficients(muSBigBig, muSBigBig, muSBigBig);
      crShear.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
      crShear.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);

      crShear.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
      crShear.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
      crShear.setBigParticlePlasticProperties(k1Big, k2MaxBig, initialKCArray[runNumber], phiBig);
      crShear.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);

      crShear.setTorqueCalibrationPoints(nDataPointsToAnalyze, pressureArray, torqueArray, pauseDuration, pauseBeforeTorqueCalibration);
      crShear.setCalibrationMaximumRelativeDeviations(maximumPressureRelativeDeviation, maximumTorqueRelativeDeviation);
      crShear.setCalibrationRoutineParameters(runNumber, nMaximumRoutineCycles, initialKCArray, deviationArray, displacementArray);
      crShear.setNTimeStepsBetweenStiffnesAdjustments(timeStepsBetweenStiffnessTunings);

      crShear.setPistonType(isPistonFlat);
      crShear.setPistonRotation(pistonAngularVelocity);
      crShear.setPistonVelocityScalingFactor(pistonVelocityScalingFactor);
      crShear.setPreShearingParameters(preShearingDuration, preShearingOscillationAmplitude, preShearingOscillationFrequency);

      crShear.solve();
      runNumber++;
   }
   while (runNumber < nMaximumRoutineCycles && !(displacementArray[runNumber - 1] < maximumStiffnessRelativeDeviation && deviationArray[runNumber - 1] < maximumStiffnessRelativeDeviation));

   return 0;
}
