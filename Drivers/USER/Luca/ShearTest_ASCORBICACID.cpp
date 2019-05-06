#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "CompressionPistonSurface.h"
#include <math.h>
#include <fstream>
#include <vector>
#include <numeric>

// VERSION: ShearTest 30.01.19 (pressure constrained)

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
      pistonAngle = 0.0;
      meanMuSBig = 0.0;
      meanKCBig = 0.0;
      isPressureTargetMet = false;
      isTorqueTargetMet = false;

      // species set and reset
      setSpecies();
      resetSpecies();

      // setup of all the rest
      makePiston();
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
      if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep())
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() <<
         ", h/h0 = " << pistonHeight/0.019 <<
         ", vZ = " << pistonVelocity <<
         ", theta = " << pistonAngle <<
         ", P = " << std::setprecision(6) << std::left << pistonPressure <<
         ", P_target = " << targetPressure <<
         ", T = " << pistonTorque <<
         ", T_target = " << targetTorque <<
         ", mu_S = " << speciesBig -> getSlidingFrictionCoefficient() <<
         ", kC = " << speciesBig -> getCohesionStiffness() <<
         std::endl;
      }
   }

   void actionsAfterSolve() override
   {
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

      for (int i = 0; i < nCalibrationDataPoints; i++)
      {
         pressureSteps.push_back(pressureArrayPointer[i]);
         torqueSteps.push_back(torqueArrayPointer[i]);
      }
   }

   void setCalibrationMaximumRelativeDeviations(double maxP, double maxT)
   {
      maxPressureRelDeviation = maxP;
      maxTorqueRelDeviation = maxT;
   }

   void setCalibrationRoutineParameters(int nRun, int nMax, double *muSArray, double *devMuSArray, double *dispMuSArray, double *kCArray, double *devKCArray, double *dispKCArray)
   {
      nMaximumRoutineRuns = nMax;
      nCurrentRoutineRun = nRun;

      pointerToInitialMuSArray = muSArray;
      pointerToMuSDeviationArray = devMuSArray;
      pointerToMuSDisplacementArray = dispMuSArray;

      pointerToInitialKCArray = kCArray;
      pointerToKCDeviationArray = devKCArray;
      pointerToKCDisplacementArray = dispKCArray;
   }

   void setNTimeStepsBetweenStiffnesAdjustments(double n)
   {
      nTimeStepsBetweenParametersAdjustments = n;
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

   void setBladedPistonParameters(int n, double height, double thickness)
   {
      numberOfBlades = n;
      dentsHeight = height;
      dentsThickness = thickness;
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

   // creates the piston
   void makePiston()
   {
      std::cout << "CREATING DENTED PISTON" << std::endl;

      pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius()) + dentsHeight;
      pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

      InfiniteWall compressionPiston;
      compressionPiston.setSpecies(speciesBig);
      compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(compressionPiston);

      CompressionPistonSurface w;
      w.setSpecies(speciesBig);
      for (int n = 0; n < numberOfBlades; n++)
      {
         w.set(pistonHeight, dentsHeight, 0.5*dentsThickness, constants::pi*n/numberOfBlades);
         dentsArray.push_back(wallHandler.copyAndAddObject(w));
      }
   }

   // analyzes simulations data of interest
   void makeDataAnalysis()
   {
      // meanCoordinationNumber = 0.0;
      // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++) meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
      // meanCoordinationNumber /= particleHandler.getNumberOfObjects();

      pistonPressure = 0.0;
      pistonTorque = 0.0;
      basePressure = 0.0;
      cylinderPressure = 0.0;
      meanTotalRelativeOverlap = 0.0;
      maxTotalRelativeOverlap = 0.0;

      int totalInteractionCounter = 0;
      int pistonInteractionCounter = 0;
      int baseInteractionCounter = 0;
      int cylinderInteractionCounter = 0;

      // computing interaction-based data (e.g. pressure, torque, overlaps, etc.)
      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

         // piston interactions -- flat wall
         if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
         {
            pistonPressure -= ((*i) -> getForce()).Z;
            pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;
            pistonInteractionCounter++;
         }

         // piston interactions -- blades
         for (int n = 0; n < numberOfBlades; n++)
         {
            if ((*i) -> getI() -> getIndex() == dentsArray[n] -> getIndex())
            {
               pistonPressure -= ((*i) -> getForce()).Z;
               pistonTorque += ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;
               pistonInteractionCounter++;
            }
         }
         // base interactions
         if ((*i) -> getI() -> getIndex() == 0)
         {
            basePressure += ((*i) -> getForce()).Z;
            baseInteractionCounter++;
         }

         // cylinder interactions
         if ((*i) -> getI() -> getIndex() == 2)
         {
            cylinderPressure -= (((*i) -> getContactPoint()).X * ((*i) -> getForce()).X + ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).Y)/casingRadius;
            cylinderInteractionCounter++;
         }

         totalInteractionCounter++;
      }

      // computes the averaged torque instead of the actual torque
      if (averagedTorque.size() == nTimeStepsBetweenParametersAdjustments) averagedTorque.erase(averagedTorque.begin());
      averagedTorque.push_back(fabs(pistonTorque));
      pistonTorque = std::accumulate(averagedTorque.begin(), averagedTorque.end(), 0.0)/nTimeStepsBetweenParametersAdjustments;

      // rescaling intrinsically averaged variables (e.g. pressure)
      pistonPressure /= constants::pi*pow(casingRadius,2.0);
      basePressure /= constants::pi*pow(casingRadius,2.0);
      if (pistonHeight > particleBedHeight) cylinderPressure /= 2.0*constants::pi*casingRadius*particleBedHeight;
      else cylinderPressure /= 2.0*constants::pi*casingRadius*pistonHeight;
      meanTotalRelativeOverlap /= totalInteractionCounter;
   }

   // sets the next data point to calibrate, and resets the height/pressure/torque flags
   void updateTargetDataPoint()
   {
      // fixes the previous target pressure to see if the piston is compressing or decompressing
      if (calibrationStep != 0) previousTargetPressure = pressureSteps[calibrationStep - 1];
      else previousTargetPressure = 0.0;

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
      for (int n = 0; n < numberOfBlades; n++) dentsArray[n] -> move(pistonVelocity*getTimeStep());
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }

   // rotates the piston
   void rotatePiston()
   {
      pistonAngle += pistonAngularVelocity*getTimeStep();
      for (int n = 0; n < numberOfBlades; n++) dentsArray[n] -> rotate(pistonAngularVelocity*getTimeStep());
      pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,pistonAngularVelocity));
      pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
   }

   // wiggles the piston to pre-shear the particle bed
   void wigglePiston(double t0)
   {
      double wigglingAngularVelocity = preShearingOscillationFrequency*preShearingOscillationAmplitude*std::sin(2.0*constants::pi*(getTime() - t0)*preShearingOscillationFrequency);

      pistonAngle += wigglingAngularVelocity*getTimeStep();
      for (int n = 0; n < numberOfBlades; n++) dentsArray[n] -> rotate(wigglingAngularVelocity*getTimeStep());
      pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,wigglingAngularVelocity));
      pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
   }

   // resets the species for the big particles, according to the new cohesion stiffness
   void refreshSpeciesBig()
   {
      // BIG-BIG
      speciesBig -> setCohesionStiffness(kCBig);
      speciesBig -> setSlidingFrictionCoefficient(bigBigSlidingFrictionCoeff);
      speciesBig -> setSlidingStiffness(speciesBig -> getLoadingStiffness()*2.0/7.0);
      speciesBig -> setSlidingDissipation(speciesBig -> getDissipation()*2.0/7.0);

      // BIG-WALL
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kCBig);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingFrictionCoefficient(bigWallSlidingFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);

      // BIG-SMALL
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kCBig + kCSmall));
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingFrictionCoefficient(bigSmallSlidingFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);
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
         setTimeMax(getTime() + getTimeStep());
         stage = 10;
      }
   }

   // adjusts the sliding friction coefficient of big particles according to the piston torque
   void tuneMUsBig()
   {
      if (pistonTorque < targetTorque) bigBigSlidingFrictionCoeff *= 1.001;
      else bigBigSlidingFrictionCoeff *= 0.999;

      refreshSpeciesBig();

      if (bigBigSlidingFrictionCoeff < 0.0)
      {
         std::cout << "\n\nFATAL ERROR: muS_BIG < 0.0. QUITTING.";
         setTimeMax(getTime() + getTimeStep());
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
            movePiston();
         }
         else
         {
            std::cout << "Height target reached.\n";
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
   void runTorqueCalibration()
   {
      if (!isPressureTargetMet && !isTorqueTargetMet)
      {
         if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
         {
            if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -0.99;
            movePiston();
            wigglePiston(t0);
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
            // stiffness or friction is changed every n time steps
            if (fmod(getTime() - t0, nTimeStepsBetweenParametersAdjustments*getTimeStep()) < getTimeStep())
            {
               // if the piston is in compression tune the friction coefficient, tine the cohesion stiffness otherwise
               if (previousTargetPressure < targetPressure)
               {
                 // sanity check to keep mu_S > 0
                 if (bigBigSlidingFrictionCoeff < 0.001 && targetTorque < pistonTorque)
                 {
                    std::cout << "FRICTION COEFFICIENT FALLING TOO LOW, SWITCHING TO NEXT CALIBRATION POINT." << std::endl;
                    std::cout << "Chillin' for " << pauseAfterDataPointMet << " seconds" << std::endl;
                    t0 = getTime();
                    isTorqueTargetMet = true;
                 }

                 // sanity check to keep mu_S < 0.50
                 if (bigBigSlidingFrictionCoeff > 0.50 && targetTorque > pistonTorque)
                 {
                    std::cout << "FRICTION COEFFICIENT INCREASING TOO MUCH, SWITCHING TO NEXT CALIBRATION POINT." << std::endl;
                    std::cout << "Chillin' for " << pauseAfterDataPointMet << " seconds" << std::endl;
                    t0 = getTime();
                    isTorqueTargetMet = true;
                 }

                  tuneMUsBig();
               }
               else
               {
                  // sanity check to keep kC > 0
                  if (kCBig < 0.001 && targetTorque < pistonTorque)
                  {
                     std::cout << "COHESION STIFFNESS FALLING TOO LOW, SWITCHING TO NEXT CALIBRATION POINT." << std::endl;
                     std::cout << "Chillin' for " << pauseAfterDataPointMet << " seconds" << std::endl;
                     t0 = getTime();
                     isTorqueTargetMet = true;
                  }

                  // sanity check to keep kC < kE
                  if (kCBig > k2MaxBig && targetTorque > pistonTorque)
                  {
                     std::cout << "COHESION STIFFNESS INCREASING TOO MUCH, SWITCHING TO NEXT CALIBRATION POINT." << std::endl;
                     std::cout << "Chillin' for " << pauseAfterDataPointMet << " seconds" << std::endl;
                     t0 = getTime();
                     isTorqueTargetMet = true;
                  }

                  tuneKCBig();
               }
            }
         }
         else
         {
            std::cout << "Torque target reached, now chillin' for " << pauseAfterDataPointMet << " seconds\n";
            t0 = getTime();
            isTorqueTargetMet = true;
         }
      }

      // switch to next calibration point if the test of the previous is passed
      if (isPressureTargetMet && isTorqueTargetMet && getTime() > t0 + pauseAfterDataPointMet)
      {
         // updates the mean value of friction and cohesion at the calibrated data points
         if (previousTargetPressure < targetPressure) muSCalibratedArray.push_back(bigBigSlidingFrictionCoeff);
         else kCCalibratedArray.push_back(kCBig);

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
            if (nCurrentRoutineRun < nMaximumRoutineRuns) updateInitialArrays();

            setTimeMax(getTime() + getTimeStep());
            stage = 10;
         }
      }
   }

   // computes averages, deviations and displacements of the calibrated values of KC
   void computeVariations()
   {
      // compute the friction variations for compression steps
      meanMuSBig = 0.0;
      meanRelativeDeviationMuSBig = 0.0;
      meanRelativeDisplacementMuSBig = 0.0;

      if (muSCalibratedArray.size() > 0)
      {
         meanMuSBig = std::accumulate(muSCalibratedArray.begin(), muSCalibratedArray.end(), 0.0)/muSCalibratedArray.size();
         for (int n = 0; n < muSCalibratedArray.size(); n++) meanRelativeDeviationMuSBig += fabs(muSCalibratedArray[n] - meanMuSBig);
         meanRelativeDeviationMuSBig /= muSCalibratedArray.size()*meanMuSBig;
         meanRelativeDisplacementMuSBig = fabs(meanMuSBig - pointerToInitialMuSArray[nCurrentRoutineRun])/pointerToInitialMuSArray[nCurrentRoutineRun];
      }

      // compute the friction variations if decompression steps were made
      meanKCBig = 0.0;
      meanRelativeDeviationKCBig = 0.0;
      meanRelativeDisplacementKCBig = 0.0;

      if (kCCalibratedArray.size() > 0)
      {
         meanKCBig = std::accumulate(kCCalibratedArray.begin(), kCCalibratedArray.end(), 0.0)/kCCalibratedArray.size();
         for (int n = 0; n < kCCalibratedArray.size(); n++) meanRelativeDeviationKCBig += fabs(kCCalibratedArray[n] - meanKCBig);
         meanRelativeDeviationKCBig /= kCCalibratedArray.size()*meanKCBig;
         meanRelativeDisplacementKCBig = fabs(meanKCBig - pointerToInitialKCArray[nCurrentRoutineRun])/pointerToInitialKCArray[nCurrentRoutineRun];
      }
   }

   // updates the external array of the deviations
   void updateVariationArrays()
   {
      pointerToMuSDeviationArray[nCurrentRoutineRun] = meanRelativeDeviationMuSBig;
      pointerToMuSDisplacementArray[nCurrentRoutineRun] = meanRelativeDisplacementMuSBig;
      pointerToKCDeviationArray[nCurrentRoutineRun] = meanRelativeDeviationKCBig;
      pointerToKCDisplacementArray[nCurrentRoutineRun] = meanRelativeDisplacementKCBig;
   }

   // updates the initial values of the next calibration run
   void updateInitialArrays()
   {
      pointerToInitialMuSArray[nCurrentRoutineRun + 1] = meanMuSBig;
      pointerToInitialKCArray[nCurrentRoutineRun + 1] = meanKCBig;
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
      kdataFile << "phi = " << phiBig << std::endl << std::endl;
      kdataFile << "mu_S_init = " << pointerToInitialMuSArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "|mu_S| = " << meanMuSBig << std::endl;
      kdataFile << "MD(mu_S)/|mu_S| = " << pointerToMuSDeviationArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "||mu_S| - mu_S_init|/mu_S_init = " << pointerToMuSDisplacementArray[nCurrentRoutineRun] << std::endl << std::endl;
      kdataFile << "kC_init = " << pointerToInitialKCArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "|kC| = " << meanKCBig << std::endl;
      kdataFile << "MD(kC)/|kC| = " << pointerToKCDeviationArray[nCurrentRoutineRun] << std::endl;
      kdataFile << "||kC| - kC_init|/kC_init = " << pointerToKCDisplacementArray[nCurrentRoutineRun] << std::endl << std::endl;
      kdataFile << "height\t torque\t mu_S\t kC" << std::endl;
      for (int i = 0; i < nCalibrationDataPoints; i++)
      {
         kdataFile << pressureSteps[i] << "\t" << torqueSteps[i] << "\t" << muSCalibratedArray[i] << "\t" << kCCalibratedArray[i] << std::endl;
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
      cdatFile << "time \t Eel \t Ekin/Eel \t h \t v \t omega \t P \t Pbase \t Pcyl \t torque \t torque_target \t |(T - Ttarg)/Ttarg| \t dTotMean \t dTotMax \t mu_S \t k1 \t k2 \t kC \t phi" << std::endl;
   }

   // writes the compression data to the output file
   void writeToCdatFile()
   {
      cdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      pistonHeight << "   " <<
      pistonAngle << "   " <<
      pistonVelocity << "   " << // this here is left not to change the overall number of columns in the cdat output
      pistonPressure << "   " <<
      basePressure << "   " <<
      cylinderPressure << "   " <<
      pistonTorque << "   " <<
      targetTorque << "   " <<
      fabs((pistonTorque - targetTorque)/targetTorque) << "   " <<
      meanTotalRelativeOverlap << "   " <<
      maxTotalRelativeOverlap << "   " <<
      bigBigSlidingFrictionCoeff << "   " <<
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
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() <<
         ", h = " << pistonHeight << ", theta = " << pistonAngle <<
         "P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << " , pressure_target = " << targetPressure <<
         ", torque = " << pistonTorque << " , torque_target = " << targetTorque << std::endl <<
         ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
         "mu_S = " << bigBigSlidingFrictionCoeff << ", kC = " << speciesBig -> getCohesionStiffness() <<
         "P_MET: " << isPressureTargetMet << ", T_MET: " << isTorqueTargetMet << std::endl << std::endl;
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

   // powder bed
   double particleBedHeight;

   // geometry
   double casingRadius;
   double densityWall;

   // piston
   double pistonHeight;
   double pistonAngle;
   double pistonVelocity;
   double pistonVelocityScalingFactor;
   double pistonMaximumVelocity;
   double pistonPressure;
   double pistonAngularVelocity;
   double pistonTorque;
   std::vector<double> averagedTorque;
   InfiniteWall* pistonPointer;

   // complex piston
   int numberOfBlades;
   double dentsHeight;
   double dentsThickness;
   std::vector<CompressionPistonSurface*> dentsArray;

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
   LinearPlasticViscoelasticFrictionSpecies *speciesBig, *speciesSmall, *speciesWall;

   // output
   double cdatOutputTimeInterval;
   std::ofstream kdataFile;
   std::ofstream cdatFile;

   // data analysis
   double basePressure;
   double cylinderPressure;
   // double meanCoordinationNumber;
   double meanTotalRelativeOverlap;
   double maxTotalRelativeOverlap;

   // calibration
   int nCalibrationDataPoints;
   double pauseAfterDataPointMet;
   double pauseBeforeTorqueCalibration;
   std::vector<double> pressureSteps;
   std::vector<double> torqueSteps;
   std::vector<double> kCCalibratedArray;
   std::vector<double> muSCalibratedArray;
   bool isPressureTargetMet;
   bool isTorqueTargetMet;
   int calibrationStep;
   double previousTargetPressure;
   double targetPressure;
   double targetTorque;
   double maxTorqueRelDeviation;
   double maxPressureRelDeviation;
   double nTimeStepsBetweenParametersAdjustments;
   double *pointerToInitialMuSArray;
   double *pointerToMuSDeviationArray;
   double *pointerToMuSDisplacementArray;
   double *pointerToInitialKCArray;
   double *pointerToKCDeviationArray;
   double *pointerToKCDisplacementArray;
   double meanMuSBig;
   double meanRelativeDeviationMuSBig;
   double meanRelativeDisplacementMuSBig;
   double meanKCBig;
   double meanRelativeDeviationKCBig;
   double meanRelativeDisplacementKCBig;

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
void initializeCalibrationArrays(int nRunsMax, double *initMuS, double *devMuSArray, double *dispMuSArray, double muSInit, double *initKC, double *devKCArray, double *dispKCArray, double kCBigInit)
{
   for (int i = 0; i < nRunsMax; i++)
   {
      initMuS[i] = 0.0;
      devMuSArray[i] = 1.0;
      dispMuSArray[i] = 1.0;

      initKC[i] = 0.0;
      devKCArray[i] = 1.0;
      dispKCArray[i] = 1.0;
   }

   initMuS[0] = muSInit;
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
   std::string material = "ASCORBICACID";
   std::string dataPointsFileNameCOMPRESSION = "Calibration_dataPoints/Calibration_dataPoints_COMPRESSION_ASCORBICACID___RESCALED_05.dat";
   std::string dataPointsFileNameSHEAR = "Calibration_dataPoints/Calibration_dataPoints_SHEAR_ASCORBICACID.dat";
   // std::string material = "ASCORBICACID";
   // std::string dataPointsFileNameCOMPRESSION = "Calibration_dataPoints/Calibration_dataPoints_COMPRESSION_ASCORBICACID___RESCALED_05.dat";
   // std::string dataPointsFileNameSHEAR = "Calibration_dataPoints/Calibration_dataPoints_SHEAR_ASCORBICACID.dat";

   // DATA POINT SELECTION FOR ANALYSIS
   int nDataPointsToAnalyze = 15; // ALWAYS SKIP THE FINAL UNLOADING DATA POINT
   int nTotalDataPoints = 16;
   double heightArray[nTotalDataPoints];
   double pressureArray[nTotalDataPoints];
   double torqueArray[nTotalDataPoints];

   // TIME STEP
   double timeStep = 2.0e-6;

   // GEOMETRY PARAMETERS
   double particleBedHeight = 0.019;
   double casingRadius = 0.0125;
   double densityWalls = 1500.0;

   // SHEARING PARAMETERS
   double pistonVelocityScalingFactor = 10.0;
   double pistonRealVelocity = 0.3*constants::pi/180.0;
   double pistonAngularVelocity = 500.0*pistonRealVelocity;

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

   // PISTON DENTS PARAMETERS
   double numberOfPistonBlades = 2;
   double pistonBladesHeight = 8.0e-4;
   double pistonBladesThickness = casingRadius/bigToCasingSizeRatio/2.0; // 100

   // COMPUTATION OF dHatMax AND cConstant
   double deltaHatMax[1];
   double cConstant[1];
   computeDeltaHatMaxAndCConstant(nTotalDataPoints, dataPointsFileNameCOMPRESSION, particleBedHeight, deltaHatMax, cConstant);

   // INTERACTION PARAMETERS
   double k2MaxBig = 5000.0;
   double k1Big = 31.0;
   double kCBig = 5.0;
   double phiBig = pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax[0]/(cConstant[0] - 1.0);

   // std::cout << "computed phi = " << pow(1.0 - k1Big/k2MaxBig, 2.0)/(k1Big/k2MaxBig)*deltaHatMax[0]/(cConstant[0] - 1.0) << std::endl;
   // std::cout << "set phi = " << phiBig << std::endl;

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
   double muSBigBig = 0.24;
   double muRBigBig = 0.30;
   double muSBigWall = muSBigBig;
   double muRBigWall = muRBigBig;

   // CALIBRATION PARAMETERS
   double maximumParametersRelativeDeviation = 0.05;
   double maximumPressureRelativeDeviation = 1.0e-4;
   double maximumTorqueRelativeDeviation = 1.0e-4;
   double timeStepsBetweenParametersTunings = 10.0;

   // CALIBRATION PARAMETERS
   int runNumber = 0;
   int nMaximumRoutineCycles = 5;
   double pauseDuration = 0.0;
   double pauseBeforeTorqueCalibration = 0.05;
   double initialMuSArray[nMaximumRoutineCycles];
   double deviationMuSArray[nMaximumRoutineCycles];
   double displacementMuSArray[nMaximumRoutineCycles];
   double initialKCArray[nMaximumRoutineCycles];
   double deviationKCArray[nMaximumRoutineCycles];
   double displacementKCArray[nMaximumRoutineCycles];

   // CALIBRATION DATA POINTS AND ARRAYS INITIALIZATION
   readExperimentData(nTotalDataPoints, dataPointsFileNameSHEAR, heightArray, pressureArray, torqueArray, particleBedHeight);
   initializeCalibrationArrays(nMaximumRoutineCycles, initialMuSArray, deviationMuSArray, displacementMuSArray, muSBigBig, initialKCArray, deviationKCArray, displacementKCArray, kCBig);

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

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
      name << "ShearTest_" << material << "_Nrun_" << runNumber <<
      "_k1_" << std::fixed << std::setprecision(0) << k1Big << "_k2max_" << k2MaxBig << "_kc_" << initialKCArray[runNumber] << "_phi_" << std::fixed << std::setprecision(4) << phiBig <<
      "_muS_" << std::fixed << std::setprecision(2) << initialMuSArray[runNumber] << "_muR_" << muRBigBig << "_nD_" << std::fixed << std::setprecision(0) << numberOfPistonBlades;
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

      crShear.setParticleParticleSlidingFrictionCoefficients(initialMuSArray[runNumber], initialMuSArray[runNumber], initialMuSArray[runNumber]);
      crShear.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
      crShear.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);

      crShear.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
      crShear.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
      crShear.setBigParticlePlasticProperties(k1Big, k2MaxBig, initialKCArray[runNumber], phiBig);
      crShear.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);

      crShear.setTorqueCalibrationPoints(nDataPointsToAnalyze, pressureArray, torqueArray, pauseDuration, pauseBeforeTorqueCalibration);
      crShear.setCalibrationMaximumRelativeDeviations(maximumPressureRelativeDeviation, maximumTorqueRelativeDeviation);
      crShear.setCalibrationRoutineParameters(runNumber, nMaximumRoutineCycles, initialMuSArray, deviationMuSArray, displacementMuSArray, initialKCArray, deviationKCArray, displacementKCArray);
      crShear.setNTimeStepsBetweenStiffnesAdjustments(timeStepsBetweenParametersTunings);

      crShear.setBladedPistonParameters(numberOfPistonBlades, pistonBladesHeight, pistonBladesThickness);
      crShear.setPistonRotation(pistonAngularVelocity);
      crShear.setPistonVelocityScalingFactor(pistonVelocityScalingFactor);
      crShear.setPreShearingParameters(preShearingDuration, preShearingOscillationAmplitude, preShearingOscillationFrequency);

      crShear.solve();
      runNumber++;
   }
   while (runNumber < nMaximumRoutineCycles && !(displacementMuSArray[runNumber - 1] < maximumParametersRelativeDeviation && deviationMuSArray[runNumber - 1] < maximumParametersRelativeDeviation && displacementKCArray[runNumber - 1] < maximumParametersRelativeDeviation && deviationKCArray[runNumber - 1] < maximumParametersRelativeDeviation));

   return 0;
}
