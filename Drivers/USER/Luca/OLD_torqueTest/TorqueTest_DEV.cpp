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
* - put an if/else in the main switching petween k1 calibration (compression) and kc (torque)
* - add a version of maxPressureRelDeviation for the torque
* - PUT THE TORQUE CALIBRATION IN ANOTHER DRIVER! DOESN-T BELONG HERE!
*/

// STAGE 10 IS THE "DEAD" STAGE FOR QUITTING AND REBOOTING

// 01.24 - Substituted kNew = kMean to kNew = (kFinal - kInitial)/2 (and included the computation of kMean)
/*
* 01.26 - Changed teh piston velocity and introduced the input parameter setPistonVelocityScalingFactor
* 01.26 - Now the system waits also after heightTargetMet, not only after pressureTargetMet
* 01.26 - changed teh whole tuning routine and implemented the new kc calibration
*
* 02.20 - debugged displacement and deviation checks
* 02.20 - introduced .cdat output and changed kdata output
* 02.20 - introduced phi auto-tuning
* 02.20 - added computeVariations() and updateVariationArrays() functions
*
* 03.12 - debugging and testing of particle piston
* 03.12 - correct the piston height check bearing in mind the piston thickness
*/

// NOW THE PRE-SHEARING TIME IS HARDCODED. FIX THIS. (IT SHOULD BE A FUNCTION OF THE ANGULAR VELOCITY AND THE PRE-SHEARING EXPERIMENTAL ANGLE)

class TorqueTest_calibrationRoutine : public Mercury3D
{

private:

  void setupInitialConditions() override
  {
    std::cout << "THIS IS NOT SUPPOSED TO BE RAN." << std::endl;
    std::cout << "QUITTING..." << std::endl;

    exit(0);
  }

  void actionsOnRestart() override
  {
    stage = 1;
    calibrationStep = 0;
    t0 = getTime();
    setTimeMax(10.0);
    setTimeMax(2.0);

    setSpecies();
    setPistonSpecies();

    isPreShearingDone = false;
    isHeightTargetMet = false;
    isPressureTargetMet = false;
    isTorqueTargetMet = false;

    // std::cout << particleHandler.getNumberOfObjects() << std::endl;
  }

  void actionsAfterTimeStep() override {}

  // void actionsAfterTimeStep() override
  // {
  //   std::cout << "TEST\n";
  //
  //   if (stage == 1)
  //   {
  //     nOfFreeParticles = particleHandler.getNumberOfObjects();
  //
  //     std::cout << nOfFreeParticles << std::endl;
  //     makeSampledParticlePiston();
  //
  //     std::cout << "TEST\n";
  //     meankcBig = 0.0;
  //     setTime(0.0);
  //
  //     // TEST
  //     for (int i=0; i < nCalibrationDataPoints; i++) std::cout << heightSteps[i] << "   " << torqueSteps[i] << std::endl;
  //
  //     updateTargetTorque();
  //     resetPistonVelocity();
  //     makeCdatFile();
  //
  //     t0 = getTime();
  //     stage++;
  //   }
  //
  //   if (stage >= 2)
  //   {
  //     if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) writeToCdatFile();
  //   }
  //
  //   if (stage == 2)
  //   {
  //     //movePiston();
  //
  //     makeDataAnalysis();
  //     runTorqueCalibration();
  //
  //     // TEST DYNAMIC OUTPUT
  //     if (fmod(getTime(), cdatOutputTimeInterval) < getTimeStep())
  //     {
  //       std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
  //       ", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight <<
  //       ", hEff = " << effectivePistonHeight << ", v = " << pistonVelocity << ", omega = " << pistonAngularVelocity << std::endl <<
  //       "P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << std::endl <<
  //       "torque = " << pistonTorque << " , torque_target = " << targetTorque << std::endl <<
  //       "cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
  //       ", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap <<
  //       ", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
  //       "k1 = " << std::setprecision(6) << std::left << std::setw(10) << speciesBig -> getLoadingStiffness() <<
  //       ", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kc = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl;
  //       std::cout << std::endl;
  //     }
  //   }
  // }

  void actionsAfterSolve() override
  {
    delete [] heightSteps;
    delete [] kcCalibratedArray;
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

  void setWallStiffnessAndRestitutionCoefficients(double kw, double eB, double eS)
  {
    k1Wall = kw;
    bigWallRestitutionCoeff = eB;
    smallWallRestitutionCoeff = eS;
  }

  void setParticleParticleRestitutionCoefficients(double bigBigE, double smallSmallE, double bigSmallE)
  {
    bigBigRestitutionCoeff = bigBigE;
    smallSmallRestitutionCoeff = smallSmallE;
    bigSmallRestitutionCoeff = bigSmallE;
  }

  void setBigParticlePlasticProperties(double k1, double k2max, double kc, double phi)
  {
    k1Big = k1;
    k2MaxBig = k2max;
    kcBig = kc;
    phiBig = phi;
  }

  void setSmallParticlePlasticProperties(double k1, double k2max, double kc, double phi)
  {
    k1Small = k1;
    k2MaxSmall = k2max;
    kcSmall = kc;
    phiSmall = phi;
  }

  void setTorqueCalibrationPoints(int n, double *heightArrayPointer, double *pressureArrayPointer, double *torqueArrayPointer, double pauseDuration)
  {
    nCalibrationDataPoints = n;
    pauseAfterDataPointMet = pauseDuration;
    kcCalibratedArray = new double[n];
    heightSteps = new double[n];
    pressureSteps = new double[n];
    torqueSteps = new double[n];

    for (int i = 0; i < n; i++)
    {
      heightSteps[i] = heightArrayPointer[i];
      pressureSteps[i] = pressureArrayPointer[i];
      torqueSteps[i] = torqueArrayPointer[i];
    }
  }

  void setCalibrationMaximumRelativeDeviations(double maxH, double maxP, double maxT)
  {
    maxHeightRelDeviation = maxH;
    maxPressureRelDeviation = maxP;
    maxTorqueRelDeviation = maxT;
  }

  void setCalibrationRoutineParameters(int nRun, int nMax, double *kcArray, double *devArray, double *dispArray)
  {
    nMaximumRoutineRuns = nMax;
    nCurrentRoutineRun = nRun;

    pointerToInitialKCArray = kcArray;

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

  void setPistonFrictionCoefficient(double pistonMu)
  {
    bigPistonSlidingFrictionCoeff = pistonMu;
  }

  void setPreShearingDuration(double t)
  {
    preShearingDuration = t;
  }




  // FUNCTIONS CALLED IN THE CLASS   ---------------------------------

  void setSpecies()
  {
    speciesHandler.clear();

    // BIG-BIG
    speciesBig = new LinearPlasticViscoelasticFrictionSpecies;
    speciesBig -> setDensity(densityBig);
    speciesBig -> setStiffnessAndRestitutionCoefficient(k1Big, bigBigRestitutionCoeff, massBig);
    speciesBig -> setUnloadingStiffnessMax(k2MaxBig);
    speciesBig -> setCohesionStiffness(kcBig);
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
    speciesSmall -> setCohesionStiffness(kcSmall);
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
    speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kcBig);
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
    speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setCohesionStiffness(kcSmall);
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
    speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kcBig + kcSmall));
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

  void setPistonSpecies()
  {
    // WALL-WALL
    speciesPiston = new LinearPlasticViscoelasticFrictionSpecies;
    speciesPiston -> setDensity(densityWall);
    speciesPiston -> setStiffnessAndRestitutionCoefficient(k1Wall, 1.0, massBig);
    speciesPiston -> setUnloadingStiffnessMax(k1Wall);
    speciesPiston -> setCohesionStiffness(kcBig);
    speciesPiston -> setPenetrationDepthMax(phiBig);

    speciesPiston -> setSlidingFrictionCoefficient(bigPistonSlidingFrictionCoeff);
    speciesPiston -> setSlidingStiffness(speciesPiston -> getLoadingStiffness()*2.0/7.0);
    speciesPiston -> setSlidingDissipation(speciesPiston -> getDissipation()*2.0/7.0);
    speciesPiston -> setRollingFrictionCoefficient(0.0);
    speciesPiston -> setRollingStiffness(speciesPiston -> getLoadingStiffness()*2.0/7.0);
    speciesPiston -> setRollingDissipation(speciesPiston -> getDissipation()*2.0/7.0);
    speciesPiston -> setTorsionFrictionCoefficient(0.0);
    speciesPiston -> setTorsionStiffness(speciesPiston -> getLoadingStiffness()*2.0/7.0);
    speciesPiston -> setTorsionDissipation(speciesPiston -> getDissipation()*2.0/7.0);
    speciesHandler.addObject(speciesPiston);

    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Big + k1Wall), bigWallRestitutionCoeff, massBig);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setUnloadingStiffnessMax(k2MaxBig);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setCohesionStiffness(kcBig);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setPenetrationDepthMax(phiBig);

    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setSlidingFrictionCoefficient(bigPistonSlidingFrictionCoeff);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getDissipation()*2.0/7.0);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setRollingFrictionCoefficient(bigWallRollingFrictionCoeff);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setRollingStiffness(speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getDissipation()*2.0/7.0);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setTorsionFrictionCoefficient(bigWallTorsionFrictionCoeff);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesBig, speciesPiston) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getDissipation()*2.0/7.0);

    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Small + k1Wall), smallWallRestitutionCoeff, massSmall);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setUnloadingStiffnessMax(k2MaxSmall);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setCohesionStiffness(kcSmall);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setPenetrationDepthMax(phiSmall);

    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setSlidingFrictionCoefficient(bigPistonSlidingFrictionCoeff);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> getDissipation()*2.0/7.0);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setRollingFrictionCoefficient(smallWallRollingFrictionCoeff);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setRollingStiffness(speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setRollingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> getDissipation()*2.0/7.0);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setTorsionFrictionCoefficient(smallWallTorsionFrictionCoeff);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesSmall, speciesPiston) -> getDissipation()*2.0/7.0);

    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setStiffnessAndRestitutionCoefficient(k1Wall, 1.0, massSmall);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setUnloadingStiffnessMax(0.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setCohesionStiffness(0.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setPenetrationDepthMax(0.0);

    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setSlidingFrictionCoefficient(0.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesWall, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesWall, speciesPiston) -> getDissipation()*2.0/7.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setRollingFrictionCoefficient(0.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setRollingStiffness(speciesHandler.getMixedObject(speciesWall, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setRollingDissipation(speciesHandler.getMixedObject(speciesWall, speciesPiston) -> getDissipation()*2.0/7.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setTorsionFrictionCoefficient(0.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesWall, speciesPiston) -> getLoadingStiffness()*2.0/7.0);
    speciesHandler.getMixedObject(speciesWall, speciesPiston) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesWall, speciesPiston) -> getDissipation()*2.0/7.0);

    std::cout << "BIG-PISTON stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getDissipation() << "\n";
    std::cout << "BIG-PISTON friction coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getSlidingFrictionCoefficient() << " " << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getRollingFrictionCoefficient() << " " << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getTorsionFrictionCoefficient() << "\n";
    std::cout << "BIG-PISTON max unloading stiffness: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getUnloadingStiffnessMax() << "\n";
    std::cout << "BIG-PISTON cohesive stiffness: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCohesionStiffness() << "\n";
    std::cout << "BIG-PISTON collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getCollisionTime(massBig) << "\n";
    std::cout << "BIG-PISTON collision time / time step: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesPiston) -> getCollisionTime(massBig)/getTimeStep() << "\n\n";
  }

  void computeNumberOfParticles()
  {
    double particleBedVolume = 1.1*constants::pi*pow(casingRadius,2.0)*particleBedHeight;
    nSmall = (particleBedPackingFraction*particleBedVolume/volumeSmall)*(smallToBigMassRatio*densityBig/(densitySmall + smallToBigMassRatio*densityBig));
    nBig = (particleBedPackingFraction*particleBedVolume - nSmall*volumeSmall)/volumeBig;
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

    pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius());
    //pistonHeight = particleBedHeight + radiusSmall;
    pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

    //compressionPiston.setSpecies(speciesWall);
    compressionPiston.setSpecies(speciesPiston);
    compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
    pistonPointer = wallHandler.copyAndAddObject(compressionPiston);
  }

  void makeParticlePiston()
  {
    std::cout << std::endl << "CREATING PARTICLE PISTON" << std::endl << std::endl;
    pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius()) + radiusBig;
    pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

    double relativeLatticeDistance = 1.0;
    double latticeWidth = 2.0*relativeLatticeDistance*radiusBig;
    int nHorizontalWells = (int)(2.0*casingRadius/latticeWidth);
    int nVerticalWells = (int)(4.0*casingRadius/(sqrt(3.0)*latticeWidth));
    if (nHorizontalWells%2) nHorizontalWells++;
    if (nVerticalWells%2) nVerticalWells++;

    Vec3D latticeParticlePosition;
    pistonLatticeParticleCounter = 0;

    SphericalParticle p0;
    p0.setRadius(radiusBig);
    p0.setVelocity(Vec3D(0.0,0.0,pistonMaximumVelocity));
    p0.setSpecies(speciesPiston);

    for (int i=-0.5*nVerticalWells; i<=0.5*nVerticalWells; i++)
    {
      for (int j=-0.5*nHorizontalWells; j<=0.5*nHorizontalWells; j++)
      {
        latticeParticlePosition.X = (j + 0.5*(i%2))*latticeWidth;
        latticeParticlePosition.Y = 0.5*i*sqrt(3.0)*latticeWidth;
        latticeParticlePosition.Z = pistonHeight;

        if (pow(latticeParticlePosition.X,2.0) + pow(latticeParticlePosition.Y,2.0) <= pow(casingRadius,2.0))
        {
          p0.setPosition(latticeParticlePosition);
          p0.fixParticle();
          particleHandler.copyAndAddObject(p0);

          pistonLatticeParticleCounter++;
        }
      }
    }

    std::cout << "PISTON CREATED; COMPOSED OF " << pistonLatticeParticleCounter << " PARTICLES." << std::endl;
    std::cout << "NUMBER OF FREE PARTICLES " << nOfFreeParticles << std::endl;
    std::cout << "REALITY CHECK: " << particleHandler.getNumberOfObjects() << " = " << nOfFreeParticles + pistonLatticeParticleCounter << " ?" << std::endl	<< std::endl;
  }

  void makeSampledParticlePiston()
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

    SphericalParticle p0;
    p0.setVelocity(Vec3D(0.0,0.0,pistonMaximumVelocity));
    p0.setSpecies(speciesPiston);

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
    std::cout << "TARGET HEIGHT " << targetHeight << std::endl;


    // this is the flat wall
    compressionPiston.setSpecies(speciesPiston);
    compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
    pistonPointer = wallHandler.copyAndAddObject(compressionPiston);

    std::cout << "PISTON CREATED; COMPOSED OF " << pistonLatticeParticleCounter << " PARTICLES." << std::endl;
    std::cout << "NUMBER OF FREE PARTICLES " << nOfFreeParticles << std::endl;
    std::cout << "REALITY CHECK: " << particleHandler.getNumberOfObjects() << " = " << nOfFreeParticles + pistonLatticeParticleCounter << " ?" << std::endl	<< std::endl;
  }

  /*
  //void makeSampledParticlePiston()
  //{
  //std::cout << std::endl << "CREATING PARTICLE PISTON FROM SAMPLED PARTICLES" << std::endl << std::endl;
  //pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.1*(particleHandler.getLargestParticle() -> getRadius()) + radiusBig;
  //pistonMaximumVelocity = 0.0001*radiusSmall*pistonVelocityScalingFactor*(1.0 - sizeDispersitySmall)/getTimeStep();

  //double samplingHeight = 0.5*(particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z;
  //double lowerEdgeSampledParticles = getZMax();
  //pistonLatticeParticleCounter = 0;

  //// counts the sampled particles and checks for the lowest edge of the former
  //for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
  //{
  //if (fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight) < particleHandler.getObject(i) -> getRadius())
  //{
  //if ((particleHandler.getObject(i) -> getPosition()).Z - (particleHandler.getObject(i) -> getRadius()) < lowerEdgeSampledParticles) lowerEdgeSampledParticles = (particleHandler.getObject(i) -> getPosition()).Z - (particleHandler.getObject(i) -> getRadius());
  //pistonLatticeParticleCounter++;
  //}
  //}

  //// needs this value tor escale all teh position of the particles accordingly
  //double particlePistonVerticalShift = samplingHeight - lowerEdgeSampledParticles;

  //SphericalParticle p0;
  //p0.setVelocity(Vec3D(0.0,0.0,pistonMaximumVelocity));
  //p0.setSpecies(speciesPiston);

  //Vec3D sampledParticleOriginalPosition;

  //for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
  //{
  //if (fabs((particleHandler.getObject(i) -> getPosition()).Z - samplingHeight) < particleHandler.getObject(i) -> getRadius())
  //{
  //sampledParticleOriginalPosition = particleHandler.getObject(i) -> getPosition();

  //p0.setRadius(particleHandler.getObject(i) -> getRadius());
  //p0.setPosition(Vec3D(sampledParticleOriginalPosition.X, sampledParticleOriginalPosition.Y, sampledParticleOriginalPosition.Z - samplingHeight + pistonHeight + particlePistonVerticalShift));
  //p0.fixParticle();
  //particleHandler.copyAndAddObject(p0);
  //}
  //}

  //// this is the flat wall
  //compressionPiston.setSpecies(speciesPiston);
  //compressionPiston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight + particlePistonVerticalShift));
  //pistonPointer = wallHandler.copyAndAddObject(compressionPiston);

  //std::cout << "PISTON CREATED; COMPOSED OF " << pistonLatticeParticleCounter << " PARTICLES." << std::endl;
  //std::cout << "NUMBER OF FREE PARTICLES " << nOfFreeParticles << std::endl;
  //std::cout << "REALITY CHECK: " << particleHandler.getNumberOfObjects() << " = " << nOfFreeParticles + pistonLatticeParticleCounter << " ?" << std::endl	<< std::endl;
  //}
  */

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

      // piston interactions
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

  void updateTargetTorque()
  {
    targetHeight = heightSteps[calibrationStep];
    targetPressure = pressureSteps[calibrationStep];
    targetTorque = torqueSteps[calibrationStep];

    if (calibrationStep > 0) previousTargetHeight = heightSteps[calibrationStep - 1];
    else previousTargetHeight = pistonHeight;

    std::cout << "PREVIOUS PISTON HEIGHT: " << previousTargetHeight << std::endl;
    std::cout << "NEXT PISTON TARGET HEIGHT: " << targetHeight << std::endl;
    std::cout << "NEXT PISTON TARGET PRESSURE: " << targetPressure << std::endl;
    std::cout << "NEXT PISTON TARGET TORQUE: " << targetTorque << std::endl << std::endl;

    isHeightTargetMet = false;
    isPressureTargetMet = false;
    isTorqueTargetMet = false;
  }

  void resetPistonVelocity()
  {
    if (previousTargetHeight < targetHeight) pistonVelocity = pistonMaximumVelocity;
    else pistonVelocity = -pistonMaximumVelocity;

    std::cout << "RESETTING PISTON VELOCITY. NEW VELOCITY: " << pistonVelocity << std::endl << std::endl;
  }

  void movePiston()
  {
    Vec3D newLatticeParticlePosition;
    pistonHeight += pistonVelocity*getTimeStep();
    effectivePistonHeight += pistonVelocity*getTimeStep();

    for (int i=particleHandler.getNumberOfObjects()-1; i>particleHandler.getNumberOfObjects()-pistonLatticeParticleCounter-1; i--)
    {
      newLatticeParticlePosition = particleHandler.getObject(i) -> getPosition();
      newLatticeParticlePosition.Z += pistonVelocity*getTimeStep();
      particleHandler.getObject(i) -> setPosition(newLatticeParticlePosition);
      particleHandler.getObject(i) -> setVelocity(Vec3D(0.0,0.0,0.0));
    }

    //pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
    //pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));
  }

  void rotatePiston()
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

    pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,pistonAngularVelocity));
    pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
  }

  void wigglePiston()
  {
    double wigglingAmplitude = radiusBig;
    double wigglingFrequency = preShearingDuration/5.0;

    Vec3D newLatticeParticlePosition;
    Vec3D oldLatticeParticlePosition;

    for (int i=particleHandler.getNumberOfObjects()-1; i>particleHandler.getNumberOfObjects()-pistonLatticeParticleCounter-1; i--)
    {
      oldLatticeParticlePosition = particleHandler.getObject(i) -> getPosition();

      newLatticeParticlePosition.X = oldLatticeParticlePosition.X*(1.0 + wigglingAmplitude*cos(2.0*constants::pi*wigglingFrequency*(getTime()-t0))) - oldLatticeParticlePosition.Y*(1.0 + wigglingAmplitude*sin(2.0*constants::pi*wigglingFrequency*(getTime()-t0)));
      newLatticeParticlePosition.Y = oldLatticeParticlePosition.X*(1.0 + wigglingAmplitude*sin(2.0*constants::pi*wigglingFrequency*(getTime()-t0))) + oldLatticeParticlePosition.Y*(1.0 + wigglingAmplitude*cos(2.0*constants::pi*wigglingFrequency*(getTime()-t0)));
      newLatticeParticlePosition.Z = oldLatticeParticlePosition.Z;

      particleHandler.getObject(i) -> setPosition(newLatticeParticlePosition);
    }

    pistonPointer -> setAngularVelocity(Vec3D(0.0,0.0,pistonAngularVelocity));
    pistonPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
  }

  void refreshSpecies()
  {
    // BIG-BIG
    speciesBig -> setCohesionStiffness(kcBig);

    // BIG-WALL
    speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(kcBig);

    // BIG-SMALL
    speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(0.5*(kcBig + kcSmall));
  }

  void tuneKC()
  {
    if (pistonTorque < targetTorque) kcBig *= 1.001;
    else kcBig *= 0.999;

    refreshSpecies();

    if (kcBig < 0.0)
    {
      std::cout << "\n\nFATAL ERROR: KC_BIG < 0.0. QUITTING.";
      setTimeMax(getTime() + 10.0*getTimeStep());
      stage = 10;
    }
  }

  void runPreShearingPhase()
  {
    // the piston needs to reach the initial target value of the height
    if (!isHeightTargetMet && ((pistonVelocity < 0.0 && previousTargetHeight < targetHeight) || (pistonVelocity > 0.0 && previousTargetHeight > targetHeight))) pistonVelocity *= -0.5;

    // upon reaching the target height update the bool and sets the velocity to 0
    if (!isHeightTargetMet && !isPressureTargetMet && !isTorqueTargetMet)
    {
      if (fabs((effectivePistonHeight - targetHeight)/targetHeight) > maxHeightRelDeviation) movePiston(); // IT NOW ACCOUNTS FOR PISTON THICKNESS
      else
      {
        std::cout << "Height target reached.\n";
        std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
        t0 = getTime();

        pistonVelocity = 0.0;
        isHeightTargetMet = true;
      }
    }

    // rotate the piston for a "bit"
    if (isHeightTargetMet && !isPressureTargetMet && !isTorqueTargetMet && getTime() < t0 + preShearingDuration)
    {
      // rotate the piston for a "bit", then move in order to set the pressure
      if (getTime() < t0 + preShearingDuration)
      {
        rotatePiston();
      }
      else
      {
        // go on setting the height to keep the pressure constant and equal to the desired value
        if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -1.0;
        movePiston();
      }
    }
  }

  void runTorqueCalibration()
  {
    // reverse and halve the velocity if the target height was passed
    if (!isHeightTargetMet && ((pistonVelocity < 0.0 && previousTargetHeight < targetHeight) || (pistonVelocity > 0.0 && previousTargetHeight > targetHeight))) pistonVelocity *= -0.5;

    // upon reaching the target height update the bool and sets the velocity to 0
    if (!isHeightTargetMet && !isPressureTargetMet && !isTorqueTargetMet)
    {
      if (fabs((effectivePistonHeight - targetHeight)/targetHeight) > maxHeightRelDeviation) movePiston(); // IT NOW ACCOUNTS FOR PISTON THICKNESS
      else
      {
        std::cout << "Height target reached.\n";
        std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
        t0 = getTime();

        std::cout << "---------------------------------------------\n";
        std::cout << pistonPressure << "\t" << targetPressure << "\t" << maxPressureRelDeviation << std::endl;
        std::cout << "---------------------------------------------\n";

        //pistonVelocity = 0.0;
        //resetPistonVelocity();
        isHeightTargetMet = true;
      }
    }

    // THE SCALING FACTOR MUST BE 1 HERE, INSTEAD OF 0.5!
    if (isHeightTargetMet && fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
    {
      if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -1.0;
      movePiston();
    }

    // when the height is fine and the prescribed time has passed start the piston rotation
    if (isHeightTargetMet && !isPressureTargetMet && !isTorqueTargetMet && getTime() > t0 + pauseAfterDataPointMet)
    {
      rotatePiston();
    }


    // BEGIN ------------------------------ 04-04-18: ignoring the pressureMet, introducing a constant movement to mantain pressure constant and equal to target pressure


    //// after the piston rotated a bit check the pressure and move the piston if needed
    //if (isHeightTargetMet && !isPressureTargetMet && !isTorqueTargetMet && getTime() > t0 + 0.2)
    //{
    //if (fabs((pistonPressure - targetPressure)/targetPressure) > maxPressureRelDeviation)
    //{
    //if ((pistonVelocity < 0.0 && pistonPressure > targetPressure) || (pistonVelocity > 0.0 && pistonPressure < targetPressure)) pistonVelocity *= -0.5;
    //movePiston();
    //}
    //else
    //{
    //std::cout << "Pressure target reached.\n";
    //std::cout << "Now chillin' for " << pauseAfterDataPointMet << " seconds...\n";
    //t0 = getTime();

    //pistonVelocity = 0.0;
    //isPressureTargetMet = true;
    //}
    //}



    //// when the height is fine and the prescribed time has passed start the piston rotation
    //if (isHeightTargetMet && isPressureTargetMet && !isTorqueTargetMet && getTime() > t0 + pauseAfterDataPointMet)
    //{
    //rotatePiston();
    //}


    // END ------------------------------ 04-04-18: ignoring the pressureMet, introducing a constant movement to mantain pressure constant and equal to target pressure










    //// after the piston rotated for a while deal with the torque
    //if (isHeightTargetMet && !isTorqueTargetMet && getTime() > t0 + 2.0*pauseAfterDataPointMet)
    //{
    //if (fabs((pistonTorque - targetTorque)/targetTorque) > maxTorqueRelDeviation)
    //{
    //// the stiffness is changed every n time steps
    ////if (fmod(getTime() - t0, nTimeStepsBetweenStiffnessAdjustments*getTimeStep()) < getTimeStep()) tuneKC();
    //}
    //else
    //{
    //std::cout << "Torque target reached, now chillin' for " << pauseAfterDataPointMet << " seconds\n";
    //t0 = getTime();
    //isTorqueTargetMet = true;
    //}
    //}

    // switch to next calibration point if the test of the previous is passed
    if (isHeightTargetMet && isTorqueTargetMet && getTime() > t0 + pauseAfterDataPointMet)
    {
      // updates the mean value of the stiffnesses at the calibrated data points
      kcCalibratedArray[calibrationStep] = kcBig;
      calibrationStep++;

      if (calibrationStep < nCalibrationDataPoints)
      {
        updateTargetTorque();
        resetPistonVelocity();
      }
      else
      {
        computeVariations();
        updateVariationArrays();

        std::cout << "Compression cycles finished. Quitting...\n";

        makeKdataFile();
        if (nCurrentRoutineRun < nMaximumRoutineRuns) updateInitialKArray();

        setTimeMax(getTime() + 10.0*getTimeStep());
        stage = 10;
      }
    }
  }

  void computeVariations()
  {
    meankcBig = 0.0;
    meanRelativeDeviationkcBig = 0.0;
    meanRelativeDisplacementekcBig = 0.0;

    for (int i = 0; i < nCalibrationDataPoints; i++) meankcBig += kcCalibratedArray[i];
    meankcBig /= nCalibrationDataPoints;

    for (int i = 0; i < nCalibrationDataPoints; i++) meanRelativeDeviationkcBig += fabs(kcCalibratedArray[i] - meankcBig);
    meanRelativeDeviationkcBig /= nCalibrationDataPoints*meankcBig;

    meanRelativeDisplacementekcBig = fabs(meankcBig - pointerToInitialKCArray[nCurrentRoutineRun])/pointerToInitialKCArray[nCurrentRoutineRun];
  }

  void updateVariationArrays()
  {
    pointerToDeviationArray[nCurrentRoutineRun] = meanRelativeDeviationkcBig;
    pointerToDisplacementArray[nCurrentRoutineRun] = meanRelativeDisplacementekcBig;
  }

  void updateInitialKArray()
  {
    pointerToInitialKCArray[nCurrentRoutineRun + 1] = 0.5*(pointerToInitialKCArray[nCurrentRoutineRun] + meankcBig);
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
    kdataFile << "|kc| = " << meankcBig << std::endl;
    kdataFile << "MD(kc)/|kc| = " << pointerToDeviationArray[nCurrentRoutineRun] << std::endl;
    kdataFile << "||kc| - kC_init|/kC_init = " << pointerToDisplacementArray[nCurrentRoutineRun] << std::endl;
    kdataFile << "height \t torque  \t kc" << std::endl;
    for (int i = 0; i < nCalibrationDataPoints; i++)
    {
      kdataFile << heightSteps[i] << "\t" << torqueSteps[i] << "\t" << kcCalibratedArray[i] << std::endl;
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
    cdatFile << "time \t Eel \t Ekin/Eel \t h \t v \t omega \t P \t Pbase \t Pcyl \t torque \t torque_target \t |(T - Ttarg)/Ttarg| \t dTotMean \t dTotMax \t dPistonMean \t dPistonMax \t dBaseMean \t dBaseMax \t k1 \t k2 \t kc \t phi" << std::endl;
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
    kcBig << "   " <<
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
      ", k2 = " << speciesBig -> getUnloadingStiffnessMax() << ", kc = " << speciesBig -> getCohesionStiffness() << ", phi = " << speciesBig -> getPenetrationDepthMax() << std::endl <<
      "HEIGHT MET: " << isHeightTargetMet << ", PRESSURE MET: " << isPressureTargetMet << ", TORQUE MET: " << isTorqueTargetMet << std::endl << std::endl;
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
  double nBig, nSmall;
  int nOfFreeParticles;

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

  // friction
  double bigWallSlidingFrictionCoeff, smallWallSlidingFrictionCoeff;
  double bigWallRollingFrictionCoeff, smallWallRollingFrictionCoeff;
  double bigWallTorsionFrictionCoeff, smallWallTorsionFrictionCoeff;
  double bigBigSlidingFrictionCoeff, smallSmallSlidingFrictionCoeff, bigSmallSlidingFrictionCoeff;
  double bigBigRollingFrictionCoeff, smallSmallRollingFrictionCoeff, bigSmallRollingFrictionCoeff;
  double bigBigTorsionFrictionCoeff, smallSmallTorsionFrictionCoeff, bigSmallTorsionFrictionCoeff;
  double bigPistonSlidingFrictionCoeff;

  // collision
  double k1Wall, k1Big, k1Small;
  double k2MaxBig, k2MaxSmall;
  double kcBig, kcSmall;
  double phiBig, phiSmall;
  double phiMin;
  double bigWallRestitutionCoeff, smallWallRestitutionCoeff;
  double bigBigRestitutionCoeff, smallSmallRestitutionCoeff, bigSmallRestitutionCoeff;

  // species
  LinearPlasticViscoelasticFrictionSpecies *speciesBig, *speciesSmall, *speciesWall, *speciesPiston;

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
  double *kcCalibratedArray;
  double *heightSteps;
  double *pressureSteps;
  double *torqueSteps;
  bool isPreShearingDone;
  bool isHeightTargetMet;
  bool isPressureTargetMet;
  bool isTorqueTargetMet;
  int calibrationStep;
  double targetHeight;
  double targetPressure;
  double targetTorque;
  double previousTargetHeight;
  double maxHeightRelDeviation;
  double maxTorqueRelDeviation;
  double maxPressureRelDeviation;
  double nTimeStepsBetweenStiffnessAdjustments;
  double *pointerToInitialKCArray;
  double *pointerToDeviationArray;
  double *pointerToDisplacementArray;
  double meankcBig;
  double meanRelativeDeviationkcBig;
  double meanRelativeDisplacementekcBig;
  double preShearingDuration;

  // global
  int stage;
  double t0;
  int nMaximumRoutineRuns;
  int nCurrentRoutineRun;
};


int main(int argc, char *argv[])
{
  // GLOBAL FLAGS

  // TIME STEP
  double timeStep = 1.0e-5;

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
  double pistonAngularVelocity = 0.1*constants::pi/180.0*100.0;

  // CALIBRATION DATA POINTS
  double pauseDuration = 0.002;
  const int nDataPoints = 2;
  double heightArray[nDataPoints] = {0.0185328, 0.0184762};
  double pressureArray[nDataPoints] = {5000.0, 6000.0};
  double torqueArray[nDataPoints] = {0.00345, 0.00526};

  // INTERACTION PARAMETERS
  double k1Big = 600.0;
  double k2MaxBig = 5000.0;
  double kcBig = 100.0;
  double phiBig = 0.05;

  double k1Small = 400.0;
  double k2MaxSmall = 5000.0;
  double kcSmall = 100.0;
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
  double muSBigWall = 0.30;
  double muRBigWall = 0.05;

  double muSPiston = muSBigBig;

  // CALIBRATION PARAMETERS
  double maximumStiffnessRelativeDeviation = 0.05;
  double maximumHeightRelativeDeviation = 1.0e-4;
  double maximumPressureRelativeDeviation = 1.0e-4;
  double maximumTorqueRelativeDeviation = 1.0e-4;
  double timeStepsBetweenStiffnessTunings = 10.0;

  double preShearingDuration = 0.5;

  int runNumber = 0;
  const int nMaximumRoutineCycles = 10;
  double initialKCArray[nMaximumRoutineCycles];
  double deviationArray[nMaximumRoutineCycles];
  double displacementArray[nMaximumRoutineCycles];

  for (int i = 0; i < nMaximumRoutineCycles; i++)
  {
    initialKCArray[i] = 0.0;
    deviationArray[i] = 1.0;
    displacementArray[i] = 1.0;
  }

  initialKCArray[0] = kcBig;

  // STRING VARIABLES
  std::ostringstream name;
  name.str("");
  name.clear();

  if (argc > 1) // THE CALIBRATION ROUTINE
  {
    do
    {
      // --- TESTING FUNCTIONS ---
      //std::cout << std::endl << std::endl;
      //std::cout << "BEFORE EXECUTION" << std::endl;
      //std::cout << "KC_init \t devKC \t dispKC" << std::endl;
      //for (int i = 0; i < nMaximumRoutineCycles; i++) std::cout << initialKCArray[i] << "\t" << deviationArray[i] << "\t" << displacementArray[i] << std::endl;
      //std::cout << std::endl << std::endl;
      // --- --- --- --- --- ---

      TorqueTest_calibrationRoutine problem;

      problem.readArguments(argc, argv);

      // INITIALIZATION
      problem.setTimeStep(timeStep);
      problem.setTimeMax(1.0);
      problem.setGravity(Vec3D(0.00,0.00,-9.81));
      problem.setSystemDimensions(3);
      //problem.setVerbose(true);
      problem.setCdatOutputTimeInterval(0.001);
      problem.setSaveCount(0.01/problem.getTimeStep());

      problem.setCasingProperties(casingRadius, casingHeight, densityWalls);
      problem.setPowderBedProperties(particleBedHeight, particleBedPackingFraction);
      problem.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);

      problem.setParticleWallSlidingFrictionCoefficients(muSBigWall, muSBigWall);
      problem.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
      problem.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);
      problem.setPistonFrictionCoefficient(muSPiston);

      problem.setParticleParticleSlidingFrictionCoefficients(muSBigBig, muSBigBig, muSBigBig);
      problem.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
      problem.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);

      problem.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
      problem.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
      problem.setBigParticlePlasticProperties(k1Big, k2MaxBig, initialKCArray[runNumber], phiBig);
      problem.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kcSmall, phiSmall);

      problem.setTorqueCalibrationPoints(nDataPoints, heightArray, pressureArray, torqueArray, pauseDuration);
      problem.setCalibrationMaximumRelativeDeviations(maximumHeightRelativeDeviation, maximumPressureRelativeDeviation, maximumTorqueRelativeDeviation);
      problem.setCalibrationRoutineParameters(runNumber, nMaximumRoutineCycles, initialKCArray, deviationArray, displacementArray);
      problem.setNTimeStepsBetweenStiffnesAdjustments(timeStepsBetweenStiffnessTunings);
      problem.setPreShearingDuration(preShearingDuration);

      problem.setPistonRotation(pistonAngularVelocity);
      problem.setPistonVelocityScalingFactor(pistonVelocityScalingFactor);

      // NAME SETTING
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      //name << "TorqueTest_NEWROUTINE_CARRIER_Nrun_" << runNumber << "_rRatio_" << std::fixed << std::setprecision(0) <<  bigToCasingSizeRatio <<
      //"_k1_" << k1Big << "_k2max_" << k2MaxBig << "_kc_" << initialKCArray[runNumber] << "_phi_" << std::fixed << std::setprecision(4) << phiBig <<
      //"_muSBB_" << std::fixed << std::setprecision(2) << muSBigBig << "_muRBB_" << muRBigBig << "_muSBW_" << muSBigWall << "_muRBW_" << muRBigWall << "_muSpiston_" << muSPiston <<
      //"___TESTING_NOTUNING_" << std::fixed << std::setprecision(6) << pistonAngularVelocity;
      name << "TorqueTest_TESTING_kC_" << std::fixed << std::setprecision(2) << kcBig << "_muSPiston_" << muSPiston;
      problem.setName(name.str());
      //problem.setName("TorqueTest_TESTING");

      //std::cout << "This file has been restarted.\nThe new name of the associated output of this simulation will be:\n" << name.str() << std::endl << std::endl;

      problem.solve();
      runNumber++;

      // --- TESTING FUNCTIONS ---
      //std::cout << std::endl << std::endl;
      //std::cout << "AFTER EXECUTION" << std::endl;
      //std::cout << "KC_init \t devKC \t dispKC" << std::endl;
      //for (int i = 0; i < nMaximumRoutineCycles; i++) std::cout << initialKCArray[i] << "\t" << deviationArray[i] << "\t" << displacementArray[i] << std::endl;
      //std::cout << std::endl << std::endl;
      // --- --- --- --- --- ---
    }
    while (runNumber < nMaximumRoutineCycles && !(displacementArray[runNumber - 1] < maximumTorqueRelativeDeviation && deviationArray[runNumber - 1] < maximumTorqueRelativeDeviation));
  }
  else
  {
    std::cout << "THIS PROGRAM IS SUPPOSED TO BE A RESTART-ONLY." << std::endl;
    std::cout << "NO RESTART FILE WAS SPECIFIED." << std::endl;
    std::cout << "QUITTING." << std::endl;
  }
  return 0;
}
