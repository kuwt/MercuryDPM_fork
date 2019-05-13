#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

class CalibrationRoutine_SETTLE : public Mercury3D
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

   void actionsOnRestart() override
   {

   }

   void actionsAfterTimeStep() override
   {
      // particle insertion loop
      if (stage == 1) insertParticles();

      // removal of the particle above the casing and final settling
      if (stage == 2 && getKineticEnergy()/getElasticEnergy() < maximumEnergyRatio)
      {
         std::cout << "Particle settle terminated." << std::endl << "Removing the exceeding particles..." << std::endl;
         removeParticles();
         std::cout << "DONE." << std::endl << "Settling the particles..." << std::endl;
         t0 = getTime();

         stage++;
      }

      if (stage == 3 && getKineticEnergy()/getElasticEnergy() < maximumEnergyRatio && getTime() - t0 > 10.0*getTimeStep())
      {
         std::cout << "Particle settle terminated." << std::endl << "Quitting..." << std::endl;
         setTimeMax(getTime() + getTimeStep());

         stage++;
      }
   }

   void actionAfterSolve()
   {

   }


public:
   // FUNCTIONS CALLED IN MAIN  ---------------------------------------
   void setCasingProperties(double radius, double height, double pf, double density)
   {
      casingRadius = radius;
      particleBedHeight = height;
      particleBedPackingFraction = pf;
      densityWall = density;

      setXMin(-1.1*casingRadius);
      setYMin(-1.1*casingRadius);
      setZMin(-0.1*casingRadius);

      setXMax(1.1*casingRadius);
      setYMax(1.1*casingRadius);
      setZMax(2.5*particleBedHeight);
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

   void setMaximumEnergyRatio(double maxRatio)
   {
      maximumEnergyRatio = maxRatio;
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
      double particleBedVolume = 1.1*particleBedHeight*constants::pi*pow(casingRadius,2.0);
      nSmall = (particleBedPackingFraction*particleBedVolume/volumeSmall)*(smallToBigMassRatio*densityBig/(densitySmall + smallToBigMassRatio*densityBig));
      nBig = (particleBedPackingFraction*particleBedVolume - nSmall*volumeSmall)/volumeBig;
   }

   // makes the geometric components
   void makeGeometry()
   {
      wallHandler.clear();

      // the basis of the compaction chamber
      base.setSpecies(speciesWall);
      base.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
      wallHandler.copyAndAddObject(base);

      // the roof of the compaction chamber
      roof.setSpecies(speciesWall);
      roof.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
      wallHandler.copyAndAddObject(roof);

      // the external casing
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
         particlePosition.Z = random.getRandomNumber(1.5*particleBedHeight + 1.01*(p0.getRadius()), getZMax() - 1.01*(p0.getRadius()));

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

      if (waitForParticleSettling && getTime() > t0 + sqrt(2.0*(getZMax() - 1.5*particleBedHeight)/9.81)) waitForParticleSettling = false;
      if (!waitForParticleSettling && nBigInserted >= nBig && nSmallInserted >= nSmall && getKineticEnergy()/getElasticEnergy() < maximumEnergyRatio)
      {
         std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
         stage++;
      }
   }

   // removes the particles laying above the casing height
   void removeParticles()
   {
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         if (particleHandler.getObject(i) -> getPosition().Z + particleHandler.getObject(i) -> getRadius() > particleBedHeight) particleHandler.removeObject(i);
      }
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

   // powder bed
   double particleBedHeight;
   double particleBedPackingFraction;

   // geometry
   double casingRadius;
   double densityWall;
   InfiniteWall base, roof;
   AxisymmetricIntersectionOfWalls casing;

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

   // global
   int stage;
   bool waitForParticleSettling;
   double t0;
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

int main(int argc, char *argv[])
{
   // immediately aborts if too many arguments were given in command line
   if (argc > 1)
   {
      std::cout << "Too many input arguments given." << std::endl;
      std::cout << "QUITTING..." << std::endl;

      exit(1);
   }

   // DATA FILENAME  SPECIFICATION
    std::string material = "CAPSULAC60";
    std::string dataPointsFileNameCOMPRESSION = "Calibration_dataPoints/Calibration_dataPoints_COMPRESSION_CAPSULAC60.dat";

   // DATA POINT SELECTION FOR ANALYSIS
   int nTotalDataPoints = 15;
   double heightArray[nTotalDataPoints];
   double pressureArray[nTotalDataPoints];

   // TIME STEP
   double timeStep = 1.0e-6;

   // GEOMETRY PARAMETERS
   double particleBedHeight = 0.019;
   double particleBedPackingFraction = 0.60;
   double casingRadius = 0.0125;
   double densityWalls = 3000.0;

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
   double muSBigBig = 0.28;
   double muRBigBig = 0.10;
   double muSBigWall = muSBigBig;
   double muRBigWall = muRBigBig;

   // CALIBRATION PARAMETERS
   double maximumEnergyRatio = 1.0e-4;

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

   // THE SIMULATION SETUP
   CalibrationRoutine_SETTLE crSettle;

   // NAME SETTING
   name.str("");
   name.clear();
   std::cout.unsetf(std::ios::floatfield);
   name << "CalibrationTest_SETTLE_" << material <<
   "_k1_" << std::fixed << std::setprecision(0) << k1Big << "_k2max_" << k2MaxBig << "_kc_" << kCBig << "_phi_" << std::fixed << std::setprecision(4) << phiBig <<
   "_muS_" << std::fixed << std::setprecision(2) << muSBigBig << "_muR_" << muRBigBig <<
   "_cConstant_" << std::fixed << std::setprecision(4) << cConstant[0] << "_DhatMax_" << deltaHatMax[0];
   crSettle.setName(name.str());

   // INITIALIZATION
   crSettle.setTimeStep(timeStep);
   crSettle.setTimeMax(100.0);
   crSettle.setGravity(Vec3D(0.00,0.00,-9.81));
   crSettle.setSystemDimensions(3);
   crSettle.setSaveCount(0.01/crSettle.getTimeStep());

   crSettle.setCasingProperties(casingRadius, particleBedHeight, particleBedPackingFraction, densityWalls);
   crSettle.setParticleProperties(casingRadius/bigToCasingSizeRatio, casingRadius/bigToCasingSizeRatio/bigToSmallSizeRatio, sizeDispersionBigParticles, sizeDispersionSmallParticles, densityBigParticles, densitySmallParticles, smallToBigMassRatio);

   crSettle.setParticleWallSlidingFrictionCoefficients(muSBigWall, muSBigWall);
   crSettle.setParticleWallRollingFrictionCoefficients(muRBigWall, muRBigWall);
   crSettle.setParticleWallTorsionFrictionCoefficients(0.0, 0.0);

   crSettle.setParticleParticleSlidingFrictionCoefficients(muSBigBig, muSBigBig, muSBigBig);
   crSettle.setParticleParticleRollingFrictionCoefficients(muRBigBig, muRBigBig, muRBigBig);
   crSettle.setParticleParticleTorsionFrictionCoefficients(0.0, 0.0, 0.0);

   crSettle.setWallStiffnessAndRestitutionCoefficients(k1Wall, eBigWall, eSmallWall);
   crSettle.setParticleParticleRestitutionCoefficients(eBigBig, eSmallSmall, eBigSmall);
   crSettle.setBigParticlePlasticProperties(k1Big, k2MaxBig, kCBig, phiBig);
   crSettle.setSmallParticlePlasticProperties(k1Small, k2MaxSmall, kCSmall, phiSmall);

   crSettle.setMaximumEnergyRatio(maximumEnergyRatio);

   crSettle.solve();

   return 0;
}
