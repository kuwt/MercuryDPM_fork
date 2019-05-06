#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

class CalibrationRoutine_ParticleSettling : public Mercury3D
{
private:

   void setupInitialConditions() override
   {
      stage = 1;
      setsInserted = 0;
      thermalizationInstancesDone = 0;
      thermalizationCyclesCounter = 0;
      maximumPackingFraction = 0.0;

      // sets particle masses
      std::cout << "Setting particles masses...\n";
      setParticleMassAndVolume();
      std::cout << "DONE\n";

      // sets the species
      std::cout << "Setting up the species...\n";
      setSpecies();
      std::cout << "DONE\n";

      // computes the number of particle sets needed
      std::cout << "Computing the number of particle sets and the total number of particles needed...\n";
      // computeNumberOfSetsAndParticles();
      computeNumberOfParticles();
      std::cout << "DONE\n";

      // computes the total volume of particles
      std::cout << "Computing the total volume of particles...\n";
      computeParticleTotalVolume();
      std::cout << "DONE\n";

      // computes the particle level and the system max height
      std::cout << "Computing the height of the settled particle bed and the total height of the system...\n";
      computeHeights();
      std::cout << "DONE\n";

      // setting the simulation domain
      std::cout << "Setting the simulation domain...\n";
      setBoundaries();
      std::cout << "DONE\n";

      // creating the geometrical components of the casing
      std::cout << "Creating the geometry...\n";
      makeGeometry();
      std::cout << "DONE\n";

      // inserting the biggest particle to please mercury 3D grid
      std::cout << "Inserting the biggest particle...\n";
      insertBiggestParticle();
      std::cout << "DONE\n";

      std::cout << "Inserting the rest of the particles...\n";
      stage++;
   }

   void actionsOnRestart() override
   {

   }

   void actionsAfterTimeStep() override
   {
      // particle insertion loop
      if (stage == 2)
      {
         insertParticles();
      }

      // tapping loop
      if (stage == 0)
      {
         makeDataAnalysis();

         // updates the .cdat output file
         if (fmod(getTime(),cdatOutputTimeInterval) < getTimeStep()) writeDataToOutptFile();

         // if the tapping instances done, the energy ratio and the time passed since last tapping are all below the respective thresholds this kicks is
         if (thermalizationInstancesDone < numberOfThermalizationInstances && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > timePlaceholder + thermalizationInstancesDone*thermalizationTimeInterval)
         {
            // checks if the packing fraction changed by an amount above the threshold. If so sets the thermalization cycles counter to 0, if not increases it by one
            if (fabs(packingFraction - maximumPackingFraction) < packingFractionDifferenceTolerance) thermalizationCyclesCounter++;
            else thermalizationCyclesCounter = 0;

            // if the thermalization cycles counter is below the threshold thermalizes the system
            if (thermalizationCyclesCounter < packingFractionCyclesThreshold)
            {
               thermalizeParticles();

               // if the packing fraction is greater than the packing fraction at the last tapping instance it updates the latter with the value of the former
               if (maximumPackingFraction < packingFraction) maximumPackingFraction = packingFraction;
               timePlaceholder = getTime();
               thermalizationInstancesDone++;
            }
         }

         // if either the max number of thermalization instances is done or the packing fraction stabilised passes to the next stage
         if (thermalizationInstancesDone == numberOfThermalizationInstances || thermalizationCyclesCounter == packingFractionCyclesThreshold)
         {
            timePlaceholder = getTime();
            stage++;
         }
      }

      // removal of the particle above the casing and final settling
      if (stage == 3 && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
      {
         std::cout << "Particle thermalization terminated.\nRemoving the exceeding particles...\n";
         removeParticles();
         std::cout << "DONE.\nSettling the particles...\n";
         timePlaceholder = getTime();
         stage++;
      }
   }

   void actionAfterSolve()
   {
      // closes the stream to the output file
      outputFile.close();
   }

public:

   // ----- FUNCTIONS CALLED IN THE MAIN -----
   // sets verbosity on or off
   void setVerbose(bool v)
   {
      verbose = v;
   }

   // set the particle density (OVERLOADED)
   void setParticleDensity(double rhoB)
   {
      particleDensityBig = rhoB;
      particleDensitySmall = rhoB;
   }
   void setParticleDensity(double rhoB, double rhoS)
   {
      particleDensityBig = rhoB;
      particleDensitySmall = rhoS;
   }

   // set the mean and the dispersity of the particles (OVERLOADED)
   void setParticleRadiusAndDispersity(double mu, double sigma)
   {
      meanRadiusBig = mu;
      meanRadiusSmall = mu;
      dispersityBig = sigma;
      dispersitySmall = sigma;

      biModal = false;
   }
   void setParticleRadiusAndDispersity(double muBig, double sigmaBig, double muSmall, double sigmaSmall)
   {
      meanRadiusBig = muBig;
      meanRadiusSmall = muSmall;
      dispersityBig = sigmaBig;
      dispersitySmall = sigmaSmall;

      biModal = true;
   }

   // set the small-to-big total mass ratio
   void setTotalMassRatio(double ratio)
   {
      smallToBigMassRatio = ratio;
   }

   // set wall stiffness
   void setWallStiffness(double kW)
   {
      wallStiffness = kW;
   }

   // set particle stiffness (OVERLOADED)
   void setParticleStiffness(double k)
   {
      particleStiffnessBig = k;
      particleStiffnessSmall = k;
   }
   void setParticleStiffness(double kBig, double kSmall)
   {
      particleStiffnessBig = kBig;
      particleStiffnessSmall = kSmall;
   }

   // set particle-wall sliding friction coefficients (OVERLOADED)
   void setParticleWallSlidingFrictionCoeff(double bigWallMu)
   {
      bigWallSlidingFrictionCoeff = bigWallMu;
      smallWallSlidingFrictionCoeff = bigWallMu;
   }
   void setParticleWallSlidingFrictionCoeff(double bigWallMu, double smallWallMu)
   {
      bigWallSlidingFrictionCoeff = bigWallMu;
      smallWallSlidingFrictionCoeff = smallWallMu;
   }

   // set particle-wall rolling friction coefficients (OVERLOADED)
   void setParticleWallRollingFrictionCoeff(double bigWallMu)
   {
      bigWallRollingFrictionCoeff = bigWallMu;
      smallWallRollingFrictionCoeff = bigWallMu;
   }
   void setParticleWallRollingFrictionCoeff(double bigWallMu, double smallWallMu)
   {
      bigWallRollingFrictionCoeff = bigWallMu;
      smallWallRollingFrictionCoeff = smallWallMu;
   }

   // set particle-wall torsion friction coefficients (OVERLOADED)
   void setParticleWallTorsionFrictionCoeff(double bigWallMu)
   {
      bigWallTorsionFrictionCoeff = bigWallMu;
      smallWallTorsionFrictionCoeff = bigWallMu;
   }
   void setParticleWallTorsionFrictionCoeff(double bigWallMu, double smallWallMu)
   {
      bigWallTorsionFrictionCoeff = bigWallMu;
      smallWallTorsionFrictionCoeff = smallWallMu;
   }

   // set particle-particle sliding friction coefficients (OVERLOADED)
   void setParticleParticleSlidingFrictionCoeff(double bigBigMu)
   {
      bigBigSlidingFrictionCoeff = bigBigMu;
      smallSmallSlidingFrictionCoeff = bigBigMu;
      bigSmallSlidingFrictionCoeff = bigBigMu;
   }
   void setParticleParticleSlidingFrictionCoeff(double bigBigMu, double smallSmallMu, double bigSmallMu)
   {
      bigBigSlidingFrictionCoeff = bigBigMu;
      smallSmallSlidingFrictionCoeff = smallSmallMu;
      bigSmallSlidingFrictionCoeff = bigSmallMu;
   }

   // set particle-particle rolling friction coefficients (OVERLOADED)
   void setParticleParticleRollingFrictionCoeff(double bigBigMu)
   {
      bigBigRollingFrictionCoeff = bigBigMu;
      smallSmallRollingFrictionCoeff = bigBigMu;
      bigSmallRollingFrictionCoeff = bigBigMu;
   }
   void setParticleParticleRollingFrictionCoeff(double bigBigMu, double smallSmallMu, double bigSmallMu)
   {
      bigBigRollingFrictionCoeff = bigBigMu;
      smallSmallRollingFrictionCoeff = smallSmallMu;
      bigSmallRollingFrictionCoeff = bigSmallMu;
   }

   // set particle-particle torsion friction coefficients (OVERLOADED)
   void setParticleParticleTorsionFrictionCoeff(double bigBigMu)
   {
      bigBigTorsionFrictionCoeff = bigBigMu;
      smallSmallTorsionFrictionCoeff = bigBigMu;
      bigSmallTorsionFrictionCoeff = bigBigMu;
   }
   void setParticleParticleTorsionFrictionCoeff(double bigBigMu, double smallSmallMu, double bigSmallMu)
   {
      bigBigTorsionFrictionCoeff = bigBigMu;
      smallSmallTorsionFrictionCoeff = smallSmallMu;
      bigSmallTorsionFrictionCoeff = bigSmallMu;
   }

   // set particle-wall restitution coefficients (OVERLOADED)
   void setParticleWallRestitutionCoeff(double bigWallE)
   {
      bigWallRestitutionCoeff = bigWallE;
      smallWallRestitutionCoeff = bigWallE;
   }
   void setParticleWallRestitutionCoeff(double bigWallE, double smallWallE)
   {
      bigWallRestitutionCoeff = bigWallE;
      smallWallRestitutionCoeff = smallWallE;
   }

   // set particle-particle restitution coefficients (OVERLOADED)
   void setParticleParticleRestitutionCoeff(double bigBigE)
   {
      bigBigRestitutionCoeff = bigBigE;
      smallSmallRestitutionCoeff = bigBigE;
      bigSmallRestitutionCoeff = bigBigE;
   }
   void setParticleParticleRestitutionCoeff(double bigBigE, double smallSmallE, double bigSmallE)
   {
      bigBigRestitutionCoeff = bigBigE;
      smallSmallRestitutionCoeff = smallSmallE;
      bigSmallRestitutionCoeff = bigSmallE;
   }

   // set particle plastic properties
   void setParticlePlasticProperties(double k1, double k2max, double kC, double phiC)
   {
      particleStiffnessBig = k1;
      particleStiffnessSmall = k1;

      particleMaxUnloadingStiffness = k2max;
      particleCohesiveStiffness = kC;
      particlePlasticityDepth = phiC;
   }

   //    // set particle adhesion properties
   //    void setParticlesAdhesionProperties(double kA, double fAmax)
   //    {
   //        particleAdhesionStiffness = kA;
   //        particleMaxAdhesiveForce = fAmax;
   //    }

   // set compression cylinder dimensions
   void setCasingProperties(double radius, double height, double density)
   {
      casingRadius = radius;
      casingHeight = height;
      densityWall = density;
      casingVolume = constants::pi*pow(casingRadius,2.0)*casingHeight;
   }

   // set bulk density of the mixture prior to compression
   void setbulkPackingFractionPreCompression(double bulkPF)
   {
      bulkPackingFractionPreCompression = bulkPF;
   }

   // sets teh thermalization parameters
   void setThermalizationParameters(int nInstances, double scalingFactor, double dt)
   {
      numberOfThermalizationInstances = nInstances;
      thermalizationEnergyScalingFactor = scalingFactor;
      thermalizationTimeInterval = dt;
   }

   // sets the energy ratio tolerance
   void setEnergyRatioTolerance(double eRatio)
   {
      energyRatioTolerance = eRatio;
   }

   // sets the tolerance for the packing fraction difference
   void setPackingFractionDifferenceTolerance(double pfRatio, int nCycles)
   {
      packingFractionDifferenceTolerance = pfRatio;
      packingFractionCyclesThreshold = nCycles;
   }

   // sets the time between printouts of data to the cdat file
   void setCdatOutputTimeInterval(double dt)
   {
      cdatOutputTimeInterval = dt;
   }


   // ----- FUNCTIONS CALLED IN THE CLASS -----
   // set the particle masses
   void setParticleMassAndVolume()
   {
      particleVolumeBig = 4.0*constants::pi*pow(meanRadiusBig,3.)/3.;
      particleVolumeSmall = 4.0*constants::pi*pow(meanRadiusSmall,3.)/3.;

      particleMassBig = particleDensityBig*particleVolumeBig;
      particleMassSmall = particleDensitySmall*particleVolumeSmall;

      if (verbose)
      {
         if (biModal)
         {
            std::cout << "\tParticle mass BIG: " << particleMassBig << "\n";
            std::cout << "\tParticle mass SMALL: " << particleMassSmall << "\n";
         }
         else std::cout << "\tParticle mass: " << particleMassBig << "\n";
      }
   }

   // set the particle species
   void setSpecies()
   {
      speciesHandler.clear();

      // BIG-BIG
      //        speciesBig = new LinearViscoelasticFrictionSpecies;
      speciesBig = new LinearPlasticViscoelasticFrictionSpecies;
      speciesBig -> setDensity(particleDensityBig);
      speciesBig -> setStiffnessAndRestitutionCoefficient(particleStiffnessBig, bigBigRestitutionCoeff, particleMassBig);
      // plastic-adhesive part
      speciesBig -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
      speciesBig -> setCohesionStiffness(particleCohesiveStiffness);
      speciesBig -> setPenetrationDepthMax(particlePlasticityDepth);
      //        speciesBig -> setAdhesionStiffness(particleAdhesionStiffness);
      //        speciesBig -> setAdhesionForceMax(particleMaxAdhesiveForce);

      speciesBig -> setSlidingFrictionCoefficient(bigBigSlidingFrictionCoeff);
      speciesBig -> setSlidingStiffness(particleStiffnessBig*2.0/7.0);
      speciesBig -> setSlidingDissipation(speciesBig -> getDissipation()*2.0/7.0);

      speciesBig -> setRollingFrictionCoefficient(bigBigRollingFrictionCoeff);
      speciesBig -> setRollingStiffness(particleStiffnessBig*2.0/7.0);
      speciesBig -> setRollingDissipation(speciesBig -> getDissipation()*2.0/7.0);

      speciesBig -> setTorsionFrictionCoefficient(bigBigTorsionFrictionCoeff);
      speciesBig -> setTorsionStiffness(particleStiffnessBig*2.0/7.0);
      speciesBig -> setTorsionDissipation(speciesBig -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesBig);

      // SMALL-SMALL
      //        speciesSmall = new LinearViscoelasticFrictionSpecies;
      speciesSmall = new LinearPlasticViscoelasticFrictionSpecies;
      speciesSmall -> setDensity(particleDensitySmall);
      speciesSmall -> setStiffnessAndRestitutionCoefficient(particleStiffnessSmall, smallSmallRestitutionCoeff, particleMassSmall);
      // plastic-adhesive part
      speciesSmall -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
      speciesSmall -> setCohesionStiffness(particleCohesiveStiffness);
      speciesSmall -> setPenetrationDepthMax(particlePlasticityDepth);
      //        speciesSmall -> setAdhesionStiffness(particleAdhesionStiffness);
      //        speciesSmall -> setAdhesionForceMax(particleMaxAdhesiveForce);

      speciesSmall -> setSlidingFrictionCoefficient(smallSmallSlidingFrictionCoeff);
      speciesSmall -> setSlidingStiffness(particleStiffnessSmall*2.0/7.0);
      speciesSmall -> setSlidingDissipation(speciesSmall -> getDissipation()*2.0/7.0);

      speciesSmall -> setRollingFrictionCoefficient(smallSmallRollingFrictionCoeff);
      speciesSmall -> setRollingStiffness(particleStiffnessSmall*2.0/7.0);
      speciesSmall -> setRollingDissipation(speciesSmall -> getDissipation()*2.0/7.0);

      speciesSmall -> setTorsionFrictionCoefficient(smallSmallTorsionFrictionCoeff);
      speciesSmall -> setTorsionStiffness(particleStiffnessSmall*2.0/7.0);
      speciesSmall -> setTorsionDissipation(speciesSmall -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesSmall);

      // WALL-WALL
      speciesWall = new LinearPlasticViscoelasticFrictionSpecies;
      speciesWall -> setDensity(densityWall);
      speciesWall -> setStiffnessAndRestitutionCoefficient(wallStiffness, 1.0, particleMassSmall);
      speciesWall -> setUnloadingStiffnessMax(wallStiffness);
      speciesWall -> setCohesionStiffness(0.0);
      speciesWall -> setPenetrationDepthMax(0.0);
      //        speciesWall -> setAdhesionStiffness(0.0);
      //        speciesWall -> setAdhesionForceMax(0.0);

      speciesWall -> setSlidingFrictionCoefficient(0.0);
      speciesWall -> setSlidingStiffness(wallStiffness*2.0/7.0);
      speciesWall -> setSlidingDissipation(speciesWall -> getDissipation()*2.0/7.0);

      speciesWall -> setRollingFrictionCoefficient(0.0);
      speciesWall -> setRollingStiffness(wallStiffness*2.0/7.0);
      speciesWall -> setRollingDissipation(speciesWall -> getDissipation()*2.0/7.0);

      speciesWall -> setTorsionFrictionCoefficient(0.0);
      speciesWall -> setTorsionStiffness(wallStiffness*2.0/7.0);
      speciesWall -> setTorsionDissipation(speciesWall -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesWall);

      // BIG-WALL
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + wallStiffness), bigWallRestitutionCoeff, particleMassBig);
      // plastic-adhesive part
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setCohesionStiffness(particleCohesiveStiffness);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setPenetrationDepthMax(particlePlasticityDepth);
      //        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setAdhesionStiffness(particleAdhesionStiffness);
      //        speciesHandler.getMixedObject(speciesBig, speciesWall) -> setAdhesionForceMax(particleMaxAdhesiveForce);

      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingFrictionCoefficient(bigWallSlidingFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);

      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingFrictionCoefficient(bigWallRollingFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);

      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionFrictionCoefficient(bigWallTorsionFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionStiffness(0.5*(particleStiffnessBig + wallStiffness)*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation()*2.0/7.0);

      // SMALL-WALL
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessSmall + wallStiffness), smallWallRestitutionCoeff, particleMassSmall);
      // plastic-adhesive part
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setCohesionStiffness(particleCohesiveStiffness);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setPenetrationDepthMax(particlePlasticityDepth);
      //        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setAdhesionStiffness(particleAdhesionStiffness);
      //        speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setAdhesionForceMax(particleMaxAdhesiveForce);

      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingFrictionCoefficient(smallWallSlidingFrictionCoeff);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);

      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingFrictionCoefficient(smallWallRollingFrictionCoeff);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);

      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionFrictionCoefficient(smallWallTorsionFrictionCoeff);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionStiffness(0.5*(particleStiffnessSmall + wallStiffness)*2.0/7.0);
      speciesHandler.getMixedObject(speciesSmall, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation()*2.0/7.0);

      // BIG-SMALL
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffnessBig + particleStiffnessSmall), bigSmallRestitutionCoeff, 0.5*(particleMassBig + particleMassSmall));
      // plastic-adhesive part
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setUnloadingStiffnessMax(particleMaxUnloadingStiffness);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setCohesionStiffness(particleCohesiveStiffness);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setPenetrationDepthMax(particlePlasticityDepth);
      //        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setAdhesionStiffness(particleAdhesionStiffness);
      //        speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setAdhesionForceMax(particleMaxAdhesiveForce);

      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingFrictionCoefficient(bigSmallSlidingFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);

      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingFrictionCoefficient(bigSmallRollingFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);

      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionFrictionCoefficient(bigSmallTorsionFrictionCoeff);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionStiffness(0.5*(particleStiffnessBig + particleStiffnessSmall)*2.0/7.0);
      speciesHandler.getMixedObject(speciesBig, speciesSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation()*2.0/7.0);

      if (verbose)
      {
         if (biModal)
         {
            std::cout << "\tBIG-BIG stiffness and dissipation coefficients: " << speciesBig -> getLoadingStiffness() << " " << speciesBig -> getDissipation() << "\n";
            std::cout << "\tBIG-BIG friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
            std::cout << "\tBIG-BIG tangential stiffnesses: " << speciesBig -> getSlidingStiffness() << " " << speciesBig -> getRollingStiffness() << " " << speciesBig -> getTorsionStiffness() << "\n";
            std::cout << "\tBIG-BIG tangential dissipation coefficients: " << speciesBig -> getSlidingDissipation() << " " << speciesBig -> getRollingDissipation() << " " << speciesBig -> getTorsionDissipation() << "\n";
            std::cout << "\tBIG-BIG collision time: " << std::setprecision(4) << speciesBig -> getCollisionTime(particleMassBig) << "\n\n";

            std::cout << "\tSMALL-SMALL stiffness and dissipation coefficients: " << speciesSmall -> getLoadingStiffness() << " " << speciesSmall -> getDissipation() << "\n";
            std::cout << "\tSMALL-SMALL friction coefficients: " << smallSmallSlidingFrictionCoeff << " " << smallSmallRollingFrictionCoeff << " " << smallSmallTorsionFrictionCoeff << "\n";
            std::cout << "\tSMALL-SMALL tangential stiffnesses: " << speciesSmall -> getSlidingStiffness() << " " << speciesSmall -> getRollingStiffness() << " " << speciesSmall -> getTorsionStiffness() << "\n";
            std::cout << "\tSMALL-SMALL tangential dissipation coefficients: " << speciesSmall -> getSlidingDissipation() << " " << speciesSmall -> getRollingDissipation() << " " << speciesSmall -> getTorsionDissipation() << "\n";
            std::cout << "\tSMALL-SMALL collision time: " << std::setprecision(4) << speciesSmall -> getCollisionTime(particleMassSmall) << "\n\n";

            std::cout << "\tBIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation() << "\n";
            std::cout << "\tBIG-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
            std::cout << "\tBIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionStiffness() << "\n";
            std::cout << "\tBIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionDissipation() << "\n";
            std::cout << "\tBIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(particleMassBig) << "\n\n";

            std::cout << "\tSMALL-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getDissipation() << "\n";
            std::cout << "\tSMALL-WALL friction coefficients: " << smallWallSlidingFrictionCoeff << " " << smallWallRollingFrictionCoeff << " " << smallWallTorsionFrictionCoeff << "\n";
            std::cout << "\tSMALL-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getTorsionStiffness() << "\n";
            std::cout << "\tSMALL-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getTorsionDissipation() << "\n";
            std::cout << "\tSMALL-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getCollisionTime(particleMassSmall) << "\n\n";

            std::cout << "\tBIG-SMALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getDissipation() << "\n";
            std::cout << "\tBIG-SMALL friction coefficients: " << bigSmallSlidingFrictionCoeff << " " << bigSmallRollingFrictionCoeff << " " << bigSmallTorsionFrictionCoeff << "\n";
            std::cout << "\tBIG-SMALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getTorsionStiffness() << "\n";
            std::cout << "\tBIG-SMALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getTorsionDissipation() << "\n";
            std::cout << "\tBIG-SMALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getCollisionTime(0.5*(particleMassBig + particleMassSmall)) << "\n\n";

            std::cout << "\tBIG-BIG collision time / TIME STEP: " << std::setprecision(4) << (speciesBig -> getCollisionTime(particleMassBig))/getTimeStep() << "\n";
            std::cout << "\tSMALL-SMALL collision time / TIME STEP: " << std::setprecision(4) << (speciesSmall -> getCollisionTime(particleMassSmall))/getTimeStep() << "\n";
            std::cout << "\tBIG-WALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(particleMassBig))/getTimeStep() << "\n";
            std::cout << "\tSMALL-WALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(speciesSmall, speciesWall) -> getCollisionTime(particleMassSmall))/getTimeStep() << "\n";
            std::cout << "\tBIG-SMALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(speciesBig, speciesSmall) -> getCollisionTime(0.5*(particleMassBig + particleMassSmall)))/getTimeStep() << "\n\n";
         }
         else
         {
            std::cout << "\tPARTICLE stiffness and dissipation coefficients: " << speciesBig -> getLoadingStiffness() << " " << speciesBig -> getDissipation() << "\n";
            std::cout << "\tPARTICLE friction coefficients: " << bigBigSlidingFrictionCoeff << " " << bigBigRollingFrictionCoeff << " " << bigBigTorsionFrictionCoeff << "\n";
            std::cout << "\tPARTICLE tangential stiffnesses: " << speciesBig -> getSlidingStiffness() << " " << speciesBig -> getRollingStiffness() << " " << speciesBig -> getTorsionStiffness() << "\n";
            std::cout << "\tPARTICLE tangential dissipation coefficients: " << speciesBig -> getSlidingDissipation() << " " << speciesBig -> getRollingDissipation() << " " << speciesBig -> getTorsionDissipation() << "\n";
            std::cout << "\tPARTICLE collision time: " << std::setprecision(4) << speciesBig -> getCollisionTime(particleMassBig) << "\n\n";

            std::cout << "\tPARTICLE-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getDissipation() << "\n";
            std::cout << "\tPARTICLE-WALL friction coefficients: " << bigWallSlidingFrictionCoeff << " " << bigWallRollingFrictionCoeff << " " << bigWallTorsionFrictionCoeff << "\n";
            std::cout << "\tPARTICLE-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionStiffness() << "\n";
            std::cout << "\tPARTICLE-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getTorsionDissipation() << "\n";
            std::cout << "\tPARTICLE-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(particleMassBig) << "\n\n";

            std::cout << "\tPARTICLE collision time / TIME STEP: " << std::setprecision(4) << (speciesBig -> getCollisionTime(particleMassBig))/getTimeStep() << "\n";
            std::cout << "\tPARTICLE-WALL collision time / TIME STEP: " << std::setprecision(4) << (speciesHandler.getMixedObject(speciesBig, speciesWall) -> getCollisionTime(particleMassBig))/getTimeStep() << "\n\n";
         }
      }
   }

   void insertBiggestParticle()
   {
      Vec3D particlePosition;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setRadius(meanRadiusBig*(1.0 + dispersityBig));
      p0.setSpecies(speciesBig);

      particlePosition.X = (casingRadius - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
      particlePosition.Y = sqrt(pow(casingRadius - 1.01*(p0.getRadius()), 2.0) - pow(particlePosition.X, 2.0))*random.getRandomNumber(-1.0,1.0);
      particlePosition.Z = random.getRandomNumber(1.01*p0.getRadius(), casingHeight - 1.01*(p0.getRadius()));
      p0.setPosition(particlePosition);

      particleHandler.copyAndAddObject(p0);

      nBigInserted++;
   }

   void computeNumberOfParticles()
   {
      nBigParticles = (int)(1.1*bulkPackingFractionPreCompression*casingVolume/particleVolumeBig);

      std::cout << "\tNumber of BIG particles needed: " << nBigParticles << "\n";
   }

   bool particleInsertionSuccessful(bool isBig)
   {
      int insertionFailCounter = 0;
      Vec3D particlePosition;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));

      if (isBig)
      {
         p0.setRadius(meanRadiusBig*(1.0 + dispersityBig*random.getRandomNumber(-1.0,1.0)));
         p0.setSpecies(speciesBig);
      }
      else
      {
         p0.setRadius(meanRadiusSmall*(1.0 + dispersitySmall*random.getRandomNumber(-1.0,1.0)));
         p0.setSpecies(speciesSmall);
      }

      do
      {
         particlePosition.X = (casingRadius - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Y = sqrt(pow(casingRadius - 1.01*(p0.getRadius()), 2.0) - pow(particlePosition.X, 2.0))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Z = random.getRandomNumber(casingHeight + 1.01*(p0.getRadius()),totalHeight - 1.01*(p0.getRadius()));

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
      while(!waitForParticleSettling && (nBigInserted < nBigParticles || nSmallInserted < nSmallParticles))
      {
         if (nBigInserted < nBigParticles && nSmallInserted < nSmallParticles)
         {
            if (random.getRandomNumber(0.0,1.0) < nBigParticles/(nBigParticles + nSmallParticles))
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
         else if (nBigInserted < nBigParticles && nSmallInserted >= nSmallParticles)
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
         else if (nBigInserted >= nBigParticles && nSmallInserted < nSmallParticles)
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

      if (waitForParticleSettling && getTime() > t0 + sqrt(2.0*(totalHeight - casingHeight)/9.81)) waitForParticleSettling = false;
      if (!waitForParticleSettling && nBigInserted >= nBigParticles && nSmallInserted >= nSmallParticles && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
      {
         std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
         stage++;
      }
   }








   // computes the total volume of all the loaded particles
   void computeParticleTotalVolume()
   {
      particleTotalVolume = nBigParticles*particleVolumeBig + nSmallParticles*particleVolumeSmall;

      if (verbose)
      {
         std::cout << "\tTotal volume of particles: " << particleTotalVolume << "\n";
         std::cout << "\tSanity check (nSets*volumePerSet)/(bulkPackingFraction*casingVolume): " << particleTotalVolume/(bulkPackingFractionPreCompression*casingVolume) << "\n";
      }
   }

   // computes the particle bed height and the filling region positioning
   void computeHeights()
   {
      double settledBedHeight = particleTotalVolume/(bulkPackingFractionPreCompression*constants::pi*pow(casingRadius,2.0));

      // the total height of the casing is set to be bed_height + 1.5*casing_height
      totalHeight = settledBedHeight + 1.5*casingHeight;

      if (verbose)
      {
         std::cout << "\tThe settled powder bed height is assumed to be: " << settledBedHeight << "\n";
         std::cout << "\tThe total height of the system is set to: " << totalHeight << "\n";
         std::cout << "\tSanity check: settledBedHeight/casingHeight = " << settledBedHeight/casingHeight << "\n";
      }
   }

   // sets the simulation domain and the walls
   void setBoundaries()
   {
      // for the X-Y values 1.1*casingRadius is used, for the Z boundaries the extrema + 0.1*casingRadius is used
      setXMin(-1.1*casingRadius);
      setYMin(-1.1*casingRadius);
      setZMin(-0.1*casingRadius);

      setXMax(1.1*casingRadius);
      setYMax(1.1*casingRadius);
      setZMax(totalHeight + 0.1*casingRadius);

      if (verbose)
      {
         std::cout << "\tSimulation domain minimum point: (" << getXMin() << " , " << getYMin() << " , " << getZMin() << ")\n";
         std::cout << "\tSimulation domain maximum point: (" << getXMax() << " , " << getYMax() << " , " << getZMax() << ")\n";
      }
   }

   // makes the geometric components
   void makeGeometry()
   {
      wallHandler.clear();

      // the basis of the compaction chamber
      basis.setSpecies(speciesWall);
      basis.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
      wallHandler.copyAndAddObject(basis);

      // the roof of the compaction chamber
      roof.setSpecies(speciesWall);
      roof.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,totalHeight));
      wallHandler.copyAndAddObject(roof);

      // the external casing
      AxisymmetricIntersectionOfWalls(casing);
      casing.setSpecies(speciesWall);
      casing.setPosition(Vec3D(0.0,0.0,0.0));
      casing.setOrientation(Vec3D(0.0,0.0,1.0));
      casing.addObject(Vec3D(1.0,0.0,0.0),Vec3D(casingRadius,0.0,0.0));
      wallHandler.copyAndAddObject(casing);
   }

   // removes the particles laying above the casing height
   void removeParticles()
   {
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         if (particleHandler.getObject(i) -> getPosition().Z + particleHandler.getObject(i) -> getRadius() > casingHeight) particleHandler.removeObject(i);
      }
   }

   // gives teh particles a random small velocity to pack the system more densely
   // the energy given is evaluated according to Ekin = constant*(Eel + Egrav)
   // Ekin = 0.5*(N_big*m_big + N_small*m_small)*v_mean^2 = C*(Egrav - Eel) -> v_mean = sqrt(2*C*(Egrav - Eel)/(N_big*m_big + N_small*m_small))
   void thermalizeParticles()
   {
      double velocityPerParticle = sqrt(2.0*thermalizationEnergyScalingFactor*fabs(getGravitationalEnergy() - getElasticEnergy())/(nBigParticles*particleMassBig + nSmallParticles*particleMassSmall));

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         Vec3D thermalVelocity;
         thermalVelocity.X = random.getRandomNumber(-1.0,1.0);
         thermalVelocity.Y = random.getRandomNumber(-1.0,1.0);
         thermalVelocity.Z = random.getRandomNumber(-1.0,1.0);

         thermalVelocity *= velocityPerParticle/thermalVelocity.getLength();
         particleHandler.getObject(i) -> setVelocity(thermalVelocity);
      }

      if (verbose)
      {
         std::cout << "\tParticle velocity modulus added = " << velocityPerParticle << "\n";
      }
   }

   // computes teh mean coordination number and teh packing fraction
   void makeDataAnalysis()
   {
      int particlesBelowCasingHeight;
      particlesBelowCasingHeight = 0;
      meanCoordinationNumber = 0.0;
      packingFraction = 0.0;

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         if (particleHandler.getObject(i) -> getPosition().Z + particleHandler.getObject(i) -> getRadius() < casingHeight)
         {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
            packingFraction += pow(particleHandler.getObject(i) -> getRadius(),3.0);
            particlesBelowCasingHeight++;
         }
      }

      meanCoordinationNumber /= particlesBelowCasingHeight;
      packingFraction /= 3.0*casingHeight*pow(casingRadius,2.0)/4.0;
   }

   // creates the data output file and writes the first row
   void makeOutputFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      outputFile.open(cdatName.str(), std::ios::out);
      outputFile << "time \t Eel \t Ekin/Eel \t packing fraction \t mean coord. number" << std::endl;
   }

   // writes the compression data to the output file
   void writeDataToOutptFile()
   {
      outputFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      packingFraction << "   " <<
      meanCoordinationNumber << "   " <<
      std::endl;
   }


   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {

      std::cout << "t = " << std::setprecision(3) << std::left << std::setw(3) << getTime() << ", tmax = " << getTimeMax() <<
      ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() <<
      ", p.f. = " << std::setprecision(3) << std::left << std::setw(4) << packingFraction << ", mean coord. number = " << meanCoordinationNumber <<
      ", nB/nBtot = " << nBigInserted/nBigParticles << ", nS/nStot = " << nSmallInserted/nSmallParticles <<
      ", n/nTot = " << particleHandler.getNumberOfObjects()/(nBigParticles + nSmallParticles) << std::endl;

      std::cout.flush();
   }

   bool continueSolve() const override
   {
      if (stage == 4 && getKineticEnergy()/getElasticEnergy() < 1.e-10) return false;
      return true;
   }


   // ----- VARIABLES -----
   // particle intrinsic properties
   double particleDensityBig, particleDensitySmall;
   double meanRadiusBig, meanRadiusSmall;
   double dispersityBig, dispersitySmall;
   double particleVolumeBig, particleVolumeSmall;
   double particleMassBig, particleMassSmall;
   double particleTotalVolume;

   // particle interaction properties
   double wallStiffness;
   double particleStiffnessBig, particleStiffnessSmall;
   double bigWallRestitutionCoeff, smallWallRestitutionCoeff;
   double bigBigRestitutionCoeff, smallSmallRestitutionCoeff, bigSmallRestitutionCoeff;
   // particle-wall friction coefficients
   double bigWallSlidingFrictionCoeff, smallWallSlidingFrictionCoeff;
   double bigWallRollingFrictionCoeff, smallWallRollingFrictionCoeff;
   double bigWallTorsionFrictionCoeff, smallWallTorsionFrictionCoeff;
   // particle-particle friction coefficients
   double bigBigSlidingFrictionCoeff, smallSmallSlidingFrictionCoeff, bigSmallSlidingFrictionCoeff;
   double bigBigRollingFrictionCoeff, smallSmallRollingFrictionCoeff, bigSmallRollingFrictionCoeff;
   double bigBigTorsionFrictionCoeff, smallSmallTorsionFrictionCoeff, bigSmallTorsionFrictionCoeff;
   // particle plastic and adhesive coefficients   // *** CARE: NO SEPARATE CASES FOR BIG AND SMALL YET! ***
   double particleMaxUnloadingStiffness;
   double particleCohesiveStiffness;
   double particlePlasticityDepth;
   //    double particleAdhesionStiffness;
   //    double particleMaxAdhesiveForce;

   // static geometry related variables
   double casingRadius;
   double casingHeight;
   double casingVolume;
   double densityWall;
   double totalHeight;

   // simulation variables
   int stage;
   double nBigInserted, nSmallInserted;
   double nBigParticles, nSmallParticles;
   double smallToBigMassRatio;
   double bulkPackingFractionPreCompression;
   double timePlaceholder;
   int numberOfThermalizationInstances;
   double thermalizationTimeInterval;
   double thermalizationEnergyScalingFactor;
   double energyRatioTolerance;

   // initialization variables
   double initRegionHeight;
   double volumeOfParticlesPerSet;
   int setsInserted;

   // parameter study variables
   double meanCoordinationNumber;
   double packingFraction;
   double maximumPackingFraction;
   int thermalizationInstancesDone;
   int thermalizationCyclesCounter;
   int packingFractionCyclesThreshold;
   double packingFractionDifferenceTolerance;

   // Mercury-specific variables
   //    LinearViscoelasticFrictionSpecies *speciesBig, *speciesSmall, *speciesWall;
   LinearPlasticViscoelasticFrictionSpecies *speciesBig, *speciesSmall, *speciesWall;
   InfiniteWall basis, roof;
   AxisymmetricIntersectionOfWalls casing;
   BaseParticle p0;

   // global variables
   bool verbose;
   bool biModal;
   bool waitForParticleSettling;
   double t0;
   double cdatOutputTimeInterval;
   std::ofstream outputFile;
};

int main(int argc, char *argv[])
{
   // geometry variables
   double cylinderRadius = 0.0125;
   double cylinderHeight = 0.019;

   // particle variables
   double casingToParticleSizeRatio = 40.0;
   double particleRadius = cylinderRadius/casingToParticleSizeRatio;
   double particleSizeDispersity = 0.10;

   // interaction parameters
   double k1Particle = 500.0;
   double k2MaxParticle = 5000.0;
   double kCParticle = 0.0;
   double phiParticle = 0.05;
   double k1Wall = k2MaxParticle;

   double particleDensity = 1452.7;
   double wallDensity = 5000.0;
   double ppRestitutionCoeff = 0.5;
   double pwRestitutionCoeff = 0.7;

   // friction parameters
   double ppSlidingFriction = 0.16;
   double ppRollingFriction = 0.05;
   double ppTorsionFriction = 0.0;

   double pwSlidingFriction = 0.30;
   double pwRollingFriction = 0.01;
   double pwTorsionFriction = 0.0;

   // setup
   CalibrationRoutine_ParticleSettling crSettling;
   crSettling.setName("CompressionTest_SETTLING");
   crSettling.setCdatOutputTimeInterval(0.005);

   // sets simulation parameters
   crSettling.setTimeStep(1.0e-6);
   crSettling.setTimeMax(100.0);
   crSettling.setGravity(Vec3D(0., 0., -9.81));
   crSettling.setSystemDimensions(3);
   crSettling.setVerbose(true);

   // sets the number of saved timesteps such that the output is printed every 0.01s
   crSettling.setSaveCount(0.01/crSettling.getTimeStep());

   // sets the particle intrinsic properties
   crSettling.setParticleDensity(particleDensity);
   crSettling.setParticleRadiusAndDispersity(particleRadius, particleSizeDispersity);

   // sets the small-to-big total mass ratio
   crSettling.setTotalMassRatio(0.0);

   // sets the stiffnesses
   crSettling.setWallStiffness(k1Wall);
   //    crSettling.setParticleStiffness(1000.0);
   crSettling.setParticlePlasticProperties(k1Particle, k2MaxParticle, kCParticle, phiParticle);
   //    crSettling.setParticlesAdhesionProperties(0.0, 0.0);

   // sets the particle-wall friction coefficients
   crSettling.setParticleWallSlidingFrictionCoeff(pwSlidingFriction);
   crSettling.setParticleWallRollingFrictionCoeff(pwRollingFriction);
   crSettling.setParticleWallTorsionFrictionCoeff(pwTorsionFriction);

   // sets the particle-particle friction coefficients
   crSettling.setParticleParticleSlidingFrictionCoeff(ppSlidingFriction);
   crSettling.setParticleParticleRollingFrictionCoeff(ppRollingFriction);
   crSettling.setParticleParticleTorsionFrictionCoeff(ppTorsionFriction);

   // sets the particle restitution coefficients
   crSettling.setParticleWallRestitutionCoeff(pwRestitutionCoeff);
   crSettling.setParticleParticleRestitutionCoeff(ppRestitutionCoeff);

   // sets the casing dimensions
   crSettling.setCasingProperties(cylinderRadius, cylinderHeight, wallDensity);

   // sets the bulk packing fraction of the powder prior to compression
   crSettling.setbulkPackingFractionPreCompression(0.66);

   // sets the energy ratio tolerance to determine if the system is static
   crSettling.setEnergyRatioTolerance(1.0e-4);

   // sets the tolerance for the recognization of the max packing fraction and teh number of cycles for that value to be held
   crSettling.setPackingFractionDifferenceTolerance(1.0e-3, 5);

   // sets additional built-in arguments for the xballs visualization
   crSettling.setXBallsAdditionalArguments("-h 800 -p 10 -o 200 -3dturn 3");

   crSettling.solve();

   return 0;
}
