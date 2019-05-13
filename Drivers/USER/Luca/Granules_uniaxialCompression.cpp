/*
 *** GRANULES UNIAXIAL COMPRESSION ***
 A single granule is created and its properties computed
 Gravity is activated and the agglomerate left to settle
 The agglomerate is compressed uniaxially up to a desired force
 */

#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <fstream>

/*
 ToDo:
 -

 27.11:
 - added the target final deformation (height-driven) compression case
 - added a standby phase instead of the unloading one

 */

class GranulesUniaxialCompression : public Mercury3D
{
private:
   void setupInitialConditions() override
   {
      stage = 0;
      nParticlesInserted = 0;
      numberOfBindingInteractions = 0;
      forceModulus = 0.0;
      t0 = getTime();

      std::cout << "SETTING SPECIES VECTOR..." << std::endl;
      setSpeciesVector();

      std::cout << "SETTING MIXED SPECIES MATRIX..." << std::endl;
      setMixedSpeciesCohesionMatrix();

      std::cout << "SETTING CUBIC LATTICE SIZE..." << std::endl;
      setBoxSize();

      std::cout << "SETTING DOMAIN LIMITS... " << std::endl;
      setDomainLimits();

      // std::cout << "CREATING BOUNDARIES..." << std::endl;
      // makeBoundaries();

      std::cout << std::endl << "PARTICLE INSERTION" << std::endl;
      insertParticles();

      std::cout << std::endl << "COMPUTING PARTICLES VOLUME" << std::endl;
      computeTotalParticleVolume();

      std::cout << "CREATING .gdat FILE" << std::endl << std::endl;
      makeGdatFile();

      std::cout << "ACTIVATING CENTRAL FORCES" << std::endl;

      stage++;
   }

   void actionsOnRestart() override
   {

   }

   void actionsAfterTimeStep() override
   {
      // DATA ANALYSIS
      if (stage >= 1) makeDataAnalysis();
      if (stage >= 1 && fmod(getTime(), 0.001) < getTimeStep()) writeToGdatFile();

      // ACTIVATION OF LOCAL FORCE
      if (stage == 1)
      {
         applyCentralForce();
         if (getTime() - t0 < gravityAndForceTuningDuration)
         {
            if (fmod(getTime() - t0,forceTuningInterval) < getTimeStep()) increaseForce();
            if (fmod(getTime() - t0,velocityDampingInterval) < getTimeStep()) dampVelocities();
         }
         else
         {
            std::cout << "DAMPING CENTRAL FORCE" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // DAMPING OF LOCAL FORCE
      if (stage == 2)
      {
         applyCentralForce();
         if (getTime() - t0 < gravityAndForceTuningDuration)
         {
            if (fmod(getTime() - t0,forceTuningInterval) < getTimeStep()) dampForce();
            if (fmod(getTime() - t0,velocityDampingInterval) < getTimeStep()) dampVelocities();
         }
         else
         {
            std::cout << "CREATING COORDINATION MATRIX" << std::endl << std::endl;
            setupCoordinationMatrix();

            std::cout << "DISSIPATING ENERGY" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // ENERGY DISSIPATION AND CREATION OF THE BASE
      if (stage == 3)
      {
         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance || getTime() - t0 > gravityAndForceTuningDuration)
         {
            std::cout << "ENERGY DISSIPATED. COMPUTING INTERNAL STRUCTURE" << std::endl << std::endl;
            // computeInternalStructure();

            std::cout << "CREATING BASE AND HOUSING" << std::endl << std::endl;
            makeBaseAndHousing();

            std::cout << "ACTIVATING GRAVITY" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // ACTIVATING GRAVITY, SETTLING AND CREATION OF THE PISTON
      if (stage == 4)
      {
         if (getTime() - t0 < gravityAndForceTuningDuration)
         {
            activateGravity();
         }

         if (getTime() - t0 > gravityAndForceTuningDuration && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
         {
            std::cout << "CREATING PISTON AND STARTING COMPRESSION" << std::endl << std::endl;
            makeCompressionPiston();

            stage++;
         }
      }

      // COMPRESSION CYCLE
      if (stage == 5)
      {
         // if (pistonForce < maximumCompressiveForce) // uncomment this for force-driven compression
         if (pistonHeight - baseHeight > initialClusterSize*minimumRelativeHeightCompression) // uncomment this for height-driven compression
         {
            loadPiston();
         }
         else
         {
            // std::cout << "MAXIMUM COMPRESSION REACHED. UNLOADING PISTON" << std::endl << std::endl;
            std::cout << "MAXIMUM COMPRESSION REACHED. STANDBY PHASE" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // UNLOAD CYCLE
      if (stage == 6)
      {
         // if (pistonForce > energyRatioTolerance)
         // {
         //    unloadPiston();
         // }
         if (getTime() - t0 < 0.2) {}
         else
         {
            std::cout << "PISTON UNLOADED. QUITTING" << std::endl << std::endl;
            setTimeMax(getTime() + getTimeStep());
            stage++;
         }
      }
   }

   void actionsAfterSolve() override
   {
      gdatFile.close();
   }


public:
   // FUNCTIONS CALLED IN MAIN ----------------------------------------
   void setCdatOutputTimeInterval(double dt)
   {
      gdatOutputTimeInterval = dt;
   }

   void setEnergyRatioTolerance(double eRatio)
   {
      energyRatioTolerance = eRatio;
   }

   void setGranulesProperties(int nP, double sizeRatio)
   {
      totalNumberOfParticles = nP;
      granuleSizeToLatticeRatio = sizeRatio;
   }

   void setParticleProperties(double rP, double dP, double rhoP)
   {
      radiusParticle = rP;
      sizeDispersityParticle = dP;
      densityParticle = rhoP;

      volumeParticle = 4.0*constants::pi*pow(radiusParticle,3.0)/3.0;
      massParticle = volumeParticle*densityParticle;
   }

   void setWallDensity(double rhoW)
   {
      densityWall = rhoW;
   }

   void setParticleParticleFrictionCoefficients(double muSliding, double muRolling, double muTorsion)
   {
      particleParticleSlidingFrictionCoeff = muSliding;
      particleParticleRollingFrictionCoeff = muRolling;
      particleParticleTorsionFrictionCoeff = muTorsion;
   }

   void setParticleWallFrictionCoefficients(double muSliding, double muRolling, double muTorsion)
   {
      particleWallSlidingFrictionCoeff = muSliding;
      particleWallRollingFrictionCoeff = muRolling;
      particleWallTorsionFrictionCoeff = muTorsion;
   }

   void setWallStiffnessAndRestitutionCoefficients(double kW, double eB)
   {
      kpWall = kW;
      particleWallRestitutionCoeff = eB;
   }

   void setParticleParticleRestitutionCoefficients(double bigBigE)
   {
      particleParticleRestitutionCoeff = bigBigE;
   }

   void setParticlePlasticProperties(double k1, double k2max, double kCGranule, double kCLoose, double phi)
   {
      kpParticle = k1;
      keParticle = k2max;
      kcParticleIntergranule = kCGranule;
      kcParticleLoose = kCLoose;
      phiParticle = phi;
   }

   void setForceAndVelocityProperties(double fmm, double ffm, double fdi, double vdm, double vdi)
   {
      maximumForceModulus = fmm;
      finalForceModulus = ffm;
      forceTuningInterval = fdi;
      velocityDampingModulus = vdm;
      velocityDampingInterval = vdi;

      std::cout << "Maximum force modulus: " << maximumForceModulus << std::endl;
      std::cout << "Final force modulus: " << finalForceModulus << std::endl;
   }

   void setGravityAndForceTuningDuration(double dt)
   {
      gravityAndForceTuningDuration = dt;
   }

   void setInternalStructureAnalysisParameters(int length, double ratio)
   {
      gridLength = length;
      fictiousGridPointRadiusRatio = ratio;
   }

   void setCompressionParameters(double maxForce, double minRatio)
   {
      maximumCompressiveForce = maxForce;
      minimumRelativeHeightCompression = minRatio;
   }


   // FUNCTIONS CALLED IN THE CLASS -----------------------------------
   void setSpecies()
   {
      speciesHandler.clear();

      // BIG-BIG
      speciesParticle = new LinearPlasticViscoelasticFrictionSpecies;
      speciesParticle -> setDensity(densityParticle);
      speciesParticle -> setStiffnessAndRestitutionCoefficient(kpParticle, particleParticleRestitutionCoeff, massParticle);
      speciesParticle -> setUnloadingStiffnessMax(keParticle);
      speciesParticle -> setCohesionStiffness(kcParticleIntergranule);
      speciesParticle -> setPenetrationDepthMax(phiParticle);

      speciesParticle -> setSlidingFrictionCoefficient(particleParticleSlidingFrictionCoeff);
      speciesParticle -> setSlidingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setSlidingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setRollingFrictionCoefficient(particleParticleRollingFrictionCoeff);
      speciesParticle -> setRollingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setRollingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setTorsionFrictionCoefficient(particleParticleTorsionFrictionCoeff);
      speciesParticle -> setTorsionStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setTorsionDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesParticle);

      // WALL-WALL
      speciesWall = new LinearPlasticViscoelasticFrictionSpecies;
      speciesWall -> setDensity(densityWall);
      speciesWall -> setStiffnessAndRestitutionCoefficient(kpWall, 1.0, massParticle);
      speciesWall -> setUnloadingStiffnessMax(kpWall);
      speciesWall -> setCohesionStiffness(0.0);
      speciesWall -> setPenetrationDepthMax(0.001);

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
      //speciesMixedBigWall = speciesHandler.getMixedObject(speciesParticle, speciesWall);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(kpParticle + kpWall), particleWallRestitutionCoeff, massParticle);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setUnloadingStiffnessMax(keParticle);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setCohesionStiffness(0.0);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setPenetrationDepthMax(0.001);

      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setSlidingFrictionCoefficient(particleWallSlidingFrictionCoeff);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setRollingFrictionCoefficient(particleWallRollingFrictionCoeff);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setTorsionFrictionCoefficient(particleWallTorsionFrictionCoeff);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getDissipation()*2.0/7.0);

      // OUTPUTS
      // std::cout << "BIG-BIG stiffness and dissipation coefficients: " << speciesParticle -> getLoadingStiffness() << " " << speciesParticle -> getDissipation() << "\n";
      // std::cout << "BIG-BIG friction coefficients: " << particleParticleSlidingFrictionCoeff << " " << particleParticleRollingFrictionCoeff << " " << particleParticleTorsionFrictionCoeff << "\n";
      // std::cout << "BIG-BIG tangential stiffnesses: " << speciesParticle -> getSlidingStiffness() << " " << speciesParticle -> getRollingStiffness() << " " << speciesParticle -> getTorsionStiffness() << "\n";
      // std::cout << "BIG-BIG tangential dissipation coefficients: " << speciesParticle -> getSlidingDissipation() << " " << speciesParticle -> getRollingDissipation() << " " << speciesParticle -> getTorsionDissipation() << "\n";
      // std::cout << "BIG-BIG collision time: " << std::setprecision(4) << speciesParticle -> getCollisionTime(massParticle) << "\n\n";
      //
      // std::cout << "BIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getDissipation() << "\n";
      // std::cout << "BIG-WALL friction coefficients: " << particleWallSlidingFrictionCoeff << " " << particleWallRollingFrictionCoeff << " " << particleWallTorsionFrictionCoeff << "\n";
      // std::cout << "BIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getTorsionStiffness() << "\n";
      // std::cout << "BIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getTorsionDissipation() << "\n";
      // std::cout << "BIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getCollisionTime(massParticle) << "\n\n";

      std::cout << "tC_BB/dt: " << std::setprecision(4) << speciesParticle -> getCollisionTime(massParticle)/getTimeStep() << "\n";
      std::cout << "tC_BW/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getCollisionTime(massParticle)/getTimeStep() << "\n";
   }

   void setSpeciesVector()
   {
      speciesHandler.clear();

      // WALL-WALL
      speciesWall = new LinearPlasticViscoelasticFrictionSpecies;
      speciesWall -> setDensity(densityWall);
      speciesWall -> setStiffnessAndRestitutionCoefficient(kpWall, 1.0, massParticle);
      speciesWall -> setUnloadingStiffnessMax(kpWall);
      speciesWall -> setCohesionStiffness(0.0);
      speciesWall -> setPenetrationDepthMax(0.001);

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

      // SINGLE_PARTICLE-SINGLE_PARTICLE
      particlePureSpeciesVector.reserve(totalNumberOfParticles);

      for (int i = 0; i < totalNumberOfParticles; i++)
      {
         particlePureSpeciesVector[i] = new LinearPlasticViscoelasticFrictionSpecies;
         particlePureSpeciesVector[i] -> setDensity(densityParticle);
         particlePureSpeciesVector[i] -> setStiffnessAndRestitutionCoefficient(kpParticle, particleParticleRestitutionCoeff, massParticle);
         particlePureSpeciesVector[i] -> setUnloadingStiffnessMax(keParticle);
         particlePureSpeciesVector[i] -> setCohesionStiffness(kcParticleIntergranule);
         particlePureSpeciesVector[i] -> setPenetrationDepthMax(phiParticle);

         particlePureSpeciesVector[i] -> setSlidingFrictionCoefficient(particleParticleSlidingFrictionCoeff);
         particlePureSpeciesVector[i] -> setSlidingStiffness(particlePureSpeciesVector[i] -> getLoadingStiffness()*2.0/7.0);
         particlePureSpeciesVector[i] -> setSlidingDissipation(particlePureSpeciesVector[i] -> getDissipation()*2.0/7.0);
         particlePureSpeciesVector[i] -> setRollingFrictionCoefficient(particleParticleRollingFrictionCoeff);
         particlePureSpeciesVector[i] -> setRollingStiffness(particlePureSpeciesVector[i] -> getLoadingStiffness()*2.0/7.0);
         particlePureSpeciesVector[i] -> setRollingDissipation(particlePureSpeciesVector[i] -> getDissipation()*2.0/7.0);
         particlePureSpeciesVector[i] -> setTorsionFrictionCoefficient(particleParticleTorsionFrictionCoeff);
         particlePureSpeciesVector[i] -> setTorsionStiffness(particlePureSpeciesVector[i] -> getLoadingStiffness()*2.0/7.0);
         particlePureSpeciesVector[i] -> setTorsionDissipation(particlePureSpeciesVector[i] -> getDissipation()*2.0/7.0);
         speciesHandler.addObject(particlePureSpeciesVector[i]);

         particlePureSpeciesVector.push_back(particlePureSpeciesVector[i]);
      }

      // OUTPUTS
      // std::cout << "BIG stiffness and dissipation coefficients: " << particlePureSpeciesVector[0] -> getLoadingStiffness() << " " << particlePureSpeciesVector[0] -> getDissipation() << "\n";
      // std::cout << "BIG friction coefficients: " << particleParticleSlidingFrictionCoeff << " " << particleParticleRollingFrictionCoeff << " " << particleParticleTorsionFrictionCoeff << "\n";
      // std::cout << "BIG tangential stiffnesses: " << particlePureSpeciesVector[0] -> getSlidingStiffness() << " " << particlePureSpeciesVector[0] -> getRollingStiffness() << " " << particlePureSpeciesVector[0] -> getTorsionStiffness() << "\n";
      // std::cout << "BIG tangential dissipation coefficients: " << particlePureSpeciesVector[0] -> getSlidingDissipation() << " " << particlePureSpeciesVector[0] -> getRollingDissipation() << " " << particlePureSpeciesVector[0] -> getTorsionDissipation() << "\n";
      // std::cout << "BIG collision time: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle) << "\n";
      std::cout << "tC_BB/dt: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle)/getTimeStep() << "\n\n";

      // std::cout << "WALL stiffness and dissipation coefficients: " << speciesWall -> getLoadingStiffness() << " " << speciesWall -> getDissipation() << "\n\n";
   }

   void setMixedSpeciesCohesionMatrix()
   {
      for (int i = 0; i < totalNumberOfParticles; i++)
      {
         std::vector<LinearPlasticViscoelasticFrictionMixedSpecies*> temporaryRowVector;
         for (int j = 0; j < totalNumberOfParticles; j++)
         {
            if (i == j)
            {
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(kpParticle + kpWall), particleWallRestitutionCoeff, massParticle);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setUnloadingStiffnessMax(keParticle);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setCohesionStiffness(0.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setPenetrationDepthMax(phiParticle);

               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setSlidingFrictionCoefficient(particleWallSlidingFrictionCoeff);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setSlidingStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> getLoadingStiffness()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setSlidingDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> getDissipation()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setRollingFrictionCoefficient(particleWallRollingFrictionCoeff);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setRollingStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> getLoadingStiffness()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setRollingDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> getDissipation()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setTorsionFrictionCoefficient(particleWallTorsionFrictionCoeff);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setTorsionStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> getLoadingStiffness()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setTorsionDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> getDissipation()*2.0/7.0);

               temporaryRowVector.push_back(speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall));
            }
            else
            {
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setStiffnessAndRestitutionCoefficient(0.5*(particlePureSpeciesVector[i] -> getLoadingStiffness() + particlePureSpeciesVector[j] -> getLoadingStiffness()), particleParticleRestitutionCoeff, massParticle);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setUnloadingStiffnessMax(keParticle);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setCohesionStiffness(kcParticleIntergranule);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setPenetrationDepthMax(phiParticle);

               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setSlidingFrictionCoefficient(particleParticleSlidingFrictionCoeff);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setSlidingStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getLoadingStiffness()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setSlidingDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getDissipation()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setRollingFrictionCoefficient(particleParticleRollingFrictionCoeff);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setRollingStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getLoadingStiffness()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setRollingDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getDissipation()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setTorsionFrictionCoefficient(particleParticleTorsionFrictionCoeff);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setTorsionStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getLoadingStiffness()*2.0/7.0);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setTorsionDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getDissipation()*2.0/7.0);

               temporaryRowVector.push_back(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]));
            }
         }

         mixedSpeciesMatrix.push_back(temporaryRowVector);
      }

      // OUTPUTS
      // std::cout << "BIG-WALL stiffness and dissipation coefficients: " << mixedSpeciesMatrix[0][0] -> getLoadingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getDissipation() << "\n";
      // std::cout << "BIG-WALL friction coefficients: " << particleWallSlidingFrictionCoeff << " " << particleWallRollingFrictionCoeff << " " << particleWallTorsionFrictionCoeff << "\n";
      // std::cout << "BIG-WALL tangential stiffnesses: " << mixedSpeciesMatrix[0][0] -> getSlidingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getRollingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getTorsionStiffness() << "\n";
      // std::cout << "BIG-WALL tangential dissipation coefficients: " << mixedSpeciesMatrix[0][0] -> getSlidingDissipation() << " " << mixedSpeciesMatrix[0][0] -> getRollingDissipation() << " " << mixedSpeciesMatrix[0][0] -> getTorsionDissipation() << "\n";
      // std::cout << "BIG-WALL collision time: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle) << "\n";
      std::cout << "tC_BW/dt: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle)/getTimeStep() << "\n\n";
   }

   void setBoxSize()
   {
      granuleRadiusPreFormation = radiusParticle*pow(totalNumberOfParticles/0.5,1.0/3.0);
      boxSize = 2.0*granuleRadiusPreFormation/granuleSizeToLatticeRatio;

      std::cout << "Granule average radius assuming packing fraction of 0.5 " << granuleRadiusPreFormation << std::endl;
      std::cout << "Cubic lattice size " << boxSize << std::endl << std::endl;
   }

   void setDomainLimits()
   {
      setXMin(-0.5*boxSize);
      setYMin(-0.5*boxSize);
      setZMin(-0.5*boxSize);

      setXMax(0.5*boxSize);
      setYMax(0.5*boxSize);
      setZMax(0.5*boxSize);

      pistonHeight = getZMax();
   }

   void makeBoundaries()
   {
      wallHandler.clear();
      wall.setSpecies(speciesWall);

      // x walls
      wall.set(Vec3D(-1.,0.,0.),Vec3D(getXMin(),0.,0.));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(1.,0.,0.),Vec3D(getXMax(),0.,0.));
      wallHandler.copyAndAddObject(wall);

      // y walls
      wall.set(Vec3D(0.,-1.,0.),Vec3D(0.,getYMin(),0.));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(0.,1.,0.),Vec3D(0.,getYMax(),0.));
      wallHandler.copyAndAddObject(wall);

      // z walls
      wall.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,getZMin()));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,getZMax()));
      wallHandler.copyAndAddObject(wall);
   }

   bool particleInsertionSuccessful()
   {
      int insertionFailCounter = 0;
      double rad, theta, phi;
      Vec3D particlePosition;
      SphericalParticle p0;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setRadius(radiusParticle*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.0,1.0)));
      p0.setSpecies(particlePureSpeciesVector[particleHandler.getNumberOfObjects()]);

      do
      {
         // particlePosition.X = (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         // particlePosition.Y = (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         // particlePosition.Z = (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);

         rad = random.getRandomNumber(p0.getRadius(),0.5*boxSize - 1.01*(p0.getRadius()));
         theta = constants::pi*random.getRandomNumber(-1.0,1.0);
         phi = 0.5*constants::pi*random.getRandomNumber(-1.0,1.0);
         particlePosition.X = rad*sin(theta)*cos(phi);
         particlePosition.Y = rad*sin(theta)*sin(phi);
         particlePosition.Z = rad*cos(theta);

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
      nParticlesInserted = 0;

      while(nParticlesInserted < totalNumberOfParticles)
      {
         if (particleInsertionSuccessful()) nParticlesInserted++;
         else t0 = getTime();
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
   }

   void applyCentralForce()
   {
      Vec3D distanceFromForceCenter;

      for (int i=0; i<totalNumberOfParticles; i++)
      {
         distanceFromForceCenter = particleHandler.getObject(i) -> getPosition();

         //particleHandler.getObject(i) -> addForce(-forceScaleModulus*Vec3D(0.2*distanceFromForceCenter.X,0.2*distanceFromForceCenter.Y,distanceFromForceCenter.Z));
         particleHandler.getObject(i) -> addForce(-forceModulus*distanceFromForceCenter);
         //particleHandler.getObject(i) -> addForce(-forceScaleModulus*distanceFromForceCenter*(particleHandler.getObject(i) -> getRadius())/(0.5*(radiusParticle + radiusSmall)));
         //particleHandler.getObject(i) -> addForce(-1.0*distanceFromForceCenter/(distanceFromForceCenter.getLength()*distanceFromForceCenter.getLength())/10000);
      }
   }

   void increaseForce()
   {
      forceModulus = maximumForceModulus*(getTime() - t0)/gravityAndForceTuningDuration;
   }

   void dampForce()
   {
      forceModulus = maximumForceModulus - (maximumForceModulus - finalForceModulus)*(getTime() - t0)/gravityAndForceTuningDuration;
   }

   void dampVelocities()
   {
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         particleHandler.getObject(i) -> setVelocity(velocityDampingModulus*(particleHandler.getObject(i) -> getVelocity()));
      }
   }

   void computeTotalParticleVolume()
   {
      totalParticleVolume = 0.0;
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--) totalParticleVolume += pow(particleHandler.getObject(i) -> getRadius(),3.0);
      totalParticleVolume *= 4.0*constants::pi/3.0;
   }

   double computeParticleVolume(int i)
   {
      return 4.0*constants::pi*pow(particleHandler.getObject(i) -> getRadius(), 3.0)/3.0;
   }

   // void computeCenterOfMass()
   // {
   //    centerOfMass.setZero();
   //    for (int i = 0; i < particleHandler.getNumberOfObjects(); i++) centerOfMass += computeParticleVolume(i)*(particleHandler.getObject(i) -> getPosition());
   //    centerOfMass /= totalParticleVolume;
   // }

   void makeDataAnalysis()
   {
      int totalInteractionCounter = 0;
      Vec3D distanceFromForceCenter;
      Vec3D localMin, localMax;
      localMin.setZero();
      localMax.setZero();
      centerOfMass.setZero();

      averageForceOnParticle = 0.0;
      meanGranuleRadius = 0.0;
      meanCoordinationNumber = 0.0;
      volumeRatio = 0.0;
      maxTotalRelativeOverlap = 0.0;
      meanTotalRelativeOverlap = 0.0;
      pistonForce = 0.0;

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
         centerOfMass += computeParticleVolume(i)*(particleHandler.getObject(i) -> getPosition());

         distanceFromForceCenter = particleHandler.getObject(i) -> getPosition();
         if (distanceFromForceCenter.X > localMax.X) localMax.X = distanceFromForceCenter.X;
         if (distanceFromForceCenter.X < localMin.X) localMin.X = distanceFromForceCenter.X;
         if (distanceFromForceCenter.Y > localMax.Y) localMax.Y = distanceFromForceCenter.Y;
         if (distanceFromForceCenter.Y < localMin.Y) localMin.Y = distanceFromForceCenter.Y;
         if (distanceFromForceCenter.Z > localMax.Z) localMax.Z = distanceFromForceCenter.Z;
         if (distanceFromForceCenter.Z < localMin.Z) localMin.Z = distanceFromForceCenter.Z;
      }
      meanCoordinationNumber /= particleHandler.getNumberOfObjects();
      meanGranuleRadius = (localMax.X - localMin.X + localMax.Y - localMin.Y + localMax.Z - localMin.Z)/6.0;
      volumeRatio = 4.0*constants::pi*pow(meanGranuleRadius,3.0)/(3.0*totalParticleVolume);
      centerOfMass /= totalParticleVolume;

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         averageForceOnParticle += ((*i) -> getForce()).getLength();

         meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

         totalInteractionCounter++;

         // piston interactions
         if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
         {
            pistonForce -= ((*i) -> getForce()).Z;
         }
      }
      averageForceOnParticle /= totalInteractionCounter;
      meanTotalRelativeOverlap /= totalInteractionCounter;
   }

   // creates the matrix containing the interaction information: "1" -> in contact, "0" -> not in contact
   void setupCoordinationMatrix()
   {
      // resets the number of binding interactions (the bounds after the granule formation)
      numberOfBindingInteractions = 0;

      // creates the matrix and fills it with zeroes
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         std::vector<int> temporaryRowVector;
         for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) temporaryRowVector.push_back(0);
         coordinationMatrix.push_back(temporaryRowVector);
         coordinationMatrixPreviousCheck.push_back(temporaryRowVector);
      }

      // now checks for interactions and allocates a "1" to corresponding interacting particle indices
      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         coordinationMatrix[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
         coordinationMatrix[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;
         coordinationMatrixPreviousCheck[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
         coordinationMatrixPreviousCheck[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;

         // updates the number of bounds
         numberOfBindingInteractions++;
      }

      // // print the coordination matrix
      // std::ofstream matrixOutput;
      // matrixOutput.open("Granules_INTERACTIONMATRIX", std::ios::out);
      // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      // {
      //    for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) matrixOutput << coordinationMatrix[i][j] << "\t";
      //    matrixOutput << std::endl;
      // }
      // matrixOutput.close();
   }

   // creates the data output file and writes the first row
   void makeGdatFile()
   {
      std::ostringstream gdatName;
      std::cout.unsetf(std::ios::floatfield);
      gdatName << getName() << ".gdat";

      gdatFile.open(gdatName.str(), std::ios::out);
      gdatFile << "time \t Eel \t Ekin/Eel \t coord_number \t n_bounds \t meanRadius \t Vgran/Vtot \t Fcomp \t Fmean \t dTotMean \t dTotMax \t cm_X \t cm_Y \t cm_Z \t Fpiston \t h_piston" << std::endl;
   }

   // writes the compression data to the output file
   void writeToGdatFile()
   {
      gdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      meanCoordinationNumber << "   " <<
      numberOfBindingInteractions << "   " <<
      meanGranuleRadius << "   " <<
      volumeRatio << "   " <<
      forceModulus << "   " <<
      averageForceOnParticle << "   " <<
      meanTotalRelativeOverlap << "   " <<
      maxTotalRelativeOverlap << "   " <<
      centerOfMass << "   " <<
      pistonForce << "   " <<
      pistonHeight << "   " <<
      std::endl;
   }

   // creates the grid file and writes the first row
   void makeGridFile()
   {
      std::ostringstream gridName;
      std::cout.unsetf(std::ios::floatfield);
      gridName << getName() << ".grid";

      gridFile.open(gridName.str(), std::ios::out);
      gridFile << "Grid_length: " << gridLength << "\t grid_resolution_length: " << gridResolutionLength << "\t fictiousGridParticleRadiusRatio: " << fictiousGridPointRadiusRatio << std::endl;
   }

   void computeInternalStructure()
   {
      Vec3D gridPoint;
      SphericalParticle p0;
      double nPointsInsideAgglomerateBoundary;
      double nPointsInsideComponents;

      gridResolutionLength = 2.0*meanGranuleRadius/(gridLength - 1.0);
      nPointsInsideAgglomerateBoundary = 0;
      nPointsInsideComponents = 0;

      // std::cout << "--- 1 ---" << std::endl;

      makeGridFile();
      p0.setRadius(radiusParticle*fictiousGridPointRadiusRatio);
      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setSpecies(speciesWall);

      // std::cout << "--- 2 ---" << std::endl;

      // creates fictious particles at every grid node
      // if they are inside the agglomerate perimeter then save them in the .grid file
      // if they are inside a particle label with a "1", with a "0" otherwise
      for(int i = 0; i < gridLength; i++)
      {
         gridPoint.X = centerOfMass.X - meanGranuleRadius + gridResolutionLength*i;

         for(int j = 0; j < gridLength; j++)
         {
            gridPoint.Y = centerOfMass.Y - meanGranuleRadius + gridResolutionLength*j;

            for(int k = 0; k < gridLength; k++)
            {
               gridPoint.Z = centerOfMass.Z - meanGranuleRadius + gridResolutionLength*k;

               // std::cout << i << "\t" << j << "\t" << k << "\t" << gridPoint << "\t" << (gridPoint - centerOfMass).getLength() << "\t" << (gridPoint - centerOfMass).getLength()/meanGranuleRadius << "\t" << nPointsInsideAgglomerateBoundary << std::endl;

               p0.setPosition(gridPoint);

               // checks if inside the agglomerate boundary
               if ((gridPoint - centerOfMass).getLength() < meanGranuleRadius)
               {
                  nPointsInsideAgglomerateBoundary++;

                  if (checkParticleForInteraction(p0)) // no collision -> the counter goes to the void fraction
                  {
                     gridFile << gridPoint << "\t" << 0 << std::endl;
                  }
                  else // collision -> the counter goes to the mass fraction
                  {
                     nPointsInsideComponents++;
                     gridFile << gridPoint << "\t" << 1 << std::endl;
                  }
               }
            }
         }
      }

      massFraction = nPointsInsideComponents/nPointsInsideAgglomerateBoundary;

      gridFile << "n_points_inside_boundary: " << nPointsInsideAgglomerateBoundary << "\t n_points_inside_components: " << nPointsInsideComponents << "\t mass_fraction: " << massFraction << std::endl;
      gridFile.close();
   }

   void activateGravity()
   {
      setGravity(Vec3D(0.00,0.00,-9.81*(getTime() - t0)/gravityAndForceTuningDuration));
   }

   void makeBaseAndHousing()
   {
      // clearing every wall
      wallHandler.clear();

      baseHeight = particleHandler.getLowestPositionComponentParticle(2) -> getPosition().Z - 1.1*(1.0 + sizeDispersityParticle)*radiusParticle;
      base.setSpecies(speciesWall);
      base.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,baseHeight));
      basePointer = wallHandler.copyAndAddObject(base);

      cylinder.setSpecies(speciesWall);
      cylinder.setPosition(Vec3D(0.0,0.0,0.0));
      cylinder.setOrientation(Vec3D(0.0,0.0,1.0));
      cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(0.5*boxSize,0.0,0.0));
      wallHandler.copyAndAddObject(cylinder);
   }

   void makeCompressionPiston()
   {
      pistonStepDisplacement = (1.0 - sizeDispersityParticle)*radiusParticle/10000.0;
      pistonHeight = particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z + 1.5*(1.0 + sizeDispersityParticle)*radiusParticle;

      initialClusterSize = 2.0*meanGranuleRadius;

      piston.setSpecies(speciesWall);
      piston.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(piston);
   }

   void loadPiston()
   {
      pistonHeight -= pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      // pistonPointer -> setVelocity(Vec3D(0.0,0.0,-pistonVelocity));
   }

   void unloadPiston()
   {
      pistonHeight += pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      // pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));
   }


   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      if (stage < 5)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(4) << getTime() << ", tmax = " << getTimeMax() <<
         ", E_ratio = " << getKineticEnergy()/getElasticEnergy() <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfBindingInteractions << ", rMean = " << meanGranuleRadius << ", vF = " << volumeRatio <<
         ", Force Modulus = " << forceModulus << ", g_z = " << getGravity().Z << ", dMean = " << meanTotalRelativeOverlap << ", dMax = " << maxTotalRelativeOverlap <<
         std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(4) << getTime() << ", tmax = " << getTimeMax() <<
         ", E_ratio = " << getKineticEnergy()/getElasticEnergy() <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfBindingInteractions <<
         ", Fpiston = " << pistonForce << ", Hpiston = " << pistonHeight << ", dMean = " << meanTotalRelativeOverlap << ", dMax = " << maxTotalRelativeOverlap <<
         std::endl;
      }

      std::cout.flush();
   }


   // VARIABLES -------------------------------------------------------
   // particles
   double radiusParticle;
   double sizeDispersityParticle;
   double densityParticle;
   double volumeParticle;
   double massParticle;
   double totalParticleVolume;
   int nParticlesInserted;

   // granules
   int totalNumberOfParticles;
   double granuleRadiusPreFormation;
   double meanGranuleRadius;
   Vec3D centerOfMass;

   // geometry
   double boxSize;
   double granuleSizeToLatticeRatio;
   InfiniteWall wall;

   // friction
   double particleWallSlidingFrictionCoeff;
   double particleWallRollingFrictionCoeff;
   double particleWallTorsionFrictionCoeff;
   double particleParticleSlidingFrictionCoeff;
   double particleParticleRollingFrictionCoeff;
   double particleParticleTorsionFrictionCoeff;

   // collision
   double kpWall, kpParticle;
   double keParticle;
   double kcParticleIntergranule;
   double kcParticleLoose;
   double phiParticle;
   double particleWallRestitutionCoeff;
   double particleParticleRestitutionCoeff;
   double densityWall;

   // central force
   double forceModulus;
   double maximumForceModulus;
   double finalForceModulus;
   double forceDampingModulus;
   double forceTuningInterval;
   double velocityDampingModulus;
   double velocityDampingInterval;

   // contact related
   int numberOfBindingInteractions;
   std::vector< std::vector<int> > coordinationMatrix;
   std::vector< std::vector<int> > coordinationMatrixPreviousCheck;

   // species
   LinearPlasticViscoelasticFrictionSpecies *speciesParticle, *speciesWall;
   std::vector<LinearPlasticViscoelasticFrictionSpecies*> particlePureSpeciesVector;
   std::vector< std::vector<LinearPlasticViscoelasticFrictionMixedSpecies*> > mixedSpeciesMatrix;

   // output
   double gdatOutputTimeInterval;
   std::ofstream gdatFile;
   std::ofstream gridFile;

   // data analysis
   double meanCoordinationNumber;
   double maxTotalRelativeOverlap;
   double meanTotalRelativeOverlap;
   double volumeRatio;
   double averageForceOnParticle;
   double massFraction;
   double gridLength;
   double gridResolutionLength;
   double fictiousGridPointRadiusRatio;

   // compression
   InfiniteWall piston;
   InfiniteWall* pistonPointer;
   InfiniteWall base;
   InfiniteWall* basePointer;
   AxisymmetricIntersectionOfWalls cylinder;
   double baseHeight;
   double maximumCompressiveForce;
   double minimumRelativeHeightCompression;
   double pistonForce;
   double pistonStepDisplacement;
   double pistonHeight;
   double initialClusterSize;

   // global
   int stage;
   double t0;
   double energyRatioTolerance;
   double gravityAndForceTuningDuration;
};


int main(int argc, char *argv[])
{
   // SPHERE AND PSD
   // possible eventual shapes: SPHERE, CYLINDER, ELLIPSOID, PROLATE, OBLATE
   std::string granuleShape = "SPHERE";
   std::string psd = "UNIFORM";

   // TIME STEP
   double timeStep = 2.0e-6;

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;
   double sizeDispersionParticles = 0.1;
   double densityParticles = 1500.0;
   double densityWall = 5000.0;

   // INTERACTION PARAMETERS
   double plasticStiffness = 750.0;
   double elasticStiffness = 1000.0;
   double intergranuleCohesionStiffness = 100.0;
   // double intergranuleCohesionStiffness[7] = {1.0, 50.0, 100.0, 150.0, 200.0, 250.0, 500.0};
   // double intergranuleCohesionStiffness[6] = {250.0, 500.0, 750.0, 1000.0, 1250.0, 1500.0};
   double looseParticleCohesionStiffness = 0.0;
   double phiParticle = 0.1;
   // double phiParticle[3] = {0.1, 0.2, 0.5};
   double particleParticleRestitutionCoefficient = 0.5;
   double particleWallRestitutionCoefficient = 0.7;

   double elasticStiffnessWalls = 1000.0;

   // FRICTION PARAMETERS
   double muSlidingParticleParticle = 0.5;
   double muRollingParticleParticle = 0.3;
   double muSlidingParticleWall = 0.5;
   double muRollingParticleWall = 0.3;

   // FORCE AND DAMPING PARAMETERS
   double maximumForceMultiplicationFactor = 200.0;
   double finalForceMultiplicationFactor = 0.0001;
   double forceTuningInterval = timeStep;
   double velocityDampingModulus = 0.9;
   double velocityDampingInterval = 0.01;

   // GRANULES PARAMETERS
   int totalNumberOfParticles = 200;
   // int totalNumberOfParticles[7] = {10, 20, 50, 100, 200, 500, 1000};
   double granuleSizeToLatticeRatio = 0.5;

   // INTERNAL STRUCTURE ANALYSIS PARAMETERS
   int gridLength = 300;
   double fictiousGridParticleRadiusRatio = 1.0e-5;

   // COMPRESSION PARAMETERS
   double maximumCompressiveForce = 1.0;
   double minimumRelativeHeightCompression = 5.0/7.0;

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 1.0e-4;
   double gravityAndForceTuningDuration = 0.2;

   // for (int j = 0; j < 3; j++)
   // {
      // for (int i = 0; i < 6; i++)
      // {
         // INITIALIZATION
         GranulesUniaxialCompression problem;

         problem.setTimeStep(timeStep);
         problem.setTimeMax(10.0);
         problem.setGravity(Vec3D(0.00,0.00,0.00));
         problem.setSystemDimensions(3);
         problem.setCdatOutputTimeInterval(0.001);
         problem.setEnergyRatioTolerance(energyRatioTolerance);
         problem.setSaveCount(0.001/problem.getTimeStep());

         problem.setGranulesProperties(totalNumberOfParticles, granuleSizeToLatticeRatio);
         problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
         problem.setWallDensity(densityWall);
         problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
         problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
         problem.setParticlePlasticProperties(plasticStiffness, elasticStiffness, intergranuleCohesionStiffness, looseParticleCohesionStiffness, phiParticle);
         problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
         problem.setWallStiffnessAndRestitutionCoefficients(elasticStiffnessWalls, particleWallRestitutionCoefficient);
         problem.setForceAndVelocityProperties(maximumForceMultiplicationFactor*particleRadius*elasticStiffness, finalForceMultiplicationFactor*particleRadius*elasticStiffness, forceTuningInterval, velocityDampingModulus, velocityDampingInterval);
         problem.setGravityAndForceTuningDuration(gravityAndForceTuningDuration);
         problem.setInternalStructureAnalysisParameters(gridLength, fictiousGridParticleRadiusRatio);
         problem.setCompressionParameters(maximumCompressiveForce, minimumRelativeHeightCompression);

         // NAME SETTING
         std::ostringstream name;
         name.str("");
         name.clear();
         std::cout.unsetf(std::ios::floatfield);
         // name << "GranuleUniaxialCompression___SPHERE___nP_" << totalNumberOfParticles << "_" << psd << "_rDisp_" << sizeDispersionParticles << "_kP_" << plasticStiffness << "_kE_" << elasticStiffness
         // << "_kCi_" << intergranuleCohesionStiffness << "_kCe_" << looseParticleCohesionStiffness << "_phi_" << phiParticle
         // << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle << "_fRatio_" << maximumForceMultiplicationFactor << "_Fmax_" << maximumCompressiveForce;
         name << "GranuleUniaxialCompression___SPHERE___nP_" << totalNumberOfParticles << "_" << psd << "_rDisp_" << sizeDispersionParticles << "_kP_" << plasticStiffness << "_kE_" << elasticStiffness
         << "_kCi_" << intergranuleCohesionStiffness << "_kCe_" << looseParticleCohesionStiffness << "_phi_" << phiParticle
         << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle << "_fRatio_" << maximumForceMultiplicationFactor << "_Hmin_" << minimumRelativeHeightCompression;
         problem.setName(name.str());

         problem.solve();
      // }
   // }

   return 0;
}
