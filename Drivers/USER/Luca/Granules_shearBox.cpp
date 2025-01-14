/*
 *** GRANULES SHEAR CELL ***
 Granules are created in a lattice, with periodic boundaries along X and Y direction.
 Particle-made flat walls are loaded at zMin and zMax and moved together to reach a prescribed packing fraction OR pressure OR distance.
 The walls are then moved in opposite directions to shear the granules in between.

 LAST UPDATE: 9.8.18
 */

#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <fstream>

/*
 ToDo:
- patch particles to walls to improve agglomerate convection
 */

class GranuleShearBox : public Mercury3D
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

      std::cout << "CREATING BOUNDARIES..." << std::endl;
      makeBoundaries();

      std::cout << "CREATING LATTICE GRID" << std::endl;
      setLatticeGrid();

      std::cout << std::endl << "PARTICLE INSERTION" << std::endl;
      insertParticles();

      std::cout << "CREATING CDAT FILE" << std::endl << std::endl;
      makeCdatFile();

      std::cout << "ACTIVATING CENTRAL FORCES" << std::endl;

      stage++;
   }

   void actionsOnRestart() override
   {

   }

   void actionsAfterTimeStep() override
   {
      // DATA ANALYSIS
      if (stage >= 3) makeDataAnalysis();

      // UPDATE OF COORDINATION MATRIX AND COHESION MATRIX AND OUTPUT TO CDAT
      if (stage >= 3 && fmod(getTime(), 0.001) < getTimeStep()) refreshCoordinationMatrixAndCohesionMatrix();
      if (stage >= 3 && fmod(getTime(), 0.001) < getTimeStep()) writeToCdatFile();

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
            std::cout << "CREATING COORDINATION MATRIX AND REFRESHING COORDINATION MATRIX" << std::endl << std::endl;
            setupCoordinationMatrix();
            refreshCoordinationMatrixAndCohesionMatrix();

            std::cout << "COMPUTING AVERAGE GRANULE RADIUS" << std::endl << std::endl;
            computeAverageGranuleRadius();

            std::cout << "TUNING GRAVITY" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // GRAVITY TUNING
      if (stage == 3)
      {
         tuneGravity();

         if (getTime() - t0 > gravityAndForceTuningDuration)
         {
            std::cout << "SETTLING" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // ENERGY DISSIPATION
      if (stage == 4)
      {
         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance || getTime() - t0 > gravityAndForceTuningDuration)
         {
            std::cout << "ENERGY DISSIPATED. STARTING COMPRESSION" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // COMPRESSION STAGE
      if (stage == 5)
      {
         wiggleWalls();

         if (isPressureDriven)
         {
            if (topWallStress.Z < targetPressure) moveTopWall();
            else
            {
               std::cout << "STARTING SHEAR" << std::endl << std::endl;
               t0 = getTime();
               stopWiggling();
               shearWalls();
               setTimeMax(t0 + shearDuration);
               stage++;
            }
         }
         else if (isPackingFractionDriven)
         {
            computeGranularPackingFraction();

            if (granularPackingFraction < targetPackingFraction)
            {
               moveTopWall();
            }
            else
            {
               std::cout << "STARTING SHEAR" << std::endl << std::endl;
               t0 = getTime();
               stopWiggling();
               shearWalls();
               setTimeMax(t0 + shearDuration);
               stage++;
            }
         }
         else
         {
            std::cout << " SOMETHING WENT WRONG WHILE DECLARING THE SHEARING CONDITIONS. EXITING..." << std::endl;
            exit(1);
         }
      }

      // SHEAR STAGE
      if (stage == 6)
      {
         if (isPressureDriven) mantainTargetPressure();
      }
   }

   void actionsAfterSolve() override
   {
      delete [] boxCentreCoordinates;
      cdatFile.close();
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

   void setGranulesProperties(int nP, double nu)
   {
      nParticlesPerGranule = nP;
      assumedInterGranulePackingFraction = nu;
   }

   void setLatticeDimensions(double nX, double nY, double nZ, double sizeRatio)
   {
      nXlayers = nX;
      nYlayers = nY;
      nZlayers = nZ;
      granuleSizeToLatticeRatio = sizeRatio;

      numberOfGranules = nXlayers*nYlayers*nZlayers;
      totalNumberOfParticles = numberOfGranules*nParticlesPerGranule;

      std::cout << "Number of granules needed " << numberOfGranules << std::endl;
      std::cout << "Total number of particles needed " << totalNumberOfParticles << std::endl << std::endl;
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

   // void setShearAndCompressionParameters(double sv, double sd, double pTarget)
   // {
   //    shearingVelocity = sv;
   //    shearDuration = sd;
   //    targetPressure = pTarget;
   //
   //    compressionVelocity = -0.001*radiusParticle/getTimeStep();
   // }

   void setShearAndPressureParameters(double sv, double sd, double pTarget)
   {
      shearingVelocity = sv;
      shearDuration = sd;
      targetPressure = pTarget;

      isPressureDriven = true;
      isPackingFractionDriven = false;

      compressionVelocity = -0.001*radiusParticle/getTimeStep();
   }

   void setShearAndPackingFractionParameters(double sv, double sd, double pfTarget)
   {
      shearingVelocity = sv;
      shearDuration = sd;
      targetPackingFraction = pfTarget;

      isPressureDriven = false;
      isPackingFractionDriven = true;

      compressionVelocity = -0.001*radiusParticle/getTimeStep();
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
      granuleRadius = radiusParticle*pow(nParticlesPerGranule/assumedInterGranulePackingFraction,1.0/3.0);
      boxSize = 2.0*granuleRadius/granuleSizeToLatticeRatio;

      std::cout << "Granule average radius " << granuleRadius << std::endl;
      std::cout << "Cubic lattice size " << boxSize << std::endl << std::endl;
   }

   void setDomainLimits()
   {
      setXMin(0.0);
      setYMin(-0.5*nYlayers*boxSize);
      setZMin(0.0);

      setXMax(nXlayers*boxSize);
      setYMax(0.5*nYlayers*boxSize);
      setZMax(nZlayers*boxSize);
   }

   void makeBoundaries()
   {
      wallHandler.clear();

      // z walls
      bottomWall.setSpecies(speciesWall);
      bottomWall.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,getZMin()));
      bottomWallPointer = wallHandler.copyAndAddObject(bottomWall);

      topWall.setSpecies(speciesWall);
      topWall.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,getZMax()));
      topWallPointer = wallHandler.copyAndAddObject(topWall);

      topWallHeight = getZMax();

      // x boundary
      xBoundary.set(Vec3D(1.0,0.0,0.0), getXMin(), getXMax());
      boundaryHandler.copyAndAddObject(xBoundary);

      // y boundary
      yBoundary.set(Vec3D(0.0,1.0,0.0), getYMin(), getYMax());
      boundaryHandler.copyAndAddObject(yBoundary);
   }

   void setLatticeGrid()
   {
      boxCentreCoordinates = new Vec3D[nXlayers*nZlayers*nYlayers];

      double offsetX = 0.5*boxSize;
      double offsetY = -0.5*(nYlayers - 1.0)*boxSize;
      double offsetZ = 0.5*boxSize;

      for (int i = 0; i < nYlayers; i++) // loop along Y
      {
         for (int j = 0; j < nZlayers; j++) // loop along Z
         {
            for (int k = 0; k < nXlayers; k++) // loop along X
            {
               boxCentreCoordinates[k + j*nXlayers + i*nXlayers*nZlayers].X = offsetX + k*boxSize;
               boxCentreCoordinates[k + j*nXlayers + i*nXlayers*nZlayers].Y = offsetY + i*boxSize;
               boxCentreCoordinates[k + j*nXlayers + i*nXlayers*nZlayers].Z = offsetZ + j*boxSize;
            }
         }
      }
   }

   bool particleInsertionSuccessful(int nGranule)
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
         // particlePosition.X = boxCentreCoordinates[nGranule].X + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         // particlePosition.Y = boxCentreCoordinates[nGranule].Y + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         // particlePosition.Z = boxCentreCoordinates[nGranule].Z + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);

         rad = random.getRandomNumber(p0.getRadius(),0.5*boxSize - 1.01*(p0.getRadius()));
         theta = constants::pi*random.getRandomNumber(-1.0,1.0);
         phi = 0.5*constants::pi*random.getRandomNumber(-1.0,1.0);
         particlePosition.X = boxCentreCoordinates[nGranule].X + rad*sin(theta)*cos(phi);
         particlePosition.Y = boxCentreCoordinates[nGranule].Y + rad*sin(theta)*sin(phi);
         particlePosition.Z = boxCentreCoordinates[nGranule].Z + rad*cos(theta);

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
      for (int i = 0; i < numberOfGranules; i++)
      {
         nParticlesInserted = 0;

         while(nParticlesInserted < nParticlesPerGranule)
         {
            if (particleInsertionSuccessful(i)) nParticlesInserted++;
            else t0 = getTime();
         }

         std::cout << "Inserted granule n. " << i + 1 << "/" << numberOfGranules << std::endl;
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
   }

   void applyCentralForce()
   {
      Vec3D distanceFromForceCenter;
      int xLatticePosition;
      int yLatticePosition;
      int zLatticePosition;

      for (int j=0; j<numberOfGranules; j++) // count the particles belonging to every single granule and point them towards their box centre
      {
         for (int i=0; i<nParticlesPerGranule; i++)
         {
            distanceFromForceCenter = particleHandler.getObject(i + nParticlesPerGranule*j) -> getPosition() - boxCentreCoordinates[j];

            //particleHandler.getObject(i) -> addForce(-forceScaleModulus*Vec3D(0.2*distanceFromForceCenter.X,0.2*distanceFromForceCenter.Y,distanceFromForceCenter.Z));
            particleHandler.getObject(i + nParticlesPerGranule*j) -> addForce(-forceModulus*distanceFromForceCenter);
            //particleHandler.getObject(i) -> addForce(-forceScaleModulus*distanceFromForceCenter*(particleHandler.getObject(i) -> getRadius())/(0.5*(radiusParticle + radiusSmall)));
            //particleHandler.getObject(i) -> addForce(-1.0*distanceFromForceCenter/(distanceFromForceCenter.getLength()*distanceFromForceCenter.getLength())/10000);
         }
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

   void tuneGravity()
   {
      if (getTime() - t0 < gravityAndForceTuningDuration) setGravity(Vec3D(0.0,0.0,-9.81*(getTime() - t0)/gravityAndForceTuningDuration));
   }

   void makeDataAnalysis()
   {
      meanCoordinationNumber = 0.0;
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
      }
      meanCoordinationNumber /= particleHandler.getNumberOfObjects();

      maxTotalRelativeOverlap = 0.0;
      meanTotalRelativeOverlap = 0.0;
      maxTopWallRelativeOverlap = 0.0;
      meanTopWallRelativeOverlap = 0.0;
      maxBottomWallRelativeOverlap = 0.0;
      meanBottomWallRelativeOverlap = 0.0;

      // topWallForce.setZero();
      topWallStress.setZero();
      // bottomWallForce.setZero();
      bottomWallStress.setZero();

      int totalInteractionCounter = 0;
      int bottomWallInteractionCounter = 0;
      int topWallInteractionCounter = 0;

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

         // bottom wall interactions
         if ((*i) -> getI() -> getIndex() == bottomWallPointer -> getIndex())
         {
            bottomWallStress -= ((*i) -> getForce());

            meanBottomWallRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxBottomWallRelativeOverlap)
            {
               maxBottomWallRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            }

            bottomWallInteractionCounter++;
         }

         // top wall interactions
         if ((*i) -> getI() -> getIndex() == topWallPointer -> getIndex())
         {
            topWallStress -= ((*i) -> getForce());

            meanTopWallRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTopWallRelativeOverlap)
            {
               maxTopWallRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            }

            topWallInteractionCounter++;
         }

         totalInteractionCounter++;
      }

      meanTotalRelativeOverlap /= totalInteractionCounter;
      meanBottomWallRelativeOverlap /= bottomWallInteractionCounter;
      meanTopWallRelativeOverlap /= topWallInteractionCounter;
      bottomWallStress /= (nXlayers + nYlayers)*pow(boxSize,2.0);
      topWallStress /= (nXlayers + nYlayers)*pow(boxSize,2.0);
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

   // refreshes the matrix containing the interaction information and updates the cohesion coefficients accordingly
   void refreshCoordinationMatrixAndCohesionMatrix()
   {
      // resets the number of binding interactions (the bounds after the granule formation)
      numberOfBindingInteractions = 0;

      // // prints the matrix before
      // std::ofstream matrixOutput2;
      // matrixOutput2.open("Granules_COHESIONMATRIX_BEFORE", std::ios::out);
      // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      // {
      //    for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) matrixOutput2 << mixedSpeciesMatrix[i][j] -> getCohesionStiffness() << "\t";
      //    matrixOutput2 << std::endl;
      // }
      // matrixOutput2.close();

      // resets the coordination matrix to 0
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) coordinationMatrix[i][j] = 0;
      }

      // re-allocates the contacts in the coordination matrix
      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         // the contact is preserved if it was in place at the previous check
         if (coordinationMatrixPreviousCheck[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()])
         {
            coordinationMatrix[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
            coordinationMatrix[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;

            // updates the number of bounds
            numberOfBindingInteractions++;
         }
      }

      // sets the old coordination matrix equal to the new one
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         for (int j = i + 1; j < particleHandler.getNumberOfObjects(); j++)
         {
            coordinationMatrixPreviousCheck[i][j] = coordinationMatrix[i][j];
            coordinationMatrixPreviousCheck[j][i] = coordinationMatrix[i][j];
         }
      }

      // updates the cohesion interactions
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         for (int j = i + 1; j < particleHandler.getNumberOfObjects(); j++)
         {
            if (coordinationMatrix[i][j]) mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleIntergranule);
            else mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleLoose);
         }
      }

      // // print the cohesive interaction matrix
      // std::ofstream matrixOutput;
      // matrixOutput.open("Granules_COHESIONMATRIX_AFTER", std::ios::out);
      // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      // {
      //    for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) matrixOutput << mixedSpeciesMatrix[i][j] -> getCohesionStiffness() << "\t";
      //    matrixOutput << std::endl;
      // }
      // matrixOutput.close();
   }

   // creates the data output file and writes the first row
   void makeCdatFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "time \t Eel \t Ekin/Eel \t h \t top_stress_X \t top_stress_Y \t top_stress_Z \t bot_stress_X \t bot_stress_Y \t bot_stress_Z \t" <<
      " coord_number \t n_bounds \t dTotMean \t dTotMax \t dTopMean \t dTopMax \t dBotMean \t dBotMax" << std::endl;
   }

   // writes the compression data to the output file
   void writeToCdatFile()
   {
      cdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      topWallHeight << "   " <<
      topWallStress.X << "   " <<
      topWallStress.Y << "   " <<
      topWallStress.Z << "   " <<
      bottomWallStress.X << "   " <<
      bottomWallStress.Y << "   " <<
      bottomWallStress.Z << "   " <<
      meanCoordinationNumber << "   " <<
      numberOfBindingInteractions << "   " <<
      meanTotalRelativeOverlap << "   " <<
      maxTotalRelativeOverlap << "   " <<
      meanTopWallRelativeOverlap << "   " <<
      maxTopWallRelativeOverlap << "   " <<
      meanBottomWallRelativeOverlap << "   " <<
      maxBottomWallRelativeOverlap << "   " <<
      std::endl;
   }

   void moveTopWall()
   {
      topWallHeight += compressionVelocity*getTimeStep();
      topWallPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,topWallHeight));
      // topWallPointer -> setVelocity(Vec3D(0.0,0.0,compressionVelocity));
   }

   void shearWalls()
   {
      topWallPointer -> setVelocity(Vec3D(shearingVelocity,0.0,0.0));
      // bottomWallPointer -> setVelocity(Vec3D(-0.5*shearingVelocity,0.0,0.0));
   }

   void mantainTargetPressure()
   {
      if ((topWallStress.Z < targetPressure && compressionVelocity > 0.0) || (topWallStress.Z > targetPressure && compressionVelocity < 0.0)) compressionVelocity *= -1.0;
      if (fabs(topWallStress.Z - targetPressure)/targetPressure > 1.0e-4)
      {
         topWallHeight += compressionVelocity*getTimeStep();
         topWallPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,topWallHeight));
      }
   }

   void computeAverageGranuleRadius()
   {
      Vec3D distanceFromForceCenter;
      Vec3D localMin, localMax;
      meanGranuleRadius = 0.0;

      for (int j=0; j<numberOfGranules; j++)
      {
         localMin.setZero();
         localMax.setZero();

         for (int i=0; i<nParticlesPerGranule; i++)
         {
            distanceFromForceCenter = particleHandler.getObject(i + nParticlesPerGranule*j) -> getPosition() - boxCentreCoordinates[j];

            if (distanceFromForceCenter.X > localMax.X) localMax.X = distanceFromForceCenter.X;
            if (distanceFromForceCenter.X < localMin.X) localMin.X = distanceFromForceCenter.X;
            if (distanceFromForceCenter.Y > localMax.Y) localMax.Y = distanceFromForceCenter.Y;
            if (distanceFromForceCenter.Y < localMin.Y) localMin.Y = distanceFromForceCenter.Y;
            if (distanceFromForceCenter.Z > localMax.Z) localMax.Z = distanceFromForceCenter.Z;
            if (distanceFromForceCenter.Z < localMin.Z) localMin.Z = distanceFromForceCenter.Z;
         }

         // meanGranuleRadius += (localMax - localMin).getLength();
         // meanGranuleRadius += localMin.getLength();
         // meanGranuleRadius += localMax.getLength();
         meanGranuleRadius += localMax.X - localMin.X + localMax.Y - localMin.Y + localMax.Z - localMin.Z;
      }

      meanGranuleRadius /= 6.0*numberOfGranules;
      std::cout << "The average granule radius is: " << std::setprecision(3) << std::left << std::setw(6) << meanGranuleRadius << std::endl;
   }

   void computeGranularPackingFraction()
   {
      granularPackingFraction =  (4.0*constants::pi*numberOfGranules*pow(meanGranuleRadius,3.0))/(3.0*(getXMax() - getXMin())*(getYMax() - getYMin())*topWallHeight);
   }

   void wiggleWalls()
   {
      double wigglingOmega = 2.0*constants::pi/0.001;
      double wigglingVelocity = 0.05*shearingVelocity;
      topWallPointer -> setVelocity(Vec3D(wigglingVelocity*sin(wigglingOmega*(getTime() - t0)),0.0,0.0));
      bottomWallPointer -> setVelocity(Vec3D(-wigglingVelocity*sin(wigglingOmega*(getTime() - t0)),0.0,0.0));
   }

   void stopWiggling()
   {
      topWallPointer -> setVelocity(Vec3D(0.0,0.0,0.0));
      bottomWallPointer -> setVelocity(Vec3D(0.0,0.0,0.0));
   }

   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      if (stage < 5)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", E_ratio = " << std::setprecision(4) << std::left << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
         ", Force Modulus = " << forceModulus << ", g_y = " << getGravity().Z << ", STAGE " << std::setprecision(0) << stage
         << std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", Eel = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Ekin/Eel = " << getKineticEnergy()/getElasticEnergy() <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfBindingInteractions << std::endl <<
         ", h = " << topWallHeight << ", Ptop = " << std::setprecision(3) << std::left << std::setw(6) << topWallStress.Z << ", Pbot = " << bottomWallStress.Z <<
         ", xShearTop = " << topWallStress.X << ", xShearBot = " << bottomWallStress.X << std::endl << std::endl;
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
   int nParticlesInserted;

   // granules
   int numberOfGranules;
   int nParticlesPerGranule;
   int totalNumberOfParticles;
   double granuleRadius;
   double assumedInterGranulePackingFraction;
   double granularPackingFraction;
   double meanGranuleRadius;
   Vec3D *boxCentreCoordinates;

   // geometry
   double boxSize;
   double granuleSizeToLatticeRatio;
   int nXlayers;
   int nYlayers;
   int nZlayers;
   InfiniteWall topWall, bottomWall;
   InfiniteWall *topWallPointer, *bottomWallPointer;

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

   // shear
   double targetPressure;
   double targetPackingFraction;
   double topWallHeight;
   double compressionVelocity;
   double shearingVelocity;
   double shearDuration;

   // boundaries
   PeriodicBoundary xBoundary, yBoundary;

   // output
   double cdatOutputTimeInterval;
   std::ofstream cdatFile;

   // data analysis
   double meanCoordinationNumber;
   // Vec3D topWallForce;
   Vec3D topWallStress;
   // Vec3D bottomWallForce;
   Vec3D bottomWallStress;
   double maxTotalRelativeOverlap;
   double meanTotalRelativeOverlap;
   double maxTopWallRelativeOverlap;
   double meanTopWallRelativeOverlap;
   double maxBottomWallRelativeOverlap;
   double meanBottomWallRelativeOverlap;

   // global
   int stage;
   double t0;
   double energyRatioTolerance;
   double gravityAndForceTuningDuration;
   bool isPressureDriven = false;
   bool isPackingFractionDriven = false;
};


int main(int argc, char *argv[])
{
   // COMPRESSION DRIVING MODE
   std::string compressionDrivingMode = "PRESSURE";

   // TIME STEP
   double timeStep = 5.0e-6;

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;
   double sizeDispersionParticles = 0.1;
   double densityParticles = 1500.0;
   double densityWall = 5000.0;

   // INTERACTION PARAMETERS
   double plasticStiffness = 500.0;
   double elasticStiffness = 1000.0;
   double intergranuleCohesionStiffness = 500.0;
   double looseParticleCohesionStiffness = 0.0;
   double phiParticle = 0.1;
   double particleParticleRestitutionCoefficient = 0.5;
   double particleWallRestitutionCoefficient = 0.7;

   double elasticStiffnessWalls = 3000.0;

   // FRICTION PARAMETERS
   double muSlidingParticleParticle = 0.5;
   double muRollingParticleParticle = 0.3;
   double muSlidingParticleWall = 0.8;
   double muRollingParticleWall = 0.5;

   // FORCE AND DAMPING PARAMETERS
   double maximumForceModulus = 1000.0*particleRadius*elasticStiffness;
   double finalForceModulus = 0.001*particleRadius*elasticStiffness;
   double forceTuningInterval = 0.002;
   double velocityDampingModulus = 0.9;
   double velocityDampingInterval = 0.01;

   // GRANULES PARAMETERS
   int nParticlesPerGranule = 100;
   int nXWells = 5;
   int nYWells = 4;
   int nZWells = 5;
   double assumedInterGranulePackingFraction = 0.80;
   double granuleSizeToLatticeRatio = 0.75;

   // SHEAR PARAMETERS
   double preShearPressure = 10000.0;
   double preShearPackingFraction = 0.6;
   double wallVelocity = 0.1;
   double shearDuration = 1.0;

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 1.0e-4;
   double gravityAndForceTuningDuration = 0.2;

   // INITIALIZATION
   GranuleShearBox problem;

   problem.setTimeStep(timeStep);
   problem.setTimeMax(100.0);
   problem.setGravity(Vec3D(0.00,0.00,0.00));
   problem.setSystemDimensions(3);
   problem.setCdatOutputTimeInterval(0.001);
   problem.setEnergyRatioTolerance(energyRatioTolerance);
   problem.setSaveCount(0.001/problem.getTimeStep());

   problem.setGranulesProperties(nParticlesPerGranule, assumedInterGranulePackingFraction);
   problem.setLatticeDimensions(nXWells, nYWells, nZWells, granuleSizeToLatticeRatio);

   problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
   problem.setWallDensity(densityWall);

   problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
   problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
   problem.setParticlePlasticProperties(plasticStiffness, elasticStiffness, intergranuleCohesionStiffness, looseParticleCohesionStiffness, phiParticle);
   problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
   problem.setWallStiffnessAndRestitutionCoefficients(elasticStiffnessWalls, particleWallRestitutionCoefficient);

   problem.setForceAndVelocityProperties(maximumForceModulus, finalForceModulus, forceTuningInterval, velocityDampingModulus, velocityDampingInterval);
   problem.setGravityAndForceTuningDuration(gravityAndForceTuningDuration);
   // problem.setShearAndCompressionParameters(wallVelocity, shearDuration, preShearPressure);

   if (compressionDrivingMode == "PACKINGFRACTION")
   {
      problem.setShearAndPackingFractionParameters(wallVelocity, shearDuration, preShearPackingFraction);

      // NAME SETTING
      std::ostringstream name;
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      name << "GranuleShearBox_nP_" << nParticlesPerGranule << "_" << nXWells << "_" << nYWells << "_" << nZWells << "_" << "_kP_" << plasticStiffness << "_kE_" << elasticStiffness
      << "_kCi_" << intergranuleCohesionStiffness << "_kCe_" << looseParticleCohesionStiffness << "_phi_" << phiParticle << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle
      << "_pf_" << preShearPackingFraction << "_v_" << wallVelocity;
      problem.setName(name.str());
   }

   if (compressionDrivingMode == "PRESSURE")
   {
      problem.setShearAndPressureParameters(wallVelocity, shearDuration, preShearPressure);

      // NAME SETTING
      std::ostringstream name;
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      name << "GranuleShearBox_nP_" << nParticlesPerGranule << "_" << nXWells << "_" << nYWells << "_" << nZWells << "_" << "_kP_" << plasticStiffness << "_kE_" << elasticStiffness
      << "_kCi_" << intergranuleCohesionStiffness << "_kCe_" << looseParticleCohesionStiffness << "_phi_" << phiParticle << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle
      << "_P_" << preShearPressure << "_v_" << wallVelocity;
      problem.setName(name.str());
   }

   problem.solve();

   return 0;
}

// XBALLS ARGUMENTS TO ADD
// -h 800 -p 1 -s 8 -v0 -3dturn 3
