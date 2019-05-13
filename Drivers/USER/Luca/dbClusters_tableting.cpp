/*
 *** CLUSTERS AOR TEST ***
 Clusters are created in a lattice, with periodic boundaries along Y direction.
 When settled, the x+ wall is removed and the clusters let to flow out of the box.
 When settled the angle of repose is computed.

 LAST UPDATE: 21.9.18
 */

#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <fstream>

/*
 ToDo:
- patch particles to walls to improve agglomerate convection
 */

class ClusterTableting : public Mercury3D
{
private:
   void setupInitialConditions() override
   {
      stage = 0;
      nParticlesInserted = 0;
      numberOfBindingInteractions = 0;
      forceModulus = 0.0;
      pistonHeight = 0.0;
      pistonForce = 0.0;
      random.randomise();
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

      std::cout << "CREATING .cdat FILE" << std::endl;
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
      if (stage >= 1 && fmod(getTime(), 0.001) < getTimeStep()) writeToCdatFile();

      // UPDATE OF COORDINATION MATRIX AND COHESION MATRIX
      if (stage >= 3 && fmod(getTime(), 0.001) < getTimeStep()) refreshCoordinationMatrixAndCohesionMatrix();

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
         if (getTime() - t0 > 3.0*gravityAndForceTuningDuration) // getKineticEnergy()/getElasticEnergy() < energyRatioTolerance ||
         {
            std::cout << "ENERGY DISSIPATED. CREATING PISTON AND STARTING COMPRESSION" << std::endl << std::endl;
            makeCompressionPiston();

            std::cout << "Initial bed height: " << initialBedHeight << std::endl;
            std::cout << "Piston height: " << pistonHeight << std::endl;
            std::cout << "Final target piston height: " << initialBedHeight*minimumRelativeHeightCompression << std::endl;

            t0 = getTime();
            stage++;
         }
      }

      // COMPRESSION STAGE
      if (stage == 5)
      {
         // if (pistonForce < maximumCompressiveForce) // uncomment this for force-driven compression
         if (pistonHeight > initialBedHeight*minimumRelativeHeightCompression) // uncomment this for height-driven compression
         {
            loadPiston();
         }
         else
         {
            std::cout << "MAXIMUM COMPRESSION REACHED. UNLOADING PISTON" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // DECOMPRESSION CYCLE
      if (stage == 6)
      {
         if (pistonForce > 1.0e-5*radiusParticle*keParticle)
         {
            unloadPiston();
         }
         else
         {
            std::cout << "Piston unloaded." << std::endl;
            std::cout << "QUITTING" << std::endl << std::endl;
            setTimeMax(getTime() + getTimeStep());
         }
      }
   }

   void actionsAfterSolve() override
   {
      delete [] boxCentreCoordinates;
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

   void setClustersProperties(int nP, double nu)
   {
      nParticlesPerCluster = nP;
      assumedInterClusterPackingFraction = nu;
   }

   void setLatticeDimensions(int gridLength, int nL, double sizeRatio)
   {
      nLayers = nL;
      gridSize = gridLength;
      clusterSizeToLatticeRatio = sizeRatio;

      numberOfClusters = gridSize*gridSize*nLayers;
      totalNumberOfParticles = numberOfClusters*nParticlesPerCluster;

      std::cout << "Number of clusters needed " << numberOfClusters << std::endl;
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

   void setParticlePlasticProperties(double k1, double k2max, double kCCluster, double kCLoose, double phi)
   {
      kpParticle = k1;
      keParticle = k2max;
      kcParticleIntercluster = kCCluster;
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

   void setCompressionParameters(double minRatio)
   {
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
      speciesParticle -> setCohesionStiffness(kcParticleIntercluster);
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
         particlePureSpeciesVector[i] -> setCohesionStiffness(kcParticleIntercluster);
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

      std::cout << "tC_BB/dt: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle)/getTimeStep() << "\n\n";
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
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setCohesionStiffness(kcParticleIntercluster);
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

      std::cout << "tC_BW/dt: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle)/getTimeStep() << "\n\n";
   }

   void setBoxSize()
   {
      clusterRadius = radiusParticle*pow(nParticlesPerCluster/assumedInterClusterPackingFraction,1.0/3.0);
      boxSize = clusterRadius/clusterSizeToLatticeRatio;
      gridSemiDiagonal = 0.5*sqrt(2.0)*gridSize*boxSize;

      std::cout << "Cluster average radius " << clusterRadius << std::endl;
      std::cout << "Cubic lattice size " << boxSize << std::endl << std::endl;
   }

   void setDomainLimits()
   {
      setXMin(-gridSemiDiagonal);
      setYMin(-gridSemiDiagonal);
      setZMin(0.0);

      setXMax(gridSemiDiagonal);
      setYMax(gridSemiDiagonal);
      setZMax(nLayers*boxSize);
   }

   void makeBoundaries()
   {
      wallHandler.clear();

      // z walls
      base.setSpecies(speciesWall);
      base.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,getZMin()));
      wallHandler.copyAndAddObject(base);

      cylinder.setSpecies(speciesWall);
      cylinder.setPosition(Vec3D(0.0,0.0,0.0));
      cylinder.setOrientation(Vec3D(0.0,0.0,1.0));
      cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(gridSemiDiagonal,0.0,0.0));
      wallHandler.copyAndAddObject(cylinder);
   }

   void setLatticeGrid()
   {
      boxCentreCoordinates = new Vec3D[gridSize*gridSize*nLayers];

      for (int k = 0; k < nLayers; k++) // loop along Z
      {
         for (int j = 0; j < gridSize; j++) // loop along y
         {
            for (int i = 0; i < gridSize; i++) // loop along x
            {
               boxCentreCoordinates[i + gridSize*j + gridSize*gridSize*k].X = (i - 0.5*(gridSize - 1.0))*boxSize;
               boxCentreCoordinates[i + gridSize*j + gridSize*gridSize*k].Y = (j - 0.5*(gridSize - 1.0))*boxSize;
               boxCentreCoordinates[i + gridSize*j + gridSize*gridSize*k].Z = (k + 0.5)*boxSize;
            }
         }
      }
   }

   bool particleInsertionSuccessful(int nCluster)
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
         // particlePosition.X = boxCentreCoordinates[nCluster].X + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         // particlePosition.Y = boxCentreCoordinates[nCluster].Y + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         // particlePosition.Z = boxCentreCoordinates[nCluster].Z + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);

         rad = random.getRandomNumber(p0.getRadius(),0.5*boxSize - 1.01*(p0.getRadius()));
         theta = constants::pi*random.getRandomNumber(-1.0,1.0);
         phi = 0.5*constants::pi*random.getRandomNumber(-1.0,1.0);
         particlePosition.X = boxCentreCoordinates[nCluster].X + rad*sin(theta)*cos(phi);
         particlePosition.Y = boxCentreCoordinates[nCluster].Y + rad*sin(theta)*sin(phi);
         particlePosition.Z = boxCentreCoordinates[nCluster].Z + rad*cos(theta);

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
      for (int i = 0; i < numberOfClusters; i++)
      {
         nParticlesInserted = 0;

         while(nParticlesInserted < nParticlesPerCluster)
         {
            if (particleInsertionSuccessful(i)) nParticlesInserted++;
            else t0 = getTime();
         }

         std::cout << "Inserted cluster n. " << i + 1 << "/" << numberOfClusters << std::endl;
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
   }

   void applyCentralForce()
   {
      Vec3D distanceFromForceCenter;
      int xLatticePosition;
      int yLatticePosition;
      int zLatticePosition;

      for (int j=0; j<numberOfClusters; j++) // count the particles belonging to every single cluster and point them towards their box centre
      {
         for (int i=0; i<nParticlesPerCluster; i++)
         {
            distanceFromForceCenter = particleHandler.getObject(i + nParticlesPerCluster*j) -> getPosition() - boxCentreCoordinates[j];

            //particleHandler.getObject(i) -> addForce(-forceScaleModulus*Vec3D(0.2*distanceFromForceCenter.X,0.2*distanceFromForceCenter.Y,distanceFromForceCenter.Z));
            particleHandler.getObject(i + nParticlesPerCluster*j) -> addForce(-forceModulus*distanceFromForceCenter);
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
      pistonForce = 0.0;

      int totalInteractionCounter = 0;

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

         // piston interactions
         if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
         {
            pistonForce -= ((*i) -> getForce()).Z;
         }

         totalInteractionCounter++;
      }

      meanTotalRelativeOverlap /= totalInteractionCounter;
   }

   // creates the matrix containing the interaction information: "1" -> in contact, "0" -> not in contact
   void setupCoordinationMatrix()
   {
      // resets the number of binding interactions (the bounds after the cluster formation)
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
      // matrixOutput.open("Clusters_INTERACTIONMATRIX", std::ios::out);
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
      // resets the number of binding interactions (the bounds after the cluster formation)
      numberOfBindingInteractions = 0;

      // // prints the matrix before
      // std::ofstream matrixOutput2;
      // matrixOutput2.open("Clusters_COHESIONMATRIX_BEFORE", std::ios::out);
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
            if (coordinationMatrix[i][j]) mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleIntercluster);
            else mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleLoose);
         }
      }

      // // print the cohesive interaction matrix
      // std::ofstream matrixOutput;
      // matrixOutput.open("Clusters_COHESIONMATRIX_AFTER", std::ios::out);
      // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      // {
      //    for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) matrixOutput << mixedSpeciesMatrix[i][j] -> getCohesionStiffness() << "\t";
      //    matrixOutput << std::endl;
      // }
      // matrixOutput.close();
   }

   void makeCompressionPiston()
   {
      pistonStepDisplacement = (1.0 - sizeDispersityParticle)*radiusParticle*1.0e-5;
      initialBedHeight = particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z + (1.0 + sizeDispersityParticle)*radiusParticle;
      pistonHeight = initialBedHeight + 0.1*radiusParticle;

      piston.setSpecies(speciesWall);
      piston.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(piston);
   }

   void loadPiston()
   {
      pistonHeight -= pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }

   void unloadPiston()
   {
      pistonHeight += pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }

   // creates the data output file and writes the first row
   void makeCdatFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "time \t Eel \t Ekin/Eel \t h_piston \t Fpiston" << std::endl;
   }

   // writes the compression data to the output file
   void writeToCdatFile()
   {
      cdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      pistonHeight << "   " <<
      pistonForce << "   " <<
      std::endl;
   }


   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      if (stage < 4)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(5) << getTime() << ", delta_mean = " << meanTotalRelativeOverlap <<
         ", Force Modulus = " << forceModulus << ", g_y = " << getGravity().Z << std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(5) << getTime() <<
         ", h = " << std::setprecision(6) << std::left << std::setw(10) << pistonHeight << ", F_piston = " << pistonForce <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfBindingInteractions << std::endl;
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

   // clusters
   int numberOfClusters;
   int nParticlesPerCluster;
   int totalNumberOfParticles;
   double clusterRadius;
   double assumedInterClusterPackingFraction;
   double granularPackingFraction;
   double meanClusterRadius;
   Vec3D *boxCentreCoordinates;

   // geometry
   double boxSize;
   double clusterSizeToLatticeRatio;
   int nLayers;
   int gridSize;
   double gridSemiDiagonal;
   InfiniteWall base;

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
   double kcParticleIntercluster;
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
   double cdatOutputTimeInterval;
   std::ofstream cdatFile;

   // data analysis
   double meanCoordinationNumber;
   double maxTotalRelativeOverlap;
   double meanTotalRelativeOverlap;

   // compression
   InfiniteWall piston;
   InfiniteWall *pistonPointer;
   AxisymmetricIntersectionOfWalls cylinder;
   double minimumRelativeHeightCompression;
   double pistonForce;
   double pistonStepDisplacement;
   double pistonHeight;
   double initialBedHeight;

   // global
   int stage;
   double t0;
   double energyRatioTolerance;
   double gravityAndForceTuningDuration;
};


int main(int argc, char *argv[])
{
   // TIME STEP
   double timeStep = 5.0e-6;

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;
   double sizeDispersionParticles = 0.0;
   double densityParticles = 1500.0;
   double densityWall = 5000.0;

   // INTERACTION PARAMETERS
   double plasticStiffness = 500.0;
   double elasticStiffness = 1000.0;
   double interclusterCohesionStiffness = 500.0;
   double looseParticleCohesionStiffness = 100.0;
   double phiParticle = 0.3;
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
   double finalForceModulus = 0.00001*particleRadius*elasticStiffness;
   double forceTuningInterval = 0.002;
   double velocityDampingModulus = 0.9;
   double velocityDampingInterval = 0.01;

   // CLUSTERS PARAMETERS
   int nParticlesPerCluster = 1000;
   int gridLength = 1;
   int nLayers = 5;
   double assumedInterClusterPackingFraction = 1.0;
   double clusterSizeToLatticeRatio = 0.25;

   // COMPRESSION PARAMETERS
   double minimumRelativeHeightCompression = 0.3; // it compresses up to an amount equal to minimumRelativeHeightCompression*initialBedHeight

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 1.0e-4;
   double gravityAndForceTuningDuration = 0.2;

   // INITIALIZATION
   ClusterTableting problem;

   problem.setTimeStep(timeStep);
   problem.setTimeMax(100.0);
   problem.setGravity(Vec3D(0.00,0.00,0.00));
   problem.setSystemDimensions(3);
   problem.setCdatOutputTimeInterval(0.001);
   problem.setEnergyRatioTolerance(energyRatioTolerance);
   problem.setSaveCount(0.001/problem.getTimeStep());

   problem.setClustersProperties(nParticlesPerCluster, assumedInterClusterPackingFraction);
   problem.setLatticeDimensions(gridLength, nLayers, clusterSizeToLatticeRatio);

   problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
   problem.setWallDensity(densityWall);

   problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
   problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
   problem.setParticlePlasticProperties(plasticStiffness, elasticStiffness, interclusterCohesionStiffness, looseParticleCohesionStiffness, phiParticle);
   problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
   problem.setWallStiffnessAndRestitutionCoefficients(elasticStiffnessWalls, particleWallRestitutionCoefficient);

   problem.setForceAndVelocityProperties(maximumForceModulus, finalForceModulus, forceTuningInterval, velocityDampingModulus, velocityDampingInterval);
   problem.setGravityAndForceTuningDuration(gravityAndForceTuningDuration);

   problem.setCompressionParameters(minimumRelativeHeightCompression);

   // NAME SETTING
   std::ostringstream name;
   name.str("");
   name.clear();
   std::cout.unsetf(std::ios::floatfield);
   name << "dbClusters_tableting___nP_" << nParticlesPerCluster << "_" << nLayers << "_" << gridLength << "_kP_" << plasticStiffness << "_kE_" << elasticStiffness
   << "_kCi_" << interclusterCohesionStiffness << "_kCe_" << looseParticleCohesionStiffness << "_phi_" << phiParticle << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle;
   problem.setName(name.str());

   problem.solve();

   return 0;
}

// XBALLS ARGUMENTS TO ADD
// -h 800 -p 1 -s 8 -v0 -3dturn 3
