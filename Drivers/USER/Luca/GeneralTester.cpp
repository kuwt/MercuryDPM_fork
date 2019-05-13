/*
 *** dbClusters AGGLOMERATION - EXPLICIT ***
 A single agglomerate is created and its properties computed.
 The EXPLICIT framework means that DEM microscopic parameters are set and the simulation run, regardless of the feasibility of the end result.
 */

#include "NewInteractions/LinearPlasticViscoelasticFrictionSpeciesExtended.h"
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <fstream>
#include <vector>

/*
   - Granules axes ratio is left inside to switch from cubic shape to parallelepipedal

   ToDo:
   - find a way to check the internal structure
   - substitute the check to switch to dissipation stage with a wall-force related one
   - change the hardcoded maximum normal angle of planes with a user-defined one (basically back-compute the max angle, keep theta \in [0,2\pi] and r \in [0,1] and change only the cone height)
   - for the internal structure create random coarse-graining points inside the volume
   - add a rand-srand function for the angle determination
   - ADD THE COHESION INTERACTION REFRESH
   - right now IT FAILS WHEN aspectRatio > 2! Design a new initial bounding box for this case.
 */

class dbClusters_polyhedron : public Mercury3D
{
private:
   void setupInitialConditions() override
   {
      stage = 0;
      nParticlesInserted = 0;
      numberOfBindingInteractions = 0;
      t0 = getTime();

      std::cout << "SETTING SPECIES VECTOR..." << std::endl;
      setSpeciesVector();

      std::cout << "SETTING MIXED SPECIES MATRIX..." << std::endl;
      setMixedSpeciesCohesionMatrix();

      std::cout << "SETTING CUBIC LATTICE SIZE..." << std::endl;
      setBoxSize();

      std::cout << "SETTING DOMAIN LIMITS... " << std::endl;
      setDomainLimits();

      std::cout << "COMPUTING WALLS PARAMETERS... " << std::endl;
      computeWallsParameters();

      std::cout << "CREATING EXTERNAL WALLS... " << std::endl;
      makeBoundaries();

      std::cout << std::endl << "PARTICLE INSERTION" << std::endl;
      insertParticles();

      std::cout << std::endl << "COMPUTING PARTICLES VOLUME..." << std::endl;
      computeTotalParticleVolume();

      std::cout << "CREATING .cdat FILE" << std::endl << std::endl;
      makeCdatFile();

      std::cout << "STARTING COMPRESSION" << std::endl;

      stage++;
   }

   void actionsOnRestart() override
   {

   }

   void actionsAfterTimeStep() override
   {
      // DATA ANALYSIS
      if (stage >= 1) makeDataAnalysis();
      if (stage >= 1 && fmod(getTime(), 0.001) < getTimeStep()) writeToCdatFile();

      // COHESION BOUNDS REFRESH
      if (stage >= 4 && fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) refreshAdjacencyAndCohesionMatrices();

      // COMPRESSION STAGE
      if (stage == 1)
      {
         closeWalls();

         if (fmod(getTime() - t0,velocityDampingInterval) < getTimeStep()) dampVelocities();

         if (meanRelativeOverlap > 1.2*phiParticle/(1.0 - kpParticle/keParticle))
         {
            std::cout << "Average critical overlap reached." << std::endl;
            std::cout << "STARTING DECOMPRESSION" << std::endl << std::endl;

            t0 = getTime();
            stage++;
         }
      }

      // DECOMPRESSION STAGE
      if (stage == 2)
      {
         openWalls();

         if (fmod(getTime() - t0,velocityDampingInterval) < getTimeStep()) dampVelocities();

         if (averageForceOnWalls < 1.0e-3*radiusParticle*keParticle) // wallXpMidpoint.X > 0.5*boxSize*b
         {
            std::cout << "Wall unload reached." << std::endl;
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
         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
         {
            std::cout << "ENERGY DISSIPATED. COMPUTING INTERNAL STRUCTURE AND QUITTING" << std::endl << std::endl;
            // if (computeMassFraction_bool) computeInternalStructure();

            std::cout << "REMOVING COMPRESSION WALLS AND CREATING BASE AND CYLINDER" << std::endl << std::endl;
            makeSettlingWalls();

            // setTimeMax(getTime() + finalStandbyTime);
            t0 = getTime();
            stage++;
         }
      }

      // ACTIVATING GRAVITY, SETTLING AND CREATION OF THE PISTON
      if (stage == 4)
      {
         if (getTime() - t0 < gravityTuningDuration)
         {
            activateGravity();
         }

         if (getTime() - t0 > gravityTuningDuration && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
         {
            std::cout << "CREATING PISTON AND STARTING COMPRESSION" << std::endl << std::endl;
            makeCompressionPiston();

            std::cout << "Base height: " << baseHeight << std::endl;
            std::cout << "Piston height: " << pistonHeight << std::endl;
            std::cout << "Final target piston height: " << baseHeight + initialClusterHeight*minimumRelativeHeightCompression << std::endl;

            // setTimeMax(getTime() + finalStandbyTime); stage = 10;
            stage++;
         }
      }

      // COMPRESSION CYCLE
      if (stage == 5)
      {
         // if (pistonForce < maximumCompressiveForce) // uncomment this for force-driven compression
         if (pistonHeight > baseHeight + initialClusterHeight*minimumRelativeHeightCompression) // uncomment this for height-driven compression
         {
            loadPiston();
         }
         else
         {
            std::cout << "MAXIMUM COMPRESSION REACHED. STANDBY PHASE AND QUITTING" << std::endl << std::endl;
            setTimeMax(getTime() + finalStandbyTime);
            t0 = getTime();
            stage++;
         }
      }
   }

   void actionsAfterSolve() override
   {
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

   void setClustersProperties(int nP, double sizeRatio, double axesRatio, bool parallelFacesFlag)
   {
      totalNumberOfParticles = nP;
      clusterSizeToLatticeRatio = sizeRatio;
      clusterAxesRatio = axesRatio;
      parallelFaces = parallelFacesFlag;

      a = pow(clusterAxesRatio, 2.0/3.0);
      b = pow(clusterAxesRatio, -1.0/3.0);
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

   void setParticlePlasticProperties(double k1, double k2max, double kCCluster, double kCLoose, double phi, double epsilon, bool irreversible)
   {
      kpParticle = k1;
      keParticle = k2max;
      kcParticleIntercluster = kCCluster;
      kcParticleLoose = kCLoose;
      phiParticle = phi;
      relativeDeformationThreshold = epsilon;
      irreversibleCohesiveInteraction = irreversible;
   }

   void setVelocityTuningProperties(double vdm, double vdi)
   {
      velocityDampingModulus = vdm;
      velocityDampingInterval = vdi;
   }

   void setFinalStandbyTimeDuration(double dt)
   {
      finalStandbyTime = dt;
   }

   void setInternalStructureAnalysisParameters(bool compute, bool print, int length, double ratio)
   {
      computeMassFraction_bool = compute;
      exportGridData_bool = print;
      gridLength = length;
      fictiousGridPointRadiusRatio = ratio;
   }

   void setGravityTuningDuration(double dt)
   {
      gravityTuningDuration = dt;
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
      speciesParticle = new LinearPlasticViscoelasticFrictionSpeciesExtended;
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
      speciesWall = new LinearPlasticViscoelasticFrictionSpeciesExtended;
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
      speciesWall = new LinearPlasticViscoelasticFrictionSpeciesExtended;
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
         particlePureSpeciesVector[i] = new LinearPlasticViscoelasticFrictionSpeciesExtended;
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
         std::vector<LinearPlasticViscoelasticFrictionMixedSpeciesExtended*> temporaryRowVector;
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
      // the box size is set by computing the volume of all the particles in a cubic lattice with a packing fraction of 0.1: Vc = N*Vp/0.10
      boxSize = pow(10.0*totalNumberOfParticles*volumeParticle, 1.0/3.0);

      std::cout << "Cubic lattice size " << boxSize << std::endl;
      std::cout << "a = " << a << ", b = " << b << std::endl << std::endl;
   }

   void setDomainLimits()
   {
      setXMin(-0.5*boxSize*b*constants::sqrt_2);
      setYMin(-0.5*boxSize*b*constants::sqrt_2);
      setZMin(-0.5*boxSize*a*constants::sqrt_2);

      setXMax(0.5*boxSize*b*constants::sqrt_2);
      setYMax(0.5*boxSize*b*constants::sqrt_2);
      setZMax(0.5*boxSize*a*constants::sqrt_2);
   }

   void computeWallsParameters()
   {
      double conicalProjectorHeight = 0.5*constants::sqrt_3*2.0;
      double rad, theta;

      // X+ wall
      // random.randomise();
      rad = random.getRandomNumber(0.0, 1.0);
      theta = 2.0*constants::pi*random.getRandomNumber(0.0, 1.0);

      wallXpMidpoint.setZero();
      wallXpMidpoint.X = 0.5*boxSize*b;

      wallXpNormal.X = conicalProjectorHeight;
      wallXpNormal.Y = rad*sin(theta);
      wallXpNormal.Z = rad*cos(theta);
      wallXpNormal /= wallXpNormal.getLength();

      // X- wall
      // if (parallelFaces) random.randomise();
      rad = random.getRandomNumber(0.0, 1.0);
      theta = 2.0*constants::pi*random.getRandomNumber(0.0, 1.0);

      wallXnMidpoint.setZero();
      wallXnMidpoint.X = -0.5*boxSize*b;

      wallXnNormal.X = -conicalProjectorHeight;
      wallXnNormal.Y = rad*cos(theta);
      wallXnNormal.Z = rad*sin(theta);
      wallXnNormal /= wallXnNormal.getLength();

      // Y+ wall
      // if (parallelFaces) random.randomise();
      rad = random.getRandomNumber(0.0, 1.0);
      theta = 2.0*constants::pi*random.getRandomNumber(0.0, 1.0);

      wallYpMidpoint.setZero();
      wallYpMidpoint.Y = 0.5*boxSize*b;

      wallYpNormal.X = rad*cos(theta);
      wallYpNormal.Y = conicalProjectorHeight;
      wallYpNormal.Z = rad*sin(theta);
      wallYpNormal /= wallYpNormal.getLength();

      // Y- wall
      // if (parallelFaces) random.randomise();
      rad = random.getRandomNumber(0.0, 1.0);
      theta = 2.0*constants::pi*random.getRandomNumber(0.0, 1.0);

      wallYnMidpoint.setZero();
      wallYnMidpoint.Y = -0.5*boxSize*b;

      wallYnNormal.X = rad*sin(theta);
      wallYnNormal.Y = -conicalProjectorHeight;
      wallYnNormal.Z = rad*cos(theta);
      wallYnNormal /= wallYnNormal.getLength();

      // Z+ wall
      // if (parallelFaces) random.randomise();
      rad = random.getRandomNumber(0.0, 1.0);
      theta = 2.0*constants::pi*random.getRandomNumber(0.0, 1.0);

      wallZpMidpoint.setZero();
      wallZpMidpoint.Z = 0.5*boxSize*a;

      wallZpNormal.X = rad*sin(theta);
      wallZpNormal.Y = rad*cos(theta);
      wallZpNormal.Z = conicalProjectorHeight;
      wallZpNormal /= wallZpNormal.getLength();

      // Z- wall
      // if (parallelFaces) random.randomise();
      rad = random.getRandomNumber(0.0, 1.0);
      theta = 2.0*constants::pi*random.getRandomNumber(0.0, 1.0);

      wallZnMidpoint.setZero();
      wallZnMidpoint.Z = -0.5*boxSize*a;

      wallZnNormal.X = rad*cos(theta);
      wallZnNormal.Y = rad*sin(theta);
      wallZnNormal.Z = -conicalProjectorHeight;
      wallZnNormal /= wallZnNormal.getLength();

      // infos such as box size and external walls normals are savet to a .wdat file
      std::ofstream wdatFile;
      std::ostringstream wdatName;
      std::cout.unsetf(std::ios::floatfield);
      wdatName << getName() << ".wdat";

      wdatFile.open(wdatName.str(), std::ios::out);
      wdatFile << "P(x)_mid : " << wallXpMidpoint << std::endl;
      wdatFile << "Wall X+ : " << wallXpNormal << std::endl;
      wdatFile << "Wall X- : " << wallXnNormal << std::endl << std::endl;
      wdatFile << "P(y)_mid : " << wallYpMidpoint << std::endl;
      wdatFile << "Wall Y+ : " << wallYpNormal << std::endl;
      wdatFile << "Wall Y- : " << wallYnNormal << std::endl << std::endl;
      wdatFile << "P(z)_mid : " << wallZpMidpoint << std::endl;
      wdatFile << "Wall Z+ : " << wallZpNormal << std::endl;
      wdatFile << "Wall Z- : " << wallZnNormal << std::endl << std::endl;
      wdatFile.close();
   }

   void makeBoundaries()
   {
      wallHandler.clear();
      InfiniteWall wall;
      wall.setSpecies(speciesWall);

      // x walls
      wall.set(wallXpNormal, wallXpMidpoint);
      wallXp = wallHandler.copyAndAddObject(wall);
      wall.set(wallXnNormal, wallXnMidpoint);
      wallXn = wallHandler.copyAndAddObject(wall);

      // y walls
      wall.set(wallYpNormal, wallYpMidpoint);
      wallYp = wallHandler.copyAndAddObject(wall);
      wall.set(wallYnNormal, wallYnMidpoint);
      wallYn = wallHandler.copyAndAddObject(wall);

      // z walls
      wall.set(wallZpNormal, wallZpMidpoint);
      wallZp = wallHandler.copyAndAddObject(wall);
      wall.set(wallZnNormal, wallZnMidpoint);
      wallZn = wallHandler.copyAndAddObject(wall);
   }

   bool particleInsertionSuccessful()
   {
      int insertionFailCounter = 0;
      SphericalParticle p0;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));

      p0.setRadius(radiusParticle*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.0,1.0)));
      p0.setSpecies(particlePureSpeciesVector[particleHandler.getNumberOfObjects()]);

      do
      {
         p0.setPosition(Vec3D(0.5*boxSize*b*random.getRandomNumber(-1.0, 1.0), 0.5*boxSize*b*random.getRandomNumber(-1.0, 1.0), 0.5*boxSize*a*random.getRandomNumber(-1.0, 1.0)));

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
         else
         {
            std::cout << std::endl << "WARNING: THE PARTICLE LOADING PHASE COULD NOT BE ULTIMATED!" << std::endl;
            exit(0);
         }
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl;
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

   void closeWalls()
   {
      wallXpMidpoint.X -= 0.00005*radiusParticle*b;
      wallXp -> set(wallXpNormal, wallXpMidpoint);
      wallXnMidpoint.X += 0.00005*radiusParticle*b;
      wallXn -> set(wallXnNormal, wallXnMidpoint);

      wallYpMidpoint.Y -= 0.00005*radiusParticle*b;
      wallYp -> set(wallYpNormal, wallYpMidpoint);
      wallYnMidpoint.Y += 0.00005*radiusParticle*b;
      wallYn -> set(wallYnNormal, wallYnMidpoint);

      wallZpMidpoint.Z -= 0.00005*radiusParticle*a;
      wallZp -> set(wallZpNormal, wallZpMidpoint);
      wallZnMidpoint.Z += 0.00005*radiusParticle*a;
      wallZn -> set(wallZnNormal, wallZnMidpoint);
   }

   void openWalls()
   {
      wallXpMidpoint.X += 0.00005*radiusParticle*b;
      wallXp -> set(wallXpNormal, wallXpMidpoint);
      wallXnMidpoint.X -= 0.00005*radiusParticle*b;
      wallXn -> set(wallXnNormal, wallXnMidpoint);

      wallYpMidpoint.Y += 0.00005*radiusParticle*b;
      wallYp -> set(wallYpNormal, wallYpMidpoint);
      wallYnMidpoint.Y -= 0.00005*radiusParticle*b;
      wallYn -> set(wallYnNormal, wallYnMidpoint);

      wallZpMidpoint.Z += 0.00005*radiusParticle*a;
      wallZp -> set(wallZpNormal, wallZpMidpoint);
      wallZnMidpoint.Z -= 0.00005*radiusParticle*a;
      wallZn -> set(wallZnNormal, wallZnMidpoint);
   }

   void makeSettlingWalls()
   {
      // clearing every wall
      wallHandler.clear();

      baseHeight = particleHandler.getLowestPositionComponentParticle(2) -> getPosition().Z - 1.5*(1.0 + sizeDispersityParticle)*radiusParticle;

      base.setSpecies(speciesWall);
      base.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,baseHeight));
      basePointer = wallHandler.copyAndAddObject(base);

      cylinder.setSpecies(speciesWall);
      cylinder.setPosition(Vec3D(0.0,0.0,0.0));
      cylinder.setOrientation(Vec3D(0.0,0.0,1.0));
      cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(0.5*boxSize*b*constants::sqrt_2,0.0,0.0));
      wallHandler.copyAndAddObject(cylinder);
   }

   void activateGravity()
   {
      setGravity(Vec3D(0.00,0.00,-9.81*(getTime() - t0)/gravityTuningDuration));
   }

   void makeDataAnalysis()
   {
      int totalInteractionCounter = 0;

      meanCoordinationNumber = 0.0;
      maximumRelativeOverlap = 0.0;
      meanRelativeOverlap = 0.0;
      averageForceOnWalls = 0.0;
      pistonForce = 0.0;

      Vec3D averageCompressionWallsForce;
      averageCompressionWallsForce.setZero();

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
      }
      meanCoordinationNumber /= particleHandler.getNumberOfObjects();

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         meanRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maximumRelativeOverlap) maximumRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

         // agglomeration walls interaction
         if (stage < 3)
         {
            if ((*i) -> getI() -> getIndex() == wallXp -> getIndex()) averageCompressionWallsForce -= (*i) -> getForce();
            if ((*i) -> getI() -> getIndex() == wallXn -> getIndex()) averageCompressionWallsForce += (*i) -> getForce();

            if ((*i) -> getI() -> getIndex() == wallYp -> getIndex()) averageCompressionWallsForce -= (*i) -> getForce();
            if ((*i) -> getI() -> getIndex() == wallYn -> getIndex()) averageCompressionWallsForce += (*i) -> getForce();

            if ((*i) -> getI() -> getIndex() == wallZp -> getIndex()) averageCompressionWallsForce -= (*i) -> getForce();
            if ((*i) -> getI() -> getIndex() == wallZn -> getIndex()) averageCompressionWallsForce += (*i) -> getForce();
         }

         // piston interactions
         if (stage >= 5)
         {
            if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex())
            {
               pistonForce -= ((*i) -> getForce()).Z;
            }
         }

         totalInteractionCounter++;
      }
      meanRelativeOverlap /= totalInteractionCounter;

      averageForceOnWalls = averageCompressionWallsForce.getLength();
      averageForceOnWalls /= 6.0;
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
   }

   // refreshes the matrix containing the interaction information and updates the cohesion coefficients accordingly
   void refreshAdjacencyAndCohesionMatrices()
   {
      // resets the number of bounds between particles (after the cluster formation)
      numberOfBindingInteractions = 0;

      // resets the matrix content to 0
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         for (int j = 0; j < particleHandler.getNumberOfObjects(); j++)
         {
            coordinationMatrix[i][j] = 0;
         }
      }

      // variable for the reduced radius computation
      double reducedRadius;

      // re-allocates the contacts in the coordination matrix
      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         reducedRadius = 0.0;

         // the contact is checked if it was in place at the previous check
         if (coordinationMatrixPreviousCheck[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()])
         {
            // // computing the reduced radius for this interaction
            // reducedRadius = 2.0*((particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius())*(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()));
            // reducedRadius /= ((particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius()) + (particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()));

            // using average particle radius for threshold computation

            // if the relative overlap with respect to the reduced radius is above the threshold the contact is preserved
            if (((*i) -> getOverlap())/radiusParticle > phiParticle*relativeDeformationThreshold)
            {
               coordinationMatrix[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
               coordinationMatrix[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;
               numberOfBindingInteractions++;
            }
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
            else
            {
               mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleLoose);

               if (irreversibleCohesiveInteraction)
               {
                  mixedSpeciesMatrix[i][j] -> setStiffnessAndRestitutionCoefficient(kpParticle, particleParticleRestitutionCoeff, massParticle);
                  mixedSpeciesMatrix[i][j] -> setSlidingStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getLoadingStiffness()*2.0/7.0);
                  mixedSpeciesMatrix[i][j] -> setSlidingDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getDissipation()*2.0/7.0);
                  mixedSpeciesMatrix[i][j] -> setRollingStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getLoadingStiffness()*2.0/7.0);
                  mixedSpeciesMatrix[i][j] -> setRollingDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getDissipation()*2.0/7.0);
                  mixedSpeciesMatrix[i][j] -> setTorsionStiffness(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getLoadingStiffness()*2.0/7.0);
                  mixedSpeciesMatrix[i][j] -> setTorsionDissipation(speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> getDissipation()*2.0/7.0);
               }
            }
         }
      }
   }

   // creates the data output file and writes the first row
   void makeCdatFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "time \t Eel \t Ekin/Eel \t coord_number \t n_bounds \t dTotMean \t dTotMax" << std::endl;
   }

   // writes the compression data to the output file
   void writeToCdatFile()
   {
      cdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      meanCoordinationNumber << "   " <<
      numberOfBindingInteractions << "   " <<
      meanRelativeOverlap << "   " <<
      maximumRelativeOverlap << "   " <<
      pistonForce << "   " <<
      std::endl;
   }

   void makeCompressionPiston()
   {
      pistonStepDisplacement = (1.0 - sizeDispersityParticle)*radiusParticle*1.0e-5;
      pistonHeight = particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z + 1.5*(1.0 + sizeDispersityParticle)*radiusParticle;

      initialClusterHeight = particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z - particleHandler.getLowestPositionComponentParticle(2) -> getPosition().Z + 2.0*radiusParticle;

      piston.setSpecies(speciesWall);
      piston.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(piston);
   }

   void loadPiston()
   {
      pistonHeight -= pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }


   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      std::cout << "t = " << std::setprecision(3) << std::left << std::setw(4) << getTime() << ", tmax = " << getTimeMax() <<
      ", cN = " << meanCoordinationNumber << ", nB = " << numberOfBindingInteractions <<
      ", dMean = " << meanRelativeOverlap << ", dMax = " << maximumRelativeOverlap <<
      ", xP = " << wallXpMidpoint.X << ", F_walls = " << averageForceOnWalls <<
      ", F = " << pistonForce <<
      std::endl;

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

   // clusters
   int totalNumberOfParticles;
   Vec3D centerOfMass;

   // cluster shape
   double clusterAxesRatio;
   bool parallelFaces;

   // geometry
   double boxSize;
   double a, b;
   double clusterSizeToLatticeRatio;
   InfiniteWall *wallXp, *wallXn;
   InfiniteWall *wallYp, *wallYn;
   InfiniteWall *wallZp, *wallZn;
   Vec3D wallXpMidpoint, wallXnMidpoint;
   Vec3D wallYpMidpoint, wallYnMidpoint;
   Vec3D wallZpMidpoint, wallZnMidpoint;
   Vec3D wallXpNormal, wallXnNormal;
   Vec3D wallYpNormal, wallYnNormal;
   Vec3D wallZpNormal, wallZnNormal;
   double averageForceOnWalls;

   // compression
   InfiniteWall piston;
   InfiniteWall* pistonPointer;
   InfiniteWall base;
   InfiniteWall* basePointer;
   AxisymmetricIntersectionOfWalls cylinder;
   double baseHeight;
   double minimumRelativeHeightCompression;
   double pistonForce;
   double pistonStepDisplacement;
   double pistonHeight;
   double initialClusterHeight;

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
   double relativeDeformationThreshold;
   double particleWallRestitutionCoeff;
   double particleParticleRestitutionCoeff;
   double densityWall;
   bool irreversibleCohesiveInteraction;

   // velocity damping
   double velocityDampingModulus;
   double velocityDampingInterval;

   // contact related
   int numberOfBindingInteractions;
   std::vector< std::vector<int> > coordinationMatrix;
   std::vector< std::vector<int> > coordinationMatrixPreviousCheck;

   // species
   LinearPlasticViscoelasticFrictionSpeciesExtended *speciesParticle, *speciesWall;
   std::vector<LinearPlasticViscoelasticFrictionSpeciesExtended*> particlePureSpeciesVector;
   std::vector< std::vector<LinearPlasticViscoelasticFrictionMixedSpeciesExtended*> > mixedSpeciesMatrix;

   // output
   double cdatOutputTimeInterval;
   std::ofstream cdatFile;
   std::ofstream gridFile;

   // data analysis
   double meanCoordinationNumber;
   double maximumRelativeOverlap;
   double meanRelativeOverlap;
   double massFraction;
   int gridLength;
   double gridResolutionLength;
   double fictiousGridPointRadiusRatio;

   // global
   int stage;
   double t0;
   double energyRatioTolerance;
   double finalStandbyTime;
   double gravityTuningDuration;
   bool computeMassFraction_bool;
   bool exportGridData_bool;
};


int main(int argc, char *argv[])
{
   // TIME STEP
   double timeStep = 5.0e-6;

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;
   double sizeDispersionParticles = 0.0;
   double densityParticles = 1500.0;
   double densityWall = 3000.0;

   // INTERACTION PARAMETERS
   double elasticStiffness = 1000.0;
   double plasticStiffness = 500.0;
   const std::vector<double> interclusterCohesionStiffness = {500.0, 1000.0, 2000.0, 5000.0}; // looping index: i
   const std::vector<double> phiParticle = {0.3, 0.2, 0.1}; // looping index: j
   const std::vector<double> relativeDeformationThreshold = {0.8, 0.85, 0.9, 0.95}; // looping index: k
   double looseParticleCohesionStiffness = 0.0;
   double particleParticleRestitutionCoefficient = 0.5;
   double particleWallRestitutionCoefficient = 0.8;
   double elasticStiffnessWalls = elasticStiffness;
   bool irreversibleInteraction = true;

   // FRICTION PARAMETERS
   double muSlidingParticleParticle = 0.5;
   double muRollingParticleParticle = 0.3;
   double muSlidingParticleWall = 0.5;
   double muRollingParticleWall = 0.3;

   // FORCE AND DAMPING PARAMETERS
   double velocityDampingModulus = 0.9;
   double velocityDampingInterval = 500.0*timeStep;

   // CLUSTER PARAMETERS
   int numberOfParticlesPerCluster = 200;
   double clusterSizeToLatticeRatio = 0.5;

   // SHAPE PARAMETERS
   double axesRatio = 1.0;
   bool parallelFaces = true;

   // INTERNAL STRUCTURE ANALYSIS PARAMETERS
   int gridLength = 300;
   double fictiousGridParticleRadiusRatio = 1.0e-5;

   // COMPRESSION PARAMETERS
   double minimumRelativeHeightCompression = 0.5; // it compresses up to an amount equal to minimumRelativeHeightCompression*cluster_size

   // GLOBAL PARAMETERS
   bool computeMassFraction = true;
   bool exportGridData = true;
   double energyRatioTolerance = 1.0e-4;
   double finalStandbyTime = 0.1;
   double gravityTuningDuration = 0.2;

   for (int i = 0; i < interclusterCohesionStiffness.size(); i++)
   {
      for (int j = 0; j < phiParticle.size(); j++)
      {
         for (int k = 0; k < relativeDeformationThreshold.size(); k++)
         {
            // INITIALIZATION
            dbClusters_polyhedron problem;

            problem.setTimeStep(timeStep);
            problem.setTimeMax(10.0);
            problem.setGravity(Vec3D(0.00,0.00,0.00));
            problem.setSystemDimensions(3);
            problem.setCdatOutputTimeInterval(0.001);
            problem.setEnergyRatioTolerance(energyRatioTolerance);
            problem.setSaveCount(0.001/problem.getTimeStep());

            problem.setClustersProperties(numberOfParticlesPerCluster, clusterSizeToLatticeRatio, axesRatio, parallelFaces);
            problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
            problem.setWallDensity(densityWall);
            problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
            problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
            problem.setParticlePlasticProperties(plasticStiffness, elasticStiffness, interclusterCohesionStiffness[i], looseParticleCohesionStiffness, phiParticle[j], relativeDeformationThreshold[k], irreversibleInteraction);
            problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
            problem.setWallStiffnessAndRestitutionCoefficients(elasticStiffnessWalls, particleWallRestitutionCoefficient);
            problem.setVelocityTuningProperties(velocityDampingModulus, velocityDampingInterval);
            problem.setFinalStandbyTimeDuration(finalStandbyTime);
            problem.setInternalStructureAnalysisParameters(computeMassFraction, exportGridData, gridLength, fictiousGridParticleRadiusRatio);
            problem.setGravityTuningDuration(gravityTuningDuration);
            problem.setCompressionParameters(minimumRelativeHeightCompression);

            // NAME SETTING
            std::ostringstream name;
            name.str("");
            name.clear();
            std::cout.unsetf(std::ios::floatfield);
            name << "general_polyhedra_tester";
            // name << "breakage_tester___kC_" << interclusterCohesionStiffness[i] << "_phi_" << phiParticle[j] << "_epsilon_" << relativeDeformationThreshold[k];
            problem.setName(name.str());

            problem.solve();

            if (i == 0 && j == 0 && k == 0) exit(0);
         }
      }
   }

   return 0;
}
