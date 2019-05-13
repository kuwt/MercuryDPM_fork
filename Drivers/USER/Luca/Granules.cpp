#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>

/*
* ToDo:
* - make the box size dependent on the particle max size
* - add force dependence on particle volume (for isotropical big/small distribution)
* - modify and correct non-isotropical forces cases
* - now that it works remove the vector and species matrix and replace all with a cohesion matrix for the mixed collision
* - provide an input for the dt between contact checks (to alter the contact-cohesive matrix)
*/

// 10.05.18 - switched to single-component particles (removed small particle entity)


class Granules : public Mercury3D
{
private:
   void setupInitialConditions() override
   {
      stage = 1;
      nParticlesInserted = 0;
      forceModulus = 0.0;
      t0 = getTime();

      std::cout << "SETTING SPECIES..." << std::endl;
      // setSpecies();

      std::cout << "SETTING SPECIES VECTOR..." << std::endl;
      setSpeciesVector();

      std::cout << "SETTING MIXED SPECIES MATRIX..." << std::endl;
      setMixedSpeciesCohesionMatrix();

      std::cout << "COMPUTING BOX SIZE... " << std::endl;
      computeBoxSize();

      std::cout << "SETTING DOMAIN LIMITS... " << std::endl;
      setDomainLimits();

      std::cout << "CREATING BOUNDARIES..." << std::endl;
      makeBoundaries();

      std::cout << std::endl << "PARTICLE INSERTION" << std::endl;
      insertParticles();

      std::cout << "SETTING COMPRESSION TARGET HEIGHT" << std::endl;
      setPistonTargetHeight();
      std::cout << "TOTAL HEIGHT " << getZMax() << std::endl;
      std::cout << "TARGET HEIGHT " << pistonTargetHeight << std::endl;
      std::cout << "TARGET/BOX_SIZE " << pistonTargetHeight/boxSize << std::endl;

      std::cout << "ACTIVATING CENTRAL FORCES" << std::endl;

      stage++;
   }

   void actionsOnRestart() override
   {

   }

   void actionsAfterTimeStep() override
   {
      // ACTIVATION OF LOCAL FORCE
      if (stage == 2)
      {
         applyCentralForce();
         if (getTime() < 0.2) increaseForce();

         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > 0.2)
         {
            std::cout << "DAMPING FORCES" << std::endl << std::endl;
            stage++;
         }
      }

      // DAMPING OF LOCAL FORCE
      if (stage == 3)
      {
         applyCentralForce();
         dampForce();

         if (forceModulus < 1.0e-2)
         {
            std::cout << "CREATING COORDINATION MATRIX AND REFRESHING COORDINATION MATRIX" << std::endl << std::endl;
            setupCoordinationMatrix();
            refreshCoordinationMatrixAndCohesionMatrix();

            std::cout << "TUNING GRAVITY" << std::endl << std::endl;
            makeFloor();
            t0 = getTime();
            stage++;
         }
      }

      // GRAVITY TUNING
      if (stage == 4)
      {
         tuneGravity();

         if (getTime() > t0 + gravityTuningDuration)
         {
            std::cout << "SETTLING" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // PISTON CREATION
      if (stage == 5)
      {
         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > t0 + 0.1)
         {
            std::cout << "CREATING PISTON" << std::endl << std::endl;
            makePiston();
            std::cout << "Piston initial height " << pistonHeight << std::endl;
            std::cout << "Piston target height " << pistonTargetHeight << std::endl;
            std::cout << "Piston velocity " << pistonVelocity << std::endl << std::endl;

            std::cout << "CREATING CDAT FILE" << std::endl << std::endl;
            makeCdatFile();

            std::cout << "COMPRESSION" << std::endl << std::endl;
            stage++;
         }
      }

      // DATA ANALYSIS AND OUTPUT TO CDAT
      if (stage >= 6) makeDataAnalysis();
      if (stage >= 6 && fmod(getTime(), 0.001) < getTimeStep()) writeToCdatFile();

      // UPDATE OF COORDINATION MATRIX AND COHESION MATRIX
      if (stage >= 6 && fmod(getTime(), 0.001) < getTimeStep()) refreshCoordinationMatrixAndCohesionMatrix();

      // COMPRESSION
      if (stage == 6)
      {
         if (pistonHeight > pistonTargetHeight)
         {
            movePiston();
         }
         else
         {
            std::cout << "PISTON TARGET REACHED. DISSIPATING ENERGY." << std::endl << std::endl;
            stage++;
         }
      }

      // ENERGY DISSIPATION
      if (stage == 7)
      {
         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
         {
            std::cout << "ENERGY DISSIPATED. UNLOADING PISTON." << std::endl << std::endl;
            pistonVelocity *= -1.0;

            stage++;
         }
      }

      // PISTON UNLOADING
      if (stage == 8)
      {
         if (pistonHeight < (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 2.0*radiusParticle)
         {
            movePiston();
         }
         else
         {
            std::cout << "PISTON UNLOADED. QUITTING." << std::endl << std::endl;
            setTimeMax(getTime() + getTimeStep());
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

   void setGranulesProperties(int nG, int nP)
   {
      numberOfGranules = nG;
      nParticlesPerGranule = nP;
      totalNumberOfParticles = numberOfGranules*nParticlesPerGranule;
   }

   //void setBoxSize(double size)
   //{
   //boxSize = size;
   //}

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
      k1Wall = kW;
      particleWallRestitutionCoeff = eB;
   }

   void setParticleParticleRestitutionCoefficients(double bigBigE)
   {
      particleParticleRestitutionCoeff = bigBigE;
   }

   void setParticlePlasticProperties(double k1, double k2max, double kC, double phi)
   {
      k1Particle = k1;
      k2MaxParticle = k2max;
      kCParticle = kC;
      phiParticle = phi;
   }

   void setForceAndVelocityProperties(double fsm, double fdm, double fdi, double vdm, double vdi)
   {
      forceScaleModulus = fsm;
      forceDampingModulus = fdm;
      forceDampingInterval = fdi;
      velocityDampingModulus = vdm;
      velocityDampingInterval = vdi;
   }

   void setGravityTuningDuration(double dt)
   {
      gravityTuningDuration = dt;
   }

   // void setCompressionRatio(double ratio)
   // {
   //    looseBedCompressionRatio = ratio;
   // }

   void setRelativeFinalDeformation(double ratio)
   {
      relativeFinalDeformation = ratio;
   }



   // FUNCTIONS CALLED IN THE CLASS -----------------------------------

   void setSpecies()
   {
      speciesHandler.clear();

      // BIG-BIG
      speciesParticle = new LinearPlasticViscoelasticFrictionSpecies;
      speciesParticle -> setDensity(densityParticle);
      speciesParticle -> setStiffnessAndRestitutionCoefficient(k1Particle, particleParticleRestitutionCoeff, massParticle);
      speciesParticle -> setUnloadingStiffnessMax(k2MaxParticle);
      speciesParticle -> setCohesionStiffness(kCParticle);
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
      speciesWall -> setStiffnessAndRestitutionCoefficient(k1Wall, 1.0, massParticle);
      speciesWall -> setUnloadingStiffnessMax(k1Wall);
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
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Particle + k1Wall), particleWallRestitutionCoeff, massParticle);
      speciesHandler.getMixedObject(speciesParticle, speciesWall) -> setUnloadingStiffnessMax(k2MaxParticle);
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
      std::cout << "BIG-BIG stiffness and dissipation coefficients: " << speciesParticle -> getLoadingStiffness() << " " << speciesParticle -> getDissipation() << "\n";
      std::cout << "BIG-BIG friction coefficients: " << particleParticleSlidingFrictionCoeff << " " << particleParticleRollingFrictionCoeff << " " << particleParticleTorsionFrictionCoeff << "\n";
      std::cout << "BIG-BIG tangential stiffnesses: " << speciesParticle -> getSlidingStiffness() << " " << speciesParticle -> getRollingStiffness() << " " << speciesParticle -> getTorsionStiffness() << "\n";
      std::cout << "BIG-BIG tangential dissipation coefficients: " << speciesParticle -> getSlidingDissipation() << " " << speciesParticle -> getRollingDissipation() << " " << speciesParticle -> getTorsionDissipation() << "\n";
      std::cout << "BIG-BIG collision time: " << std::setprecision(4) << speciesParticle -> getCollisionTime(massParticle) << "\n\n";

      std::cout << "BIG-WALL stiffness and dissipation coefficients: " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getLoadingStiffness() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getDissipation() << "\n";
      std::cout << "BIG-WALL friction coefficients: " << particleWallSlidingFrictionCoeff << " " << particleWallRollingFrictionCoeff << " " << particleWallTorsionFrictionCoeff << "\n";
      std::cout << "BIG-WALL tangential stiffnesses: " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getSlidingStiffness() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getRollingStiffness() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getTorsionStiffness() << "\n";
      std::cout << "BIG-WALL tangential dissipation coefficients: " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getSlidingDissipation() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getRollingDissipation() << " " << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getTorsionDissipation() << "\n";
      std::cout << "BIG-WALL collision time: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getCollisionTime(massParticle) << "\n\n";

      std::cout << "tC_BB/dt: " << std::setprecision(4) << speciesParticle -> getCollisionTime(massParticle)/getTimeStep() << "\n";
      std::cout << "tC_BW/dt: " << std::setprecision(4) << speciesHandler.getMixedObject(speciesParticle, speciesWall) -> getCollisionTime(massParticle)/getTimeStep() << "\n";
   }

   void setSpeciesVector()
   {
      speciesHandler.clear();

      // speciesParticle = new LinearPlasticViscoelasticFrictionSpecies;
      // speciesParticle -> setDensity(densityParticle);
      // speciesParticle -> setStiffnessAndRestitutionCoefficient(k1Particle, particleParticleRestitutionCoeff, massParticle);
      // speciesParticle -> setUnloadingStiffnessMax(k2MaxParticle);
      // speciesParticle -> setCohesionStiffness(kCParticle);
      // speciesParticle -> setPenetrationDepthMax(phiParticle);
      //
      // speciesParticle -> setSlidingFrictionCoefficient(particleParticleSlidingFrictionCoeff);
      // speciesParticle -> setSlidingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      // speciesParticle -> setSlidingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      // speciesParticle -> setRollingFrictionCoefficient(particleParticleRollingFrictionCoeff);
      // speciesParticle -> setRollingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      // speciesParticle -> setRollingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      // speciesParticle -> setTorsionFrictionCoefficient(particleParticleTorsionFrictionCoeff);
      // speciesParticle -> setTorsionStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      // speciesParticle -> setTorsionDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      // speciesHandler.addObject(speciesParticle);

      // WALL-WALL
      speciesWall = new LinearPlasticViscoelasticFrictionSpecies;
      speciesWall -> setDensity(densityWall);
      speciesWall -> setStiffnessAndRestitutionCoefficient(k1Wall, 1.0, massParticle);
      speciesWall -> setUnloadingStiffnessMax(k1Wall);
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
         particlePureSpeciesVector[i] -> setStiffnessAndRestitutionCoefficient(k1Particle, particleParticleRestitutionCoeff, massParticle);
         particlePureSpeciesVector[i] -> setUnloadingStiffnessMax(k2MaxParticle);
         particlePureSpeciesVector[i] -> setCohesionStiffness(kCParticle);
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
      std::cout << "BIG stiffness and dissipation coefficients: " << particlePureSpeciesVector[0] -> getLoadingStiffness() << " " << particlePureSpeciesVector[0] -> getDissipation() << "\n";
      std::cout << "BIG friction coefficients: " << particleParticleSlidingFrictionCoeff << " " << particleParticleRollingFrictionCoeff << " " << particleParticleTorsionFrictionCoeff << "\n";
      std::cout << "BIG tangential stiffnesses: " << particlePureSpeciesVector[0] -> getSlidingStiffness() << " " << particlePureSpeciesVector[0] -> getRollingStiffness() << " " << particlePureSpeciesVector[0] -> getTorsionStiffness() << "\n";
      std::cout << "BIG tangential dissipation coefficients: " << particlePureSpeciesVector[0] -> getSlidingDissipation() << " " << particlePureSpeciesVector[0] -> getRollingDissipation() << " " << particlePureSpeciesVector[0] -> getTorsionDissipation() << "\n";
      std::cout << "BIG collision time: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle) << "\n";
      std::cout << "tC_BB/dt: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle)/getTimeStep() << "\n\n";

      std::cout << "WALL stiffness and dissipation coefficients: " << speciesWall -> getLoadingStiffness() << " " << speciesWall -> getDissipation() << "\n\n";
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
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setStiffnessAndRestitutionCoefficient(0.5*(k1Particle + k1Wall), particleWallRestitutionCoeff, massParticle);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], speciesWall) -> setUnloadingStiffnessMax(k2MaxParticle);
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
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setUnloadingStiffnessMax(k2MaxParticle);
               speciesHandler.getMixedObject(particlePureSpeciesVector[i], particlePureSpeciesVector[j]) -> setCohesionStiffness(kCParticle);
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
      std::cout << "BIG-WALL stiffness and dissipation coefficients: " << mixedSpeciesMatrix[0][0] -> getLoadingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getDissipation() << "\n";
      std::cout << "BIG-WALL friction coefficients: " << particleWallSlidingFrictionCoeff << " " << particleWallRollingFrictionCoeff << " " << particleWallTorsionFrictionCoeff << "\n";
      std::cout << "BIG-WALL tangential stiffnesses: " << mixedSpeciesMatrix[0][0] -> getSlidingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getRollingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getTorsionStiffness() << "\n";
      std::cout << "BIG-WALL tangential dissipation coefficients: " << mixedSpeciesMatrix[0][0] -> getSlidingDissipation() << " " << mixedSpeciesMatrix[0][0] -> getRollingDissipation() << " " << mixedSpeciesMatrix[0][0] -> getTorsionDissipation() << "\n";
      std::cout << "BIG-WALL collision time: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle) << "\n";
      std::cout << "tC_BW/dt: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle)/getTimeStep() << "\n\n";
   }

   void computeBoxSize()
   {
      boxVolume = 20.0*nParticlesPerGranule*volumeParticle;
      boxSize = pow(boxVolume,1.0/3.0);
   }

   void setDomainLimits()
   {
      setXMin(-0.5*boxSize);
      setYMin(-0.5*boxSize);
      setZMin(0.0);

      setXMax(0.5*boxSize);
      setYMax(0.5*boxSize);
      setZMax(3.0*numberOfGranules*boxSize);
   }

   void makeBoundaries()
   {
      wallHandler.clear();

      xBoundary.set(Vec3D(1.0,0.0,0.0), getXMin(), getXMax());
      boundaryHandler.copyAndAddObject(xBoundary);

      yBoundary.set(Vec3D(0.0,1.0,0.0), getYMin(), getYMax());
      boundaryHandler.copyAndAddObject(yBoundary);

      zBoundary.set(Vec3D(0.0,0.0,1.0), getZMin(), getZMax());
      boundaryHandler.copyAndAddObject(zBoundary);
   }

   bool particleInsertionSuccessful(int nGranule)
   {
      int insertionFailCounter = 0;
      Vec3D particlePosition;
      SphericalParticle p0;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));

      p0.setRadius(radiusParticle*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.0,1.0)));
      p0.setSpecies(particlePureSpeciesVector[particleHandler.getNumberOfObjects()]);
      // p0.setSpecies(speciesParticle);

      do
      {
         particlePosition.X = (getXMax() - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Y = (getYMax() - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Z = (nGranule + 0.5)*boxSize + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);

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
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
   }

   void applyCentralForce()
   {
      Vec3D distanceFromForceCenter;
      int zLatticePosition;

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         zLatticePosition = (int)((particleHandler.getObject(i) -> getPosition().Z)/boxSize);

         distanceFromForceCenter.X = (particleHandler.getObject(i) -> getPosition().X);
         distanceFromForceCenter.Y = (particleHandler.getObject(i) -> getPosition().Y);
         distanceFromForceCenter.Z = (particleHandler.getObject(i) -> getPosition().Z) - boxSize*(0.5 + zLatticePosition);

         //particleHandler.getObject(i) -> addForce(-forceScaleModulus*Vec3D(0.2*distanceFromForceCenter.X,0.2*distanceFromForceCenter.Y,distanceFromForceCenter.Z));
         particleHandler.getObject(i) -> addForce(-forceModulus*distanceFromForceCenter);
         //particleHandler.getObject(i) -> addForce(-forceScaleModulus*distanceFromForceCenter*(particleHandler.getObject(i) -> getRadius())/(0.5*(radiusParticle + radiusSmall)));
         //particleHandler.getObject(i) -> addForce(-1.0*distanceFromForceCenter/(distanceFromForceCenter.getLength()*distanceFromForceCenter.getLength())/10000);

         if (fmod(getTime(),velocityDampingInterval) < getTimeStep())
         {
            particleHandler.getObject(i) -> setVelocity(velocityDampingModulus*(particleHandler.getObject(i) -> getVelocity()));
         }
      }
   }

   void increaseForce()
   {
      forceModulus = forceScaleModulus*getTime()/0.2;
   }

   void dampForce()
   {
      if (fmod(getTime(),forceDampingInterval) < getTimeStep())
      {
         forceModulus *= forceDampingModulus;
      }
   }

   void makeFloor()
   {
      wallHandler.clear();
      base.setSpecies(speciesWall);
      base.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,getZMin()));
      wallHandler.copyAndAddObject(base);

      top.setSpecies(speciesWall);
      top.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,getZMax()));
      wallHandler.copyAndAddObject(top);

      //wall.setSpecies(speciesWall);
      //wall.set(Vec3D(-1.,0.,0.),Vec3D(0.8*getXMin(),0.,0.));
      //wallHandler.copyAndAddObject(wall);
      //wall.set(Vec3D(1.,0.,0.),Vec3D(0.8*getXMax(),0.,0.));
      //wallHandler.copyAndAddObject(wall);

      //wall.set(Vec3D(0,-1.,0.),Vec3D(0.,0.8*getYMin(),0.));
      //wallHandler.copyAndAddObject(wall);
      //wall.set(Vec3D(0.,1.,0.),Vec3D(0.,0.8*getYMax(),0.));
      //wallHandler.copyAndAddObject(wall);

      AxisymmetricIntersectionOfWalls(cylinder);
      cylinder.setSpecies(speciesWall);
      cylinder.setPosition(Vec3D(0.0,0.0,0.0));
      cylinder.setOrientation(Vec3D(0.0,0.0,1.0));
      cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(0.5*boxSize,0.0,0.0));
      wallHandler.copyAndAddObject(cylinder);
   }

   void tuneGravity()
   {
      if (getTime() < t0 + gravityTuningDuration) setGravity(Vec3D(0.0,0.0,-9.81*(getTime() - t0)/gravityTuningDuration));
   }

   void setPistonTargetHeight()
   {
      // this was the old setting: compression based on volume
      // double totalParticleVolume = 0.0;
      // for (int i=0; i<particleHandler.getNumberOfObjects(); i++)
      // {
      //    totalParticleVolume += pow(particleHandler.getObject(i) -> getRadius(),3.0);
      // }
      // totalParticleVolume *= 4.0*constants::pi/3.0;
      // pistonTargetHeight = looseBedCompressionRatio*totalParticleVolume/(0.6*pow(boxSize,2.0));

      // this is the new setting: compression based on prescribed final relative deformation
      pistonTargetHeight = relativeFinalDeformation*(particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z;
   }

   void makePiston()
   {
      std::cout << std::endl << "CREATING PISTON" << std::endl << std::endl;

      pistonHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z + 1.5*radiusParticle;
      pistonVelocity = -radiusParticle/10000.0/getTimeStep();

      piston.setSpecies(speciesWall);
      piston.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(piston);

      setPistonTargetHeight();
   }

   void movePiston()
   {
      pistonHeight += pistonVelocity*getTimeStep();
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
      pistonPointer -> setVelocity(Vec3D(0.0,0.0,pistonVelocity));
   }

   void makeDataAnalysis()
   {
      meanCoordinationNumber = 0.0;

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
      }

      meanCoordinationNumber /= particleHandler.getNumberOfObjects();

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
         if ((*i) -> getI() -> getIndex() == base.getIndex())
         {
            basePressure += ((*i) -> getForce()).Z;

            meanBaseRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxBaseRelativeOverlap) maxBaseRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            baseInteractionCounter++;
         }

         // cylinder interactions
         if ((*i) -> getI() -> getIndex() == cylinder.getIndex())
         {
            cylinderPressure -= (((*i) -> getContactPoint()).X * ((*i) -> getForce()).X + ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).Y)/(0.5*boxSize);

            meanCylinderRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxCylinderRelativeOverlap) maxCylinderRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

            cylinderInteractionCounter++;
         }

         totalInteractionCounter++;
      }

      pistonPressure /= constants::pi*pow(0.5*boxSize,2.0);
      basePressure /= constants::pi*pow(0.5*boxSize,2.0);
      cylinderPressure /= 2.0*constants::pi*(0.5*boxSize)*pistonHeight;
      meanTotalRelativeOverlap /= totalInteractionCounter;
      meanPistonRelativeOverlap /= pistonInteractionCounter;
      meanBaseRelativeOverlap /= baseInteractionCounter;
      meanCylinderRelativeOverlap /= cylinderInteractionCounter;
   }

   // creates the matrix containing the interaction information: "1" -> in contact, "0" -> not in contact
   void setupCoordinationMatrix()
   {
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
      numberOfCohesiveInteractions = 0;
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
            numberOfCohesiveInteractions++;
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
            if (coordinationMatrix[i][j]) mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kCParticle);
            else mixedSpeciesMatrix[i][j] -> setCohesionStiffness(0.0);
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
      cdatFile << "time \t Eel \t Ekin/Eel \t h \t v \t P \t Pbase \t Pcyl \t coord. number \t dTotMean \t dTotMax \t dPistonMean \t dPistonMax \t dBaseMean \t dBaseMax \t k1 \t k2 \t kC \t phi \t STAGE" << std::endl;
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
      cylinderPressure << "   " <<
      meanCoordinationNumber << "   " <<
      meanTotalRelativeOverlap << "   " <<
      maxTotalRelativeOverlap << "   " <<
      meanPistonRelativeOverlap << "   " <<
      maxPistonRelativeOverlap << "   " <<
      meanBaseRelativeOverlap << "   " <<
      maxBaseRelativeOverlap << "   " <<
      k1Particle << "   " <<
      k2MaxParticle << "   " <<
      kCParticle << "   " <<
      phiParticle << "   " <<
      stage <<
      std::endl;
   }



   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      if (stage < 6)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", Eratio = " << std::setprecision(4) << std::left << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
         ", Force Modulus = " << forceModulus << ", g_z = " << getGravity().Z << ", STAGE " << std::setprecision(0) << stage
         << std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", Elastic Energy = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Eratio = " << getKineticEnergy()/getElasticEnergy() << ", h = " << pistonHeight <<
         ", v = " << pistonVelocity << std::endl <<
         "# cohesive interactions = " << numberOfCohesiveInteractions << std::endl <<
         "P = " << std::setprecision(6) << std::left << std::setw(10) << pistonPressure << ", Pbase = " << basePressure << ", Pcylinder = " << cylinderPressure << std::endl <<
         "cN = " << std::setprecision(6) << std::left << std::setw(10) << meanCoordinationNumber << ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
         ", <d_pist> = " << meanPistonRelativeOverlap << ", max(d_pist) = " << maxPistonRelativeOverlap << ", <d_base> = " << meanBaseRelativeOverlap <<
         ", max(d_base) = " << maxBaseRelativeOverlap << ", <d_cyl> = " << meanCylinderRelativeOverlap << ", max(d_cyl) = " << maxCylinderRelativeOverlap << std::endl <<
         "k1 = " << std::setprecision(6) << std::left << k1Particle << ", k2 = " << k2MaxParticle << ", kC = " << kCParticle << ", phi = " << phiParticle << std::endl << std::endl;
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

   // geometry
   double boxSize;
   double boxVolume;
   //InfiniteWall wall;
   InfiniteWall base;
   InfiniteWall top;
   InfiniteWall piston;
   InfiniteWall *pistonPointer;
   AxisymmetricIntersectionOfWalls cylinder;

   // friction
   double particleWallSlidingFrictionCoeff;
   double particleWallRollingFrictionCoeff;
   double particleWallTorsionFrictionCoeff;
   double particleParticleSlidingFrictionCoeff;
   double particleParticleRollingFrictionCoeff;
   double particleParticleTorsionFrictionCoeff;

   // collision
   double k1Wall, k1Particle;
   double k2MaxParticle;
   double kCParticle;
   double phiParticle;
   double particleWallRestitutionCoeff;
   double particleParticleRestitutionCoeff;
   double densityWall;

   // central force
   double forceModulus;
   double forceScaleModulus;
   double forceDampingModulus;
   double forceDampingInterval;
   double velocityDampingModulus;
   double velocityDampingInterval;

   // compression
   double pistonHeight;
   double pistonVelocity;
   double pistonTargetHeight;
   // double looseBedCompressionRatio;
   double relativeFinalDeformation;

   // contact related
   int numberOfCohesiveInteractions;
   std::vector< std::vector<int> > coordinationMatrix;
   std::vector< std::vector<int> > coordinationMatrixPreviousCheck;

   // species
   LinearPlasticViscoelasticFrictionSpecies *speciesParticle, *speciesWall;
   std::vector<LinearPlasticViscoelasticFrictionSpecies*> particlePureSpeciesVector;
   std::vector< std::vector<LinearPlasticViscoelasticFrictionMixedSpecies*> > mixedSpeciesMatrix;

   // boundaries
   PeriodicBoundary xBoundary, yBoundary, zBoundary;

   // output
   double cdatOutputTimeInterval;
   std::ofstream cdatFile;

   // data analysis
   double pistonPressure;
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

   // global
   int stage;
   double energyRatioTolerance;
   double t0;
   double gravityTuningDuration;
};


int main(int argc, char *argv[])
{
   // TIME STEP
   double timeStep = 1.0e-5;

   // SETUP PARAMETERS
   int numberOfGranules = 1;
   int nParticlesPerGranule = 1000;
   double boxSize = 0.05;
   double particleToBoxSizeRatio = 30.0;
   double sizeDispersionParticles = 0.1;
   double densityParticles = 1452.7;
   double densityWall = 2000.0;

   // INTERACTION PARAMETERS
   double k1Particle = 1000.0;
   double k2MaxParticle = 3000.0;
   double kCParticle = 100.0;
   // double kCParticle[5] = {100.0,200.0,500.0,1000.0,2000.0};
   double phiParticle = 0.1;
   double particleParticleRestitutionCoefficient = 0.5;
   double particleWallRestitutionCoefficient = 0.7;

   double k1Wall = 3000.0;

   // FRICTION PARAMETERS
   double muSlidingParticleParticle = 0.5;
   double muRollingParticleParticle = 0.1;
   double muSlidingParticleWall = 0.5;
   double muRollingParticleWall = 0.1;

   // FORCE AND DAMPING PARAMETERS
   double forceScaleModulus = 1000.0;
   double forceDampingModulus = 0.9;
   double forceDampingInterval = 0.005;
   double velocityDampingModulus = 0.9;
   double velocityDampingInterval = 0.01;

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 5.0e-5;
   double gravityTuningDuration = 0.2;

   // COMPRESSION PARAMETERS
   // double looseBedCompressionRatio = 0.30;
   double relativeFinalDeformation = 0.50;

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

   // for (int i=0; i<5; i++)
   // {
      Granules problem;

      // INITIALIZATION
      problem.setTimeStep(timeStep);
      problem.setTimeMax(10.0);
      problem.setGravity(Vec3D(0.00,0.00,0.00));
      problem.setSystemDimensions(3);
      problem.setCdatOutputTimeInterval(0.001);
      problem.setEnergyRatioTolerance(energyRatioTolerance);
      problem.setSaveCount(0.001/problem.getTimeStep());

      problem.setGranulesProperties(numberOfGranules, nParticlesPerGranule);

      //problem.setBoxSize(boxSize);
      problem.setParticleProperties(boxSize/particleToBoxSizeRatio, sizeDispersionParticles, densityParticles);
      problem.setWallDensity(densityWall);

      problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
      problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
      problem.setParticlePlasticProperties(k1Particle, k2MaxParticle, kCParticle, phiParticle);
      problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
      problem.setWallStiffnessAndRestitutionCoefficients(k1Wall, particleWallRestitutionCoefficient);

      problem.setForceAndVelocityProperties(forceScaleModulus, forceDampingModulus, forceDampingInterval, velocityDampingModulus, velocityDampingInterval);
      problem.setGravityTuningDuration(gravityTuningDuration);

      // problem.setCompressionRatio(looseBedCompressionRatio);
      problem.setRelativeFinalDeformation(relativeFinalDeformation);

      // NAME SETTING
      std::cout.unsetf(std::ios::floatfield);
      name << "Granules___TESTING_nP_" << nParticlesPerGranule << "_nG_" << numberOfGranules << "_relDeformation_" << relativeFinalDeformation <<
      "_F_" << forceScaleModulus << "_k1_" << k1Particle << "_k2_" << k2MaxParticle << "_kC_" << kCParticle << "_phi_" << phiParticle << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle;
      problem.setName(name.str());

      problem.solve();


   return 0;
}


// XBALLS ARGUMENTS TO ADD
// -h 800 -p 1 -s 8 -v0 -3dturn 3
