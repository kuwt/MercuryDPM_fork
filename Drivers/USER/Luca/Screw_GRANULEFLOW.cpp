/*
 *** PERIODIC SCREW FEEDER - GRANULE FLOW ***
 Granules are loaded into a 1-flighted periodic screw section
 Granules are allowed to break via the proper interaction matrix

 *** IMPORTANT ***
 The packing fraction is obtained in teh following way:
 - Computation of teh free volume inside of the screw
 - Computation of the volume of particles inside concentric cylindrical shells of radius R
 - The particle laying such that rP > R when vParticles/vFree == packingFraction are removed
 - The casing is inserted

 LAST UPDATE: 3.7.17
 */

#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Helicoid05.h"
#include "ScrewBottom.h"
#include <fstream>

/*
 ToDo:

 */

class Screw_granularFlow : public Mercury3D
{

private:

   void setupInitialConditions() override
   {
      stage = 1;
      nParticlesInserted = 0;
      forceModulus = 0.0;
      t0 = getTime();

      std::cout << "COMPUTING NUMBER OF GRANULES AND PARTICLES NEEDED..." << std::endl;
      computeNumberOfGranulesNeeded();

      std::cout << "COMPUTING SYSTEM GEOMETRY..." << std::endl;
      computeScrewAndGranuleBoxSize();

      std::cout << "SETTING SPECIES VECTOR..." << std::endl;
      setSpeciesVector();

      std::cout << "SETTING MIXED SPECIES MATRIX..." << std::endl;
      setMixedSpeciesCohesionMatrix();

      std::cout << "SETTING DOMAIN LIMITS... " << std::endl;
      setDomainLimits();

      std::cout << "CREATING BOUNDARIES..." << std::endl;
      makeBoundaries();

      std::cout << "CREATING LATTICE GRID" << std::endl;
      setLatticeGrid();

      std::cout << std::endl << "PARTICLE INSERTION" << std::endl;
      insertParticles();

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
         else
         {
            std::cout << "DAMPING FORCES" << std::endl << std::endl;
            stage++;
         }

         // if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > 0.2)
         // {
         //    std::cout << "DAMPING FORCES" << std::endl << std::endl;
         //    stage++;
         // }
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

            std::cout << "CREATING THE SCREW" << std::endl << std::endl;
            makeScrew();

            std::cout << "TUNING GRAVITY" << std::endl << std::endl;
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

      // ENERGY DISSIPATION
      if (stage == 5)
      {
         if (getKineticEnergy()/getElasticEnergy() < 0.00001 || getTime() > t0 + 0.2)
         {
            std::cout << "ENERGY DISSIPATED. CREATING DATA FILES." << std::endl << std::endl;
            makeCdatFile();

            std::cout << "CREATING SCREW CASING." << std::endl << std::endl;
            makeScrewCasing();

            std::cout << "INITIATING SCREW ROTATION" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // DATA ANALYSIS, UPDATE OF COORDINATION MATRIX AND COHESION MATRIX AND OUTPUT TO CDAT
      if (stage >= 6) makeDataAnalysis();
      if (stage >= 6 && fmod(getTime(), 0.001) < getTimeStep()) refreshCoordinationMatrixAndCohesionMatrix();
      if (stage >= 6 && fmod(getTime(), 0.001) < getTimeStep()) writeToCdatFile();

      // SCREW ROTATION
      if (stage == 6)
      {
         if (getTime() < screwRunningTime + t0) rotateScrew();
         else
         {
            std::cout << "SCREW ROTATION TERMINATED. QUITTING..." << std::endl << std::endl;
            setTimeMax(getTime() + getTimeStep());

            stage++;
         }
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

   void setGranulesProperties(int nP, double ratio, double nu)
   {
      nParticlesPerGranule = nP;
      shaftToGranuleSizeRatio = ratio;
      assumedInterGranulePackingFraction = nu;
   }

   void setScrewVariables(double omega, double t, double f, double nu)
   {
      screwAngularVelocity = omega;
      screwRunningTime = t;
      screwFillLevel = f;
      assumedScrewPackingFraction = nu;
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

   void computeNumberOfGranulesNeeded()
   {
      numberOfGranules = (int)(9.0*screwFillLevel*assumedScrewPackingFraction*pow(shaftToGranuleSizeRatio,3.0));
      totalNumberOfParticles = numberOfGranules*nParticlesPerGranule;

      std::cout << "Number of granules needed " << numberOfGranules << std::endl;
      std::cout << "Total number of particles needed " << totalNumberOfParticles << std::endl << std::endl;
   }

   void computeScrewAndGranuleBoxSize()
   {
      std::cout << "For this simulation the geometry ratios are:" << std::endl;
      std::cout << " - box_size = shaft_radius " << std::endl;
      std::cout << " - blade_radius = 2*shaft_radius " << std::endl;
      std::cout << " - screw_length = 4*shaft_radius " << std::endl;
      std::cout << " - granule_box_size = 2*shaft_radius " << std::endl << std::endl;

      granuleRadius = radiusParticle*pow(nParticlesPerGranule/assumedInterGranulePackingFraction,1.0/3.0);

      std::cout << "Granule average radius " << granuleRadius << std::endl;

      screwShaftRadius = shaftToGranuleSizeRatio*granuleRadius;
      screwBladeRadius = 3.0*screwShaftRadius;
      screwCasingRadius = screwBladeRadius;
      screwLength = 2.0*screwBladeRadius;
      screwBladeThickness = 2.0*radiusParticle;
      numberOfScrewTurns = 1;

      boxSize = screwShaftRadius;

      std::cout << "Screw shaft radius " << screwShaftRadius << std::endl;
      std::cout << "Screw blade radius " << screwBladeRadius << std::endl;
      std::cout << "Screw casing radius " << screwCasingRadius << std::endl;
      std::cout << "Screw length " << screwLength << std::endl;
      std::cout << "Screw thickness " << screwBladeThickness << std::endl << std::endl;
   }

   void setDomainLimits()
   {
      setXMin(-1.1*screwCasingRadius);
      setYMin(-1.1*screwCasingRadius);
      setZMin(0.0);

      setXMax(1.1*screwCasingRadius);
      setYMax(1.1*screwCasingRadius + 1.1*numberOfGranules/4*boxSize);
      setZMax(screwLength);
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

      // z boundaries
      zBoundary.set(Vec3D(0.0,0.0,1.0), getZMin(), getZMax());
      boundaryHandler.copyAndAddObject(zBoundary);
   }

   void setLatticeGrid()
   {
      nXlayers = (int)(2*screwCasingRadius/boxSize);
      nZlayers = (int)(screwLength/boxSize);
      nYlayers = (int)(numberOfGranules/(nXlayers*nZlayers)) + 1;
      boxCentreCoordinates = new Vec3D[nXlayers*nZlayers*nYlayers];

      double offsetX = -screwCasingRadius + 0.5*boxSize;
      double offsetY = screwCasingRadius + 0.5*boxSize;
      double offsetZ = 0.5*boxSize;

      std::cout << "n. of X layers " << nXlayers << std::endl;
      std::cout << "n. of Z layers " << nZlayers << std::endl;
      std::cout << "n. of Y layers " << nYlayers << std::endl;
      std::cout << "total n. of boxes " << nXlayers*nZlayers*nYlayers << std::endl;

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
      Vec3D particlePosition;
      SphericalParticle p0;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));

      p0.setRadius(radiusParticle*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.0,1.0)));
      p0.setSpecies(particlePureSpeciesVector[particleHandler.getNumberOfObjects()]);
      // p0.setSpecies(speciesParticle);

      do
      {
         particlePosition.X = boxCentreCoordinates[nGranule].X + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Y = boxCentreCoordinates[nGranule].Y + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);
         particlePosition.Z = boxCentreCoordinates[nGranule].Z + (0.5*boxSize - 1.01*(p0.getRadius()))*random.getRandomNumber(-1.0,1.0);

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

            if (fmod(getTime(),velocityDampingInterval) < getTimeStep())
            {
               particleHandler.getObject(i) -> setVelocity(velocityDampingModulus*(particleHandler.getObject(i) -> getVelocity()));
            }
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

   void tuneGravity()
   {
      if (getTime() < t0 + gravityTuningDuration) setGravity(Vec3D(0.0,-9.81*(getTime() - t0)/gravityTuningDuration,0.0));
   }

   void makeDataAnalysis()
   {
      meanCoordinationNumber = 0.0;

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
      }

      meanCoordinationNumber /= particleHandler.getNumberOfObjects();

      screwForce = 0.0;
      screwTorque = 0.0;
      meanTotalRelativeOverlap = 0.0;
      maxTotalRelativeOverlap = 0.0;
      meanScrewRelativeOverlap = 0.0;
      maxScrewRelativeOverlap = 0.0;

      int totalInteractionCounter = 0;
      int screwInteractionCounter = 0;

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         meanTotalRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxTotalRelativeOverlap) maxTotalRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());

         // screw interactions
         if ((*i) -> getI() -> getIndex() == screwPointer -> getIndex())
         {
            screwForce -= ((*i) -> getForce()).Z;

            meanScrewRelativeOverlap += ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            if (((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) > maxScrewRelativeOverlap) maxScrewRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
            screwTorque += ((*i) -> getTorque()).Z;

            screwInteractionCounter++;
         }

         totalInteractionCounter++;
      }

      meanTotalRelativeOverlap /= totalInteractionCounter;
      meanScrewRelativeOverlap /= screwInteractionCounter;
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

   // makes the screw
   void makeScrew()
   {
      screwOrigin.X = 0.0;
      screwOrigin.Y = 0.0;
      screwOrigin.Z = 0.0;

      // screw parameters
      screw.setSpecies(speciesWall);
      screw.set(screwOrigin, screwLength, screwBladeRadius, screwShaftRadius, numberOfScrewTurns, screwAngularVelocity, screwBladeThickness, true);
      screwPointer = wallHandler.copyAndAddObject(screw);

      bottomCasingSection.setSpecies(speciesWall);
      bottomCasingSection.set(screwOrigin, screwCasingRadius, screwLength);
      bottomCasingSectionPointer = wallHandler.copyAndAddObject(bottomCasingSection);
   }

   // makes the screw casing
   void makeScrewCasing()
   {
      screwCasing.setSpecies(speciesWall);
      screwCasing.setPosition(screwOrigin);
      screwCasing.setOrientation(Vec3D(0.0,0.0,1.0));
      screwCasing.addObject(Vec3D(1.0,0.0,0.0),Vec3D(screwCasingRadius,0.0,0.0));
      screwCasing.setAngularVelocity(Vec3D(0.,0.,0.));
      wallHandler.copyAndAddObject(screwCasing);
   }

   // rotates the screw
   void rotateScrew()
   {
       // the actual rotation of the blade
       screwPointer -> rotate(getTimeStep());

       // applies the proper angular velocity to the screw blade (for a correct collision computation)
       // IMPORTANT: the sign of the angular velocity depends on the rotation verse as well as if the screw is right/left handed
       screwPointer -> setAngularVelocity(Vec3D(0.,0.,-screwAngularVelocity));
       screwPointer -> setOrientation(Vec3D(0.0,0.0,1.0));
   }

   // creates the data output file and writes the first row
   void makeCdatFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "time \t Eel \t Ekin/Eel \t screw_force \t screw_torque \t coord_number \t n_bounds \t dTotMean \t dTotMax \t dScrewMean \t dScrewMax" << std::endl;
   }

   // writes the compression data to the output file
   void writeToCdatFile()
   {
      cdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      screwForce << "   " <<
      screwTorque << "   " <<
      meanCoordinationNumber << "   " <<
      numberOfBindingInteractions << "   " <<
      meanTotalRelativeOverlap << "   " <<
      maxTotalRelativeOverlap << "   " <<
      meanScrewRelativeOverlap << "   " <<
      maxScrewRelativeOverlap << "   " <<
      std::endl;
   }



   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      if (stage < 6)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", E_ratio = " << std::setprecision(4) << std::left << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
         ", Force Modulus = " << forceModulus << ", g_y = " << getGravity().Y << ", STAGE " << std::setprecision(0) << stage
         << std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
         ", Eel = " << std::setprecision(6) << std::left << std::setw(10) << getElasticEnergy() << ", Ekin/Eel = " << getKineticEnergy()/getElasticEnergy() <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfBindingInteractions << std::endl <<
         ", Fs = " << std::setprecision(6) << std::left << std::setw(10) << screwForce << ", Ts = " << screwTorque <<
         ", <d_tot> = " << meanTotalRelativeOverlap << ", max(d_tot) = " << maxTotalRelativeOverlap <<
         ", <d_screw> = " << meanScrewRelativeOverlap << ", max(d_screw) = " << maxScrewRelativeOverlap << std::endl << std::endl;
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
   double shaftToGranuleSizeRatio;
   Vec3D *boxCentreCoordinates;

   // geometry
   double boxSize;
   int nXlayers;
   int nYlayers;
   int nZlayers;
   InfiniteWall wall;

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

   // screw
   double screwShaftRadius;
   double screwBladeRadius;
   double screwCasingRadius;
   double screwLength;
   double screwBladeThickness;
   int numberOfScrewTurns;
   double screwAngularVelocity;
   double screwRunningTime;
   double screwFillLevel;
   double assumedScrewPackingFraction;
   Vec3D screwOrigin;
   Helicoid05 screw;
   Helicoid05 *screwPointer;
   ScrewBottom bottomCasingSection;
   ScrewBottom *bottomCasingSectionPointer;
   AxisymmetricIntersectionOfWalls screwCasing;

   // contact related
   int numberOfBindingInteractions;
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
   double meanCoordinationNumber;
   double screwForce;
   double screwTorque;
   double maxTotalRelativeOverlap;
   double meanTotalRelativeOverlap;
   double maxScrewRelativeOverlap;
   double meanScrewRelativeOverlap;

   // global
   int stage;
   double energyRatioTolerance;
   double t0;
   double gravityTuningDuration;
};


int main(int argc, char *argv[])
{
   // TIME STEP
   double timeStep = 5.0e-6;

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;
   double sizeDispersionParticles = 0.1;
   double densityParticles = 1452.7;
   double densityWall = 2000.0;

   // INTERACTION PARAMETERS
   double k1Particle = 1000.0;
   double k1Wall = 3000.0;
   double k2MaxParticle = 3000.0;
   double kCParticle = 500.0;
   double phiParticle = 0.1;
   double particleParticleRestitutionCoefficient = 0.5;
   double particleWallRestitutionCoefficient = 0.7;

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

   // GRANULES PARAMETERS
   int nParticlesPerGranule = 50;
   double shaftToGranuleSizeRatio = 4.0;
   double assumedInterGranulePackingFraction = 0.80;

   // SCREW PARAMETERS
   double screwAngularVelocity = 4.0*constants::pi;
   double screwFillLevel = 0.5;
   double assumedScrewPackingFraction = 0.60;
   double screwRunningTime = 10.0;

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 5.0e-5;
   double gravityTuningDuration = 0.2;

   // STRING VARIABLES
   std::ostringstream name;
   name.str("");
   name.clear();

   Screw_granularFlow problem;

   // INITIALIZATION
   problem.setTimeStep(timeStep);
   problem.setTimeMax(10.0);
   problem.setGravity(Vec3D(0.00,0.00,0.00));
   problem.setSystemDimensions(3);
   problem.setCdatOutputTimeInterval(0.001);
   problem.setEnergyRatioTolerance(energyRatioTolerance);
   problem.setSaveCount(0.001/problem.getTimeStep());

   problem.setGranulesProperties(nParticlesPerGranule, shaftToGranuleSizeRatio, assumedInterGranulePackingFraction);
   problem.setScrewVariables(screwAngularVelocity, screwRunningTime, screwFillLevel, assumedScrewPackingFraction);

   //problem.setBoxSize(boxSize);
   problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
   problem.setWallDensity(densityWall);

   problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
   problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
   problem.setParticlePlasticProperties(k1Particle, k2MaxParticle, kCParticle, phiParticle);
   problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
   problem.setWallStiffnessAndRestitutionCoefficients(k1Wall, particleWallRestitutionCoefficient);

   problem.setForceAndVelocityProperties(forceScaleModulus, forceDampingModulus, forceDampingInterval, velocityDampingModulus, velocityDampingInterval);
   problem.setGravityTuningDuration(gravityTuningDuration);

   // NAME SETTING
   std::cout.unsetf(std::ios::floatfield);
   name << "Screw_granularFlow___TESTING_nP_" << nParticlesPerGranule << "_lambda_" << shaftToGranuleSizeRatio << "_k1_" << k1Particle << "_k2_" << k2MaxParticle
   << "_kC_" << kCParticle << "_phi_" << phiParticle << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle
   << "_f_" << screwFillLevel << "_omega_" << screwAngularVelocity;
   problem.setName(name.str());

   problem.solve();

   return 0;
}

// XBALLS ARGUMENTS TO ADD
// -h 800 -p 1 -s 8 -v0 -3dturn 3
