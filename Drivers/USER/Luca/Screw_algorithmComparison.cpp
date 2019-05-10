#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "ScrewAuger.h"
#include <fstream>
#include <chrono>
#include <ctime>

/*
Last update : 6.5.19

ToDo:
MAJOR

MINOR

WHATEVER

*/

class ScrewAlgorithmComparison : public Mercury3D
{
private:

   void setupInitialConditions() override
   {
      stage = 0;
      screwVelocityIncrementControlVariable = 0;
      t0 = getTime();

      std::cout << "Setting up species..." << std::endl;
      setSpecies();

      std::cout << "Setting up the filling box dimensions..." << std::endl;
      setFillingBoxDimensions();

      std::cout << "Setting up the number of particles needed (accounts for the shaft volume)..." << std::endl;
      setNumberOfParticles();

      std::cout << "Setting up the loading lattice properties..." << std::endl;
      setloadingBoxGridProperties();

      std::cout << "Setting up the loading box dimensions..." << std::endl;
      setLoadingBoxDimensions();

      std::cout << "Setting up the simulation boundaries..." << std::endl;
      setBoundaries();

      std::cout << "Creating boundaries... " << std::endl;
      makeBoundaries();

      std::cout << "Creating the screw... ";
      makeScrew();

      std::cout << "Creating particles..." << std::endl;
      makeParticles();

      stage++;
   }

   void actionsBeforeTimeStep() override
   {
      if (stage == 1 && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > timeBuffer)
      {
         std::cout << "Creating the screw casing..." << std::endl;
         makeCasing();
         std::cout << "Relaxation of the system..." << std::endl;

         t0 = getTime();
         stage++;
      }

      if (stage == 2 && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > t0 + timeBuffer)
      {
         std::cout << "Cutting the particles according to the filling ratio..." << std::endl;
         makeFillingRatio();
         std::cout << "Relaxation of the system..." << std::endl;

         t0 = getTime();
         stage++;
      }

      if (stage == 3 && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance && getTime() > t0 + timeBuffer)
      {
         std::cout << "Resetting time max...";
         setTimeMax(nVelocityIncrements*runtimePerScrewVelocity + getTime());
         std::cout << "New time max: " << getTimeMax() << std::endl;
         std::cout << "Creating the particle list array..." << std::endl;
         makeParticleList();
         std::cout << "Creating the .cdat output file..." << std::endl;
         makeOutputFile();
         std::cout << "Starting the screw rotation..." << std::endl;

         t0 = getTime();
         stage++;
      }

      if (stage == 4)
      {
         makeRotation();

         if (fmod(getTime(),cdatOutputTimeInterval) < getTimeStep())
         {
            makeDataAnalysis();
            writeDataToOutptFile();
            makeParticleList();
         }
      }
   }

   void actionAfterSolve()
   {
      cdatFile.close();
   }


public:
   // FUNCTIONS CALLED IN MAIN ----------------------------------------
   void setCdatOutputTimeInterval(double dt)
   {
      cdatOutputTimeInterval = dt;
   }

   void setEnergyRatioTolerance(double dE)
   {
      energyRatioTolerance = dE;
   }

   void setTimeBuffer(double t)
   {
      timeBuffer = t;
   }

   void setParticleProperties(double pR, double dR, double rhoParticle)
   {
      particleRadius = pR;
      sizeDispersityParticle = dR;
      densityParticles = rhoParticle;

      particleVolume = 4.*constants::pi*pow(particleRadius,3.)/3.;
      particleMass = densityParticles*particleVolume;
   }

   void setWallDensity(double rhoWall)
   {
      densityWalls = rhoWall;
   }

   void setParticleParticleFrictionCoefficients(double muS, double muR, double muT)
   {
      particleParticleSlidingFriction = muS;
      particleParticleRollingFriction = muR;
      particleParticleTorsionFriction = muT;
   }

   void setParticleWallFrictionCoefficients(double muS, double muR, double muT)
   {
      particleWallSlidingFriction = muS;
      particleWallRollingFriction = muR;
      particleWallTorsionFriction = muT;
   }

   void setInteractionProperties(double kP, double kW, double ePP, double ePW)
   {
      particleStiffness = kP;
      wallsStiffness = kW;
      particleParticleRestitutionCoefficient = ePP;
      particleWallRestitutionCoefficient = ePW;
   }

   void setScrewRelativeGeometryParameters(Vec3D origin, double Rc, double Rb, double Rs, double t, double p, int nPitches)
   {
      screwOrigin = origin;
      rescaledCasingRadius = Rc;
      rescaledBladeRadius = Rb;
      rescaledShaftRadius = Rs;
      rescaledScrewThickness = t;
      rescaledPitchLength = p;
      numberOfScrewPitches = nPitches;

      screwCasingRadius = rescaledCasingRadius*particleRadius;
      screwBladeRadius = rescaledBladeRadius*particleRadius;
      screwShaftRadius = rescaledShaftRadius*particleRadius;
      screwThickness = rescaledScrewThickness*particleRadius;
      screwLength = rescaledPitchLength*numberOfScrewPitches*particleRadius;
   }

   void setScrewOperationalParameters(double time, int n, double *omega, double fR)
   {
      runtimePerScrewVelocity = time;
      screwFillingRatio = fR;
      nVelocityIncrements = n;

      for (int i = 0; i < nVelocityIncrements; i++)
      {
         screwAngularVelocity.push_back(omega[i]);
      }
   }

   void setSetupParameters(double zetaI, double sizeRatio, double deltaP)
   {
      assumedInitialPackingFraction = zetaI;
      settlingBoxToScrewSizeRatio = sizeRatio;
      interParticleLoadingLatticeRelativeDistance = deltaP;

      interParticleLoadingDistance = interParticleLoadingLatticeRelativeDistance*particleRadius*(1.0 + sizeDispersityParticle);
   }

   void developmentHackParticles(bool hackFlag, int nLayers)
   {
      particleHackFlag = hackFlag;
      particleHackLayers = nLayers;
   }


   // FUNCTIONS CALLED IN THE CLASS ----------------------------------------
   void setSpecies()
   {
      speciesHandler.clear();

      // PARTICLE-PARTICLE
      speciesParticle = new LinearViscoelasticFrictionSpecies;
      speciesParticle -> setDensity(densityParticles);
      speciesParticle -> setStiffnessAndRestitutionCoefficient(particleStiffness, particleParticleRestitutionCoefficient, particleMass);

      speciesParticle -> setSlidingFrictionCoefficient(particleParticleSlidingFriction);
      speciesParticle -> setSlidingStiffness(speciesParticle -> getStiffness()*2.0/7.0);
      speciesParticle -> setSlidingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setRollingFrictionCoefficient(particleParticleRollingFriction);
      speciesParticle -> setRollingStiffness(speciesParticle -> getStiffness()*2.0/7.0);
      speciesParticle -> setRollingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setTorsionFrictionCoefficient(particleParticleTorsionFriction);
      speciesParticle -> setTorsionStiffness(speciesParticle -> getStiffness()*2.0/7.0);
      speciesParticle -> setTorsionDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesParticle);

      // WALL-WALL
      speciesWall = new LinearViscoelasticFrictionSpecies;
      speciesWall -> setDensity(densityWalls);
      speciesWall -> setStiffnessAndRestitutionCoefficient(wallsStiffness, 1.0, particleMass);

      speciesWall -> setSlidingFrictionCoefficient(0.0);
      speciesWall -> setSlidingStiffness(speciesWall -> getStiffness()*2.0/7.0);
      speciesWall -> setSlidingDissipation(speciesWall -> getDissipation()*2.0/7.0);
      speciesWall -> setRollingFrictionCoefficient(0.0);
      speciesWall -> setRollingStiffness(speciesWall -> getStiffness()*2.0/7.0);
      speciesWall -> setRollingDissipation(speciesWall -> getDissipation()*2.0/7.0);
      speciesWall -> setTorsionFrictionCoefficient(0.0);
      speciesWall -> setTorsionStiffness(speciesWall -> getStiffness()*2.0/7.0);
      speciesWall -> setTorsionDissipation(speciesWall -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesWall);

      // PARTICLE-WALL
      auto speciesMixedParticleWall = speciesHandler.getMixedObject(speciesParticle, speciesWall);
      speciesMixedParticleWall -> setStiffnessAndRestitutionCoefficient(0.5*(particleStiffness + wallsStiffness), particleWallRestitutionCoefficient, particleMass);

      speciesMixedParticleWall -> setSlidingFrictionCoefficient(particleWallSlidingFriction);
      speciesMixedParticleWall -> setSlidingStiffness(speciesMixedParticleWall -> getStiffness()*2.0/7.0);
      speciesMixedParticleWall -> setSlidingDissipation(speciesMixedParticleWall -> getDissipation()*2.0/7.0);

      speciesMixedParticleWall -> setRollingFrictionCoefficient(particleWallRollingFriction);
      speciesMixedParticleWall -> setRollingStiffness(speciesMixedParticleWall -> getStiffness()*2.0/7.0);
      speciesMixedParticleWall -> setRollingDissipation(speciesMixedParticleWall -> getDissipation()*2.0/7.0);

      speciesMixedParticleWall -> setTorsionFrictionCoefficient(particleWallTorsionFriction);
      speciesMixedParticleWall -> setTorsionStiffness(speciesMixedParticleWall -> getStiffness()*2.0/7.0);
      speciesMixedParticleWall -> setTorsionDissipation(speciesMixedParticleWall -> getDissipation()*2.0/7.0);

      // OUTPUTS
      std::cout << "tC_BB/dt: " << std::setprecision(4) << speciesParticle -> getCollisionTime(particleMass)/getTimeStep() << std::endl;
      std::cout << "tC_BW/dt: " << std::setprecision(4) << speciesMixedParticleWall -> getCollisionTime(particleMass)/getTimeStep() << std::endl << std::endl;
   }

   void setFillingBoxDimensions()
   {
      fillingBoxMin.X = screwOrigin.X - settlingBoxToScrewSizeRatio*screwCasingRadius;
      fillingBoxMin.Y = screwOrigin.Y - settlingBoxToScrewSizeRatio*screwCasingRadius;
      fillingBoxMin.Z = screwOrigin.Z;

      fillingBoxMax.X = screwOrigin.X + settlingBoxToScrewSizeRatio*screwCasingRadius;
      fillingBoxMax.Y = screwOrigin.Y + settlingBoxToScrewSizeRatio*screwCasingRadius;
      fillingBoxMax.Z = screwOrigin.Z + rescaledPitchLength*numberOfScrewPitches*particleRadius;

      fillingBoxVolume = (fillingBoxMax.X - fillingBoxMin.X)*(fillingBoxMax.Y - fillingBoxMin.Y)*(fillingBoxMax.Z - fillingBoxMin.Z);

      std::cout << "Filling box dimensions: " << fillingBoxMax - fillingBoxMin << std::endl;
      std::cout << "Filling box volume: " << fillingBoxVolume << std::endl << std::endl;
   }

   void setNumberOfParticles() // overestimated packing fraction * (box volume - shaft volume) / mean particle volume
   {
      numberOfParticles = (int)(0.7*3.0/(4.0*constants::pi)*rescaledPitchLength*numberOfScrewPitches*(4.0*pow(settlingBoxToScrewSizeRatio*rescaledCasingRadius,2.0) - constants::pi*pow(rescaledShaftRadius,2.0)));

      std::cout << "Number of particles needed: " << numberOfParticles << std::endl;
      std::cout << "Estimated number of particles at full filling: " << (int)(0.7*3.0/4.0*rescaledPitchLength*numberOfScrewPitches*(pow(rescaledCasingRadius,2.0) - pow(rescaledShaftRadius,2.0))) << std::endl << std::endl;
   }

   void setloadingBoxGridProperties()      // explain how this thing works
   {
      loadingLatticeSitesDistanceHorizontal = 2.0*particleRadius*(1.0 + sizeDispersityParticle) + interParticleLoadingDistance;
      loadingLatticeSitesDistanceVertical = 0.5*sqrt(3.0)*loadingLatticeSitesDistanceHorizontal;

      loadingLatticeGridSize.X = (int)((fillingBoxMax.X - fillingBoxMin.X - interParticleLoadingDistance)/loadingLatticeSitesDistanceHorizontal);
      loadingLatticeGridSize.Z = (int)((fillingBoxMax.Z - fillingBoxMin.Z - interParticleLoadingDistance)/loadingLatticeSitesDistanceHorizontal);

      loadingLatticeOffset.X = 0.5*(fillingBoxMax.X - fillingBoxMin.X - (loadingLatticeSitesDistanceHorizontal*loadingLatticeGridSize.X + interParticleLoadingDistance));
      loadingLatticeOffset.Y = 0.0;
      loadingLatticeOffset.Z = 0.5*(fillingBoxMax.Z - fillingBoxMin.Z - (loadingLatticeSitesDistanceHorizontal*loadingLatticeGridSize.Z + interParticleLoadingDistance));

      int nParticlesPerOddLatticePlane = loadingLatticeGridSize.X*loadingLatticeGridSize.Z;
      int nParticlesPerEvenLatticePlane = (loadingLatticeGridSize.X - 1)*(loadingLatticeGridSize.Z - 1);

      loadingLatticeGridSize.Y = (int)(numberOfParticles/(nParticlesPerEvenLatticePlane + nParticlesPerOddLatticePlane));
      if (loadingLatticeGridSize.Y*(nParticlesPerEvenLatticePlane + nParticlesPerOddLatticePlane) - numberOfParticles > nParticlesPerOddLatticePlane)
      {
         loadingLatticeGridSize.Y = 2*(loadingLatticeGridSize.Y + 1);
      }
      else
      {
         loadingLatticeGridSize.Y = 2*loadingLatticeGridSize.Y + 1;
      }
   }

   void setLoadingBoxDimensions()
   {
      loadingBoxMin.X = fillingBoxMin.X;
      loadingBoxMin.Y = fillingBoxMax.Y;
      loadingBoxMin.Z = fillingBoxMin.Z;

      loadingBoxMax.X = fillingBoxMax.X;
      loadingBoxMax.Y = fillingBoxMax.Y + loadingLatticeSitesDistanceVertical*(loadingLatticeGridSize.Y + 1);
      loadingBoxMax.Z = fillingBoxMax.Z;
   }

   void setBoundaries()
   {
      setXMin(fillingBoxMin.X);
      setYMin(fillingBoxMin.Y);
      setZMin(fillingBoxMin.Z);

      setXMax(loadingBoxMax.X);
      setYMax(loadingBoxMax.Y);
      setZMax(loadingBoxMax.Z);
   }

   void makeBoundaries()
   {
      wallHandler.clear();
      InfiniteWall wall;
      PeriodicBoundary boundary;

      wall.setSpecies(speciesWall);

      wall.set(Vec3D(-1.,0.,0.),Vec3D(getXMin(),0.,0.));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(1.,0.,0.),Vec3D(getXMax(),0.,0.));
      wallHandler.copyAndAddObject(wall);

      wall.set(Vec3D(0,-1.,0.),Vec3D(0.,getYMin(),0.));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(0.,1.,0.),Vec3D(0.,getYMax(),0.));
      wallHandler.copyAndAddObject(wall);

      boundary.set(Vec3D(0.,0.,1.), getZMin(), getZMax());
      boundaryHandler.copyAndAddObject(boundary);
   }

   void makeScrew()
   {
      ScrewAuger screw;

      screw.setSpecies(speciesWall);
      screw.set(screwOrigin, screwLength, screwBladeRadius, screwShaftRadius, numberOfScrewPitches, screwAngularVelocity[screwVelocityIncrementControlVariable], screwThickness, true);
      screwPointer = wallHandler.copyAndAddObject(screw);

      std::cout << "Angular velocity array: ";
      for (int i = 0; i < nVelocityIncrements; i++) std::cout << screwAngularVelocity[i] << " ";
      std::cout << std::endl;
   }

   void makeParticles()
   {
      Vec3D particlePosition;
      BaseParticle p0;

      p0.setSpecies(speciesParticle);
      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));

      if (particleHackFlag) loadingLatticeGridSize.Y = particleHackLayers;

      for (int i = 0; i < loadingLatticeGridSize.Y; i++)
      {
         particlePosition.Y = loadingBoxMin.Y + loadingLatticeOffset.Y + i*loadingLatticeSitesDistanceVertical;

         if (!(i%2))
         {
            for (int j = 0; j < loadingLatticeGridSize.X; j++)
            {
               particlePosition.X = loadingBoxMin.X + loadingLatticeOffset.X + particleRadius*(1.0 + sizeDispersityParticle) + interParticleLoadingDistance + j*loadingLatticeSitesDistanceHorizontal;
               for (int k = 0; k < loadingLatticeGridSize.Z; k++)
               {
                  particlePosition.Z = loadingBoxMin.Z + loadingLatticeOffset.Z + particleRadius*(1.0 + sizeDispersityParticle) + interParticleLoadingDistance + k*loadingLatticeSitesDistanceHorizontal;

                  p0.setPosition(particlePosition);
                  p0.setRadius(particleRadius*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.0,1.0)));
                  particleHandler.copyAndAddObject(p0);
               }
            }
         }
         else
         {
            for (int j = 0; j < loadingLatticeGridSize.X - 1; j++)
            {
               particlePosition.X = loadingBoxMin.X + loadingLatticeOffset.X + 0.5*interParticleLoadingDistance + (j+1)*loadingLatticeSitesDistanceHorizontal;
               for (int k = 0; k < loadingLatticeGridSize.Z - 1; k++)
               {
                  particlePosition.Z = loadingBoxMin.Z + loadingLatticeOffset.Z + 0.5*interParticleLoadingDistance + (k+1)*loadingLatticeSitesDistanceHorizontal;

                  p0.setPosition(particlePosition);
                  p0.setRadius(particleRadius*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.,1.)));
                  particleHandler.copyAndAddObject(p0);
               }
            }
         }
      }
   }

   void makeCasing()
   {
      std::cout << "Cutting particles..." << std::endl;
      for (int i = particleHandler.getNumberOfObjects() - 1; i >= 0; i--)
      {
         if (pow(particleHandler.getObject(i) -> getPosition().X - screwOrigin.X,2.) +
         pow(particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y,2.) >
         pow(screwCasingRadius - 0.95*(particleHandler.getObject(i) -> getRadius()),2.)) particleHandler.removeObject(i);
      }

      numberOfParticles = particleHandler.getNumberOfObjects();
      std::cout << "New number of particles: " << numberOfParticles << std::endl;

      AxisymmetricIntersectionOfWalls casing;

      casing.setSpecies(speciesWall);
      casing.setPosition(screwOrigin);
      casing.setOrientation(Vec3D(0.0,0.0,1.0));
      casing.addObject(Vec3D(1.0,0.0,0.0),Vec3D(screwCasingRadius,0.0,0.0));
      casing.setAngularVelocity(Vec3D(0.,0.,0.));
      casingPointer = wallHandler.copyAndAddObject(casing);
   }

   void makeFillingRatio()
   {
      std::cout << "Cutting particles..." << std::endl;
      // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      for (int i = particleHandler.getNumberOfObjects() - 1; i >= 0; i--)
      {
         if (particleHandler.getObject(i) -> getPosition().Y - screwOrigin.Y > (2.0*screwFillingRatio-1.0)*screwCasingRadius) particleHandler.removeObject(i);
      }

      numberOfParticles = particleHandler.getNumberOfObjects();
      std::cout << "New number of particles: " << numberOfParticles << std::endl;
   }

   void makeParticleList()
   {
      particlesAxialPositionList.clear();
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++) particlesAxialPositionList.push_back(particleHandler.getObject(i) -> getPosition().Z);
   }

   void makeOutputFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "1.time\t 2.theta\t 3.omega\t 4.Rmean\t 5.PhiMean\t 6.zMean\t 7.vRmean\t 8.vPhiMean\t 9.vZmean\t 10.dV/dt\t 11.torque\t 12.meanOverlap\t 14.maxOverlap" << std::endl;
   }

   void makeRotation()
   {
      screwPointer -> rotate(getTimeStep());
      screwPointer -> setAngularVelocity(Vec3D(0.,0.,-screwAngularVelocity[screwVelocityIncrementControlVariable]));
      screwPointer -> setOrientation(Vec3D(0.0,0.0,1.0));

      if (screwVelocityIncrementControlVariable < nVelocityIncrements && getTime() > runtimePerScrewVelocity*(screwVelocityIncrementControlVariable + 1) + t0)
      {
         screwVelocityIncrementControlVariable++;
         screwPointer -> setOmega(screwAngularVelocity[screwVelocityIncrementControlVariable]);
         std::cout << std::endl << "NEW SCREW VELOCITY: " << screwAngularVelocity[screwVelocityIncrementControlVariable] << std::endl << std::endl;
      }
   }

   void makeDataAnalysis()
   {
      double particleRadialPosition;
      double particleAngularPosition;

      meanParticleCylindricalPosition.setZero();
      meanParticleCylindricalVelocity.setZero();
      volumetricThroughput = 0.0;

      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         particleRadialPosition = sqrt(pow(particleHandler.getObject(i) -> getPosition().X,2.0) + pow(particleHandler.getObject(i) -> getPosition().Y,2.0));
         particleAngularPosition = atan2(particleHandler.getObject(i) -> getPosition().Y, particleHandler.getObject(i) -> getPosition().X);
         if (particleAngularPosition < 0.0) particleAngularPosition += 2.0*constants::pi;

         meanParticleCylindricalPosition.X += particleRadialPosition;
         meanParticleCylindricalPosition.Y += particleAngularPosition;
         meanParticleCylindricalPosition.Z += particleHandler.getObject(i) -> getPosition().Z;

         meanParticleCylindricalVelocity.X += ((particleHandler.getObject(i) -> getPosition().X)*(particleHandler.getObject(i) -> getVelocity().X) + (particleHandler.getObject(i) -> getPosition().Y)*(particleHandler.getObject(i) -> getVelocity().Y))/particleRadialPosition;
         meanParticleCylindricalVelocity.Y += ((particleHandler.getObject(i) -> getPosition().Y)*(particleHandler.getObject(i) -> getVelocity().X) - (particleHandler.getObject(i) -> getPosition().X)*(particleHandler.getObject(i) -> getVelocity().Y))/pow(particleRadialPosition,2.0);
         meanParticleCylindricalVelocity.Z += particleHandler.getObject(i) -> getVelocity().Z;

         // volumetric throughput
         if (fabs(particleHandler.getObject(i) -> getPosition().Z - particlesAxialPositionList[i]) > 0.5*screwLength)
         {
            (particleHandler.getObject(i) -> getPosition().Z - particlesAxialPositionList[i]) < 0.0 ?
               (volumetricThroughput += pow(particleHandler.getObject(i) -> getRadius(),3.0)) :
               (volumetricThroughput -= pow(particleHandler.getObject(i) -> getRadius(),3.0));
         }
      }

      meanParticleCylindricalPosition /= particleHandler.getNumberOfObjects();
      meanParticleCylindricalVelocity /= particleHandler.getNumberOfObjects();
      volumetricThroughput *= 4.0*constants::pi/(3.0*cdatOutputTimeInterval);

      double particleRelativeOverlap;
      meanRelativeOverlap = 0.0;
      maxRelativeOverlap = 0.0;
      totalTorque = 0.0;
      int interactionCounter = 0;

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         particleRelativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius());
         meanRelativeOverlap += particleRelativeOverlap;
         if (particleRelativeOverlap > maxRelativeOverlap) maxRelativeOverlap = particleRelativeOverlap;

         // torque computation
         if ((*i) -> getI() -> getIndex() == screwPointer -> getIndex()) totalTorque -= ((*i) -> getContactPoint()).X * ((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y * ((*i) -> getForce()).X;

         interactionCounter++;
      }

      meanRelativeOverlap /= interactionCounter;
   }

   void writeDataToOutptFile()
   {
      cdatFile <<
      getTime() << "   " <<
      screwPointer -> getAngularOffset() << "   " <<
      screwAngularVelocity[screwVelocityIncrementControlVariable] << "   " <<
      meanParticleCylindricalPosition.X << "   " <<
      meanParticleCylindricalPosition.Y << "   " <<
      meanParticleCylindricalPosition.Z << "   " <<
      meanParticleCylindricalVelocity.X << "   " <<
      meanParticleCylindricalVelocity.Y << "   " <<
      meanParticleCylindricalVelocity.Z << "   " <<
      volumetricThroughput << "   " <<
      totalTorque << "   " <<
      meanRelativeOverlap << "   " <<
      maxRelativeOverlap << "   " <<
      std::endl;
   }


   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      if (stage < 4)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() <<
         ", Eratio = " << std::setprecision(6) << std::left << std::setw(10) << getKineticEnergy()/getElasticEnergy() << std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() <<
         ", theta = " << screwPointer -> getAngularOffset() << ", omega = " << screwAngularVelocity[screwVelocityIncrementControlVariable] <<
         ", vCyl = " << meanParticleCylindricalVelocity << ", dV/dt = " << volumetricThroughput << ", torque = " << totalTorque << std::endl;
      }
   }


   // CLASS VARIABLES ----------------------------------------
   // particle properties
   int numberOfParticles;
   double particleRadius;
   double densityParticles;
   double sizeDispersityParticle;
   double particleMass;
   double particleVolume;
   double totalParticleVolume;

   // interaction properties
   double densityWalls;
   double particleParticleSlidingFriction;
   double particleParticleRollingFriction;
   double particleParticleTorsionFriction;
   double particleWallSlidingFriction;
   double particleWallRollingFriction;
   double particleWallTorsionFriction;
   double particleStiffness;
   double wallsStiffness;
   double particleParticleRestitutionCoefficient;
   double particleWallRestitutionCoefficient;

   // screw geometry parameters
   double rescaledCasingRadius;
   double rescaledBladeRadius;
   double rescaledShaftRadius;
   double rescaledScrewThickness;
   double rescaledPitchLength;
   double screwCasingRadius;
   double screwBladeRadius;
   double screwShaftRadius;
   double screwThickness;
   double screwLength;
   int numberOfScrewPitches;
   Vec3D screwOrigin;
   AxisymmetricIntersectionOfWalls *casingPointer;
   ScrewAuger *screwPointer;

   // screw operational parameters
   double runtimePerScrewVelocity;
   int nVelocityIncrements;
   std::vector<double> screwAngularVelocity;
   double screwFillingRatio;
   int screwVelocityIncrementControlVariable;

   // initial loading parameters
   double assumedInitialPackingFraction;
   double settlingBoxToScrewSizeRatio;
   double interParticleLoadingLatticeRelativeDistance;
   double interParticleLoadingDistance;
   Vec3D fillingBoxMin;
   Vec3D fillingBoxMax;
   double fillingBoxVolume;
   Vec3D loadingLatticeOffset;
   Vec3D loadingLatticeGridSize;
   double loadingLatticeSitesDistanceHorizontal;
   double loadingLatticeSitesDistanceVertical;
   Vec3D loadingBoxMin;
   Vec3D loadingBoxMax;

   // species
   LinearViscoelasticFrictionSpecies *speciesParticle;
   LinearViscoelasticFrictionSpecies *speciesWall;
   LinearViscoelasticFrictionSpecies *speciesMixedParticleWall;

   // data analysis
   std::vector<double> particlesAxialPositionList;
   Vec3D meanParticleCylindricalPosition;
   Vec3D meanParticleCylindricalVelocity;
   double volumetricThroughput;
   double totalTorque;
   double meanRelativeOverlap;
   double maxRelativeOverlap;

   // global parameters
   std::ofstream cdatFile;
   double cdatOutputTimeInterval;
   double energyRatioTolerance;
   double timeBuffer;
   double t0;
   int stage;
   bool particleHackFlag;
   int particleHackLayers;
};


int main(int argc, char *argv[])
{
   // TIME STEP
   double timeStep = 2.0e-6;

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;
   double sizeDispersionParticles = 0.0;
   double densityParticles = 1500.0;
   double densityWall = 3000.0;

   // INTERACTION PARAMETERS
   double particleStiffness = 1000.0;
   double wallsStiffness = 2000.0;
   double particleParticleRestitutionCoefficient = 0.5;
   double particleWallRestitutionCoefficient = 0.8;

   // FRICTION PARAMETERS
   double muSlidingParticleParticle = 0.5;
   double muRollingParticleParticle = 0.05;
   double muSlidingParticleWall = 0.3;
   double muRollingParticleWall = 0.05;

   // SCREW GEOMETRY PARAMETERS
   double rescaledCasingRadius = 22.0;
   double rescaledBladeRadius = 17.0;
   double rescaledShaftRadius = 5.0;
   double rescaledScrewThickness = 2.0;
   double rescaledPitchLength = 44.0;
   int numberOfScrewPitches = 2;
   Vec3D screwOrigin = {0.0, 0.0, 0.0};

   // SCREW OPERATIONAL PARAMETERS
   double runtimePerScrewVelocity = 10.0;
   const int nVelocityIncrements = 3;
   double screwAngularVelocity[nVelocityIncrements] = {0.5*constants::pi, constants::pi, 2.0*constants::pi};
   const std::vector<double> screwFillingRatio = {0.25, 0.5, 0.75}; // loop index i

   // SYSTEM SETUP PARAMETERS
   double assumedInitialPackingFraction = 0.7;
   double settlingBoxToScrewSizeRatio = 1.05;
   double interParticleLoadingLatticeRelativeDistance = 0.1;

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 1.0e-4;
   double timeBuffer = 0.01;

   // EXECUTION LOOP
   for (int i = 0; i < screwFillingRatio.size(); i++)
   {
      // TIME COMPUTATION SETUP
      std::clock_t clockStart;
      clockStart = std::clock();
      auto chronoStart = std::chrono::high_resolution_clock::now();

      // INITIALIZATION
      ScrewAlgorithmComparison problem;
      problem.setXBallsAdditionalArguments("-v0 -h 800 -p 10 -o 200 -3dturn 1");

      problem.setTimeStep(timeStep);
      problem.setTimeMax(runtimePerScrewVelocity*nVelocityIncrements);
      problem.setGravity(Vec3D(0.00,-9.81,0.00));
      problem.setSystemDimensions(3);
      problem.setCdatOutputTimeInterval(0.001);
      problem.setEnergyRatioTolerance(energyRatioTolerance);
      problem.setTimeBuffer(timeBuffer);
      problem.setSaveCount(0.001/problem.getTimeStep());

      problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
      problem.setWallDensity(densityWall);
      problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
      problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
      problem.setInteractionProperties(particleStiffness, wallsStiffness, particleParticleRestitutionCoefficient, particleWallRestitutionCoefficient);

      problem.setScrewRelativeGeometryParameters(screwOrigin, rescaledCasingRadius, rescaledBladeRadius, rescaledShaftRadius, rescaledScrewThickness, rescaledPitchLength, numberOfScrewPitches);
      problem.setScrewOperationalParameters(runtimePerScrewVelocity, nVelocityIncrements, screwAngularVelocity, screwFillingRatio[i]);

      problem.setSetupParameters(assumedInitialPackingFraction, settlingBoxToScrewSizeRatio, interParticleLoadingLatticeRelativeDistance);
      problem.developmentHackParticles(true, 5);

      // NAME SETTING
      std::ostringstream name;
      name.str("");
      name.clear();
      std::cout.unsetf(std::ios::floatfield);
      name << "ScrewAlgorithmComparison_analytical_rP_" << particleRadius << "_rD_" << sizeDispersionParticles
      << "___" << rescaledCasingRadius << "_" << rescaledBladeRadius << "_" << rescaledShaftRadius << "_" << rescaledPitchLength << "_" << numberOfScrewPitches
      << "_fR_" << screwFillingRatio[i];
      problem.setName(name.str());

      problem.solve();

      // COMPUTATION OF EXECUTION TIME
      std::ofstream runningTime;
      std::ostringstream runningTimeFileName;
      std::cout.unsetf(std::ios::floatfield);
      runningTimeFileName << problem.getName() << ".time";

      runningTime.open(runningTimeFileName.str(), std::ios::out);
      runningTime << "CPU time used: " << (std::clock() - clockStart) / CLOCKS_PER_SEC << "s" << std::endl;
      runningTime << "Wall clock time passed: " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - chronoStart).count() << "s" << std::endl;

      runningTime.close();
   }

   return 0;
}
