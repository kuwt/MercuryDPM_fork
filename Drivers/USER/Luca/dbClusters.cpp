/*
 *** dbClusters - DEFORMABLE BREAKABLE CLUSTERS ***
 1) A single cluster is created and its properties computed
 2) Gravity is activated and the cluster left to settle
 3) The cluster is compressed uniaxially up to a desired force
 4) The piston unloaded and the simulation terminated

 LAST UPDATE: 2.10.18

  ToDo:
  - implement verbose flag
  - implement cluster PSD
  - standardise and re-add the cluster shape factor
  - fix the internal structure analysis for non-spherical clusters
  - DONE *** implement detailed cluster bounds analysis and output
  - DONE *** implement flags for optional data analysis and print (e.g. the .grid file)
  - re-add the multiple cluster loading lattice (and rename some variables accordingly)
  - add a separate phi for loose particles
  - implement user defined parameters check function
  - get rid of the matrix implementation for something faster
  - DONE *** reorder the functions alphabetically
  - check for "..." outputs and newlines, and standardize the output
  */

#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <fstream>

class dbClusters_compression : public Mercury3D
{
private: //  ----- SIMULATION FLOW ---------------------------------------------
   void setupInitialConditions() override
   {
      // resetting the counters
      stage = 0;
      nParticlesInserted = 0;
      numberOfIntraClusterBonds = 0;
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

      std::cout << "SETTING INITIAL CLUSTER LATTICE SITES... " << std::endl;
      setInitialClusterLatticeSites();

      // std::cout << "CREATING BOUNDARIES..." << std::endl;
      // makeBoundaries();

      std::cout << std::endl << "PARTICLE INSERTION" << std::endl;
      insertParticles();

      std::cout << std::endl << "COMPUTING PARTICLES VOLUME" << std::endl;
      computeTotalParticleVolume();

      if (isCdatOutputON)
      {
         std::cout << "CREATING .cdat FILE" << std::endl;
         makeCdatFile();
      }

      std::cout << "ACTIVATING CENTRAL FORCES" << std::endl;

      stage++;
   }

   void actionsOnRestart() override
   {
      std::cout << "NO RESTART BRANCH IMPLEMENTED YET" << std::endl;
      std::cout << "QUITTING..." << std::endl;

      exit(1);
   }

   void actionsAfterTimeStep() override
   {
      // DATA ANALYSIS
      // at every stage make data analysis at each time step
      if (stage >= 1) makeDataAnalysis();

      // at every stage print the results from the data analysis every specific time interval
      if (stage >= 1 && isCdatOutputON && fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) writeToCdatFile();

      // ACTIVATION OF LOCAL FORCE
      if (stage == 1)
      {
         // applies the central force to each particle
         applyCentralForce();

         // force increase and velocity dampening
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
         // applies the central force to each particle
         applyCentralForce();

         // force decrease and velocity dampening
         if (getTime() - t0 < gravityAndForceTuningDuration)
         {
            if (fmod(getTime() - t0,forceTuningInterval) < getTimeStep()) dampForce();
            if (fmod(getTime() - t0,velocityDampingInterval) < getTimeStep()) dampVelocities();
         }
         else
         {
            // the adjacency matrix is created
            std::cout << "CREATING ADJACENCY MATRIX" << std::endl << std::endl;
            setupAdjacencyMatrices();

            std::cout << "DISSIPATING ENERGY" << std::endl << std::endl;
            t0 = getTime();
            stage++;
         }
      }

      // ENERGY DISSIPATION AND CREATION OF THE BASE
      if (stage == 3)
      {
         // wait till the system is static
         if (getKineticEnergy()/getElasticEnergy() < energyRatioTolerance || getTime() - t0 > gravityAndForceTuningDuration)
         {
            std::cout << "ENERGY DISSIPATED" << std::endl << std::endl;

            // computes the internal structure of the cluster by setting up the grid
            // *** VERY COMPUTATIONALLY INTENSIVE! ***
            if (isGridOutputON)
            {
               std::cout << "COMPUTING INTERNAL STRUCTURE" << std::endl << std::endl;
               computeInternalStructure();
            }

            // generates the adjacency matrix output file
            if (isAmatOutputON)
            {
               std::cout << "CREATING .amat FILE" << std::endl << std::endl;
               makeAmatFile();
            }

            // creates a cylindrical casing for the compression test
            std::cout << "CREATING BASE AND HOUSING" << std::endl << std::endl;
            makeBaseAndHousing();

            std::cout << "ACTIVATING GRAVITY" << std::endl << std::endl;
            t0 = getTime();
            stage++;

            // // prints the adjacency matrix
            // std::ofstream matrixOutput;
            // matrixOutput.open("Clusters_ADJACENCYMATRIX_BEGIN", std::ios::out);
            // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            // {
            //    for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) matrixOutput << adjacencyMatrix[i][j] << "\t";
            //    matrixOutput << std::endl;
            // }
            // matrixOutput.close();
         }
      }

      // from now on the structure of the cluster might change, so refresh the adjacency matrix every specific time interval
      if (stage >= 4 && fmod(getTime(), cdatOutputTimeInterval) < getTimeStep()) refreshAdjacencyAndCohesionMatrices();

      // ACTIVATING GRAVITY, SETTLING AND CREATION OF THE PISTON
      if (stage == 4)
      {
         // slowly increases gravity to the real (Earth) value
         if (getTime() - t0 < gravityAndForceTuningDuration)
         {
            activateGravity();
         }

         // if the systerm is static then start with the compression test
         if (getTime() - t0 > gravityAndForceTuningDuration && getKineticEnergy()/getElasticEnergy() < energyRatioTolerance)
         {
            // creates the compression piston directly above the particle bed
            std::cout << "CREATING PISTON AND STARTING COMPRESSION" << std::endl << std::endl;
            makeCompressionPiston();

            stage++;
         }
      }

      // COMPRESSION CYCLE
      if (stage == 5)
      {
         // compresses until the maximum force is met
         if (pistonForce < maximumCompressiveForce)
         {
            loadPiston();
         }
         else
         {
            std::cout << "MAXIMUM COMPRESSION FORCE REACHED. UNLOADING PISTON" << std::endl << std::endl;
            stage++;
         }
      }

      // UNLOAD CYCLE AND SHUTDOWN
      if (stage == 6)
      {
         // decompresses until the measured force against the piston is almost zero
         if (pistonForce > energyRatioTolerance)
         {
            unloadPiston();
         }
         else
         {
            // sets the time max to the next time step, therefore terminating the simulation at the next time step
            std::cout << "PISTON UNLOADED. QUITTING" << std::endl << std::endl;
            setTimeMax(getTime() + getTimeStep());
            stage++;

            // // prints the adjacency matrix
            // std::ofstream matrixOutput;
            // matrixOutput.open("Clusters_ADJACENCYMATRIX_END", std::ios::out);
            // for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            // {
            //    for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) matrixOutput << adjacencyMatrix[i][j] << "\t";
            //    matrixOutput << std::endl;
            // }
            // matrixOutput.close();
         }
      }
   }

   void actionsAfterSolve() override
   {
      // close every stream and de-allocate every dynamic pointer
      if (isCdatOutputON) cdatFile.close();
      if (isAmatOutputON) amatFile.close();
      delete [] boxCentres;
   }


public: //  ----- FUNCTIONS AND VARIABLES --------------------------------------
   //  ----- FUNCTIONS CALLED IN MAIN ------------------------------------------
   // sets the output of additional data files ON/OFF
   void setAdditionalDataOutput(bool doCdat, bool doGrid, bool doAmat)
   {
      isCdatOutputON = doCdat;
      isGridOutputON = doGrid;
      isAmatOutputON = doAmat;
   }

   // sets the additional data analysis output time interval
   void setCdatOutputTimeInterval(double dt)
   {
      cdatOutputTimeInterval = dt;
   }

   // sets the cluster main properties
   void setClustersProperties(int nP, double sizeRatio, int nX, int nY, int nZ)
   {
      numberOfParticlesPerCluster = nP;
      clusterSizeToLatticeRatio = sizeRatio;
      nXLatticeSites = nX;
      nYLatticeSites = nY;
      nZLatticeSites = nZ;
      numberOfClusters = nX*nY*nZ;
      totalNumberOfParticles = numberOfParticlesPerCluster*numberOfClusters;
   }

   // sets the compression parameters (right now only the maximum compressive force to be reached)
   void setCompressionParameters(double maxForce, double sizeRatio, double compressionRatio)
   {
      maximumCompressiveForce = maxForce;
      cylinderToClusterSizeRatio = sizeRatio;
      particleRadiusToPistonDisplacementPerStepRatio = compressionRatio;
   }

   // sets the energy ratio tolerance (for settling and compression)
   void setEnergyRatioTolerance(double eRatio)
   {
      energyRatioTolerance = eRatio;
   }

   // sets the agglomeration force properties and the velocity damping parameters
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

   // sets the time interval for force and gravity tuning
   void setGravityAndForceTuningDuration(double dt)
   {
      gravityAndForceTuningDuration = dt;
   }

   // sets some parameters for the cluster internal structure analysis
   void setInternalStructureAnalysisParameters(int length, double ratio)
   {
      gridLength = length;
      fictiousGridPointRadiusRatio = ratio;
   }

   // sets the elasto-plastic particle-particle interaction parameters
   void setParticleElastoPlasticProperties(double kp, double ke, double kcCluster, double kcLoose, double phi)
   {
      kpParticle = kp;
      keParticle = ke;
      kcParticleIntercluster = kcCluster;
      kcParticleLoose = kcLoose;
      phiParticle = phi;
   }

   // sets the particle-particle friction coefficients
   void setParticleParticleFrictionCoefficients(double muSliding, double muRolling, double muTorsion)
   {
      particleParticleSlidingFrictionCoeff = muSliding;
      particleParticleRollingFrictionCoeff = muRolling;
      particleParticleTorsionFrictionCoeff = muTorsion;
   }

   // sets the particle-particle restitution coefficient
   void setParticleParticleRestitutionCoefficients(double bigBigE)
   {
      particleParticleRestitutionCoeff = bigBigE;
   }

   // sets the elementary particles main properties
   void setParticleProperties(double rP, double dP, double rhoP)
   {
      radiusParticle = rP;
      sizeDispersityParticle = dP;
      densityParticle = rhoP;

      // computes particle average volume and mass
      volumeParticle = 4.0*constants::pi*pow(radiusParticle,3.0)/3.0;
      massParticle = volumeParticle*densityParticle;
   }

   // sets the particle-wall friction coefficients
   void setParticleWallFrictionCoefficients(double muSliding, double muRolling, double muTorsion)
   {
      particleWallSlidingFrictionCoeff = muSliding;
      particleWallRollingFrictionCoeff = muRolling;
      particleWallTorsionFrictionCoeff = muTorsion;
   }

   // sets the walls density
   void setWallDensity(double rhoW)
   {
      densityWall = rhoW;
   }

   // sets the wall stiffness and restitution coefficients
   void setWallStiffnessAndRestitutionCoefficients(double kp, double eB)
   {
      kpWall = kp;
      particleWallRestitutionCoeff = eB;
   }


   //  ----- FUNCTIONS CALLED IN THE CLASS -------------------------------------
   // tunes gravity continuously from 0 to g
   void activateGravity()
   {
      setGravity(Vec3D(0.00,0.00,-9.81*(getTime() - t0)/gravityAndForceTuningDuration));
   }

   // applies the central force to every particles of the cluster
   void applyCentralForce()
   {
      Vec3D distanceFromForceCenter;

      // count the particles belonging to every single cluster and point them towards their box centre
      for (int j=0; j<numberOfClusters; j++)
      {
         // computes the distance of each particle from the force center and applies a proportional force
         for (int i=0; i<numberOfParticlesPerCluster; i++)
         {
            distanceFromForceCenter = particleHandler.getObject(i + numberOfParticlesPerCluster*j) -> getPosition() - boxCentres[j];
            particleHandler.getObject(i + numberOfParticlesPerCluster*j) -> addForce(-forceModulus*distanceFromForceCenter);
         }
      }
   }

   // computes the internal structure of the cluster
   void computeInternalStructure()
   {
      // initialization of important parameters
      Vec3D gridPoint;
      SphericalParticle p0;
      double nPointsInsideAgglomerateBoundary;
      double nPointsInsideComponents;

      gridResolutionLength = 2.0*meanClusterRadius/(gridLength - 1.0);
      nPointsInsideAgglomerateBoundary = 0;
      nPointsInsideComponents = 0;

      // creates the output file
      makeGridFile();

      // setup of the fictious particles for structure computation
      p0.setRadius(radiusParticle*fictiousGridPointRadiusRatio);
      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setSpecies(speciesWall);

      // creates fictious particles at every grid node
      // if they are inside the cluster perimeter then save them in the .grid file
      // if they are inside a particle label with a "1", with a "0" otherwise
      for(int i = 0; i < gridLength; i++)
      {
         gridPoint.X = centerOfMass.X - meanClusterRadius + gridResolutionLength*i;

         for(int j = 0; j < gridLength; j++)
         {
            gridPoint.Y = centerOfMass.Y - meanClusterRadius + gridResolutionLength*j;

            for(int k = 0; k < gridLength; k++)
            {
               gridPoint.Z = centerOfMass.Z - meanClusterRadius + gridResolutionLength*k;

               // sets the position of the fictious particle to the grid point
               p0.setPosition(gridPoint);

               // checks if inside the agglomerate boundary
               if ((gridPoint - centerOfMass).getLength() < meanClusterRadius)
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

      // computes the mass fraction by taking the ratio between nodes inside particles and nodes inside the perimeter
      massFraction = nPointsInsideComponents/nPointsInsideAgglomerateBoundary;

      // writes the fractions at the bottom of the grid file
      gridFile << "n_points_inside_boundary: " << nPointsInsideAgglomerateBoundary << std::endl;
      gridFile << "n_points_inside_components: " << nPointsInsideComponents << std::endl;
      gridFile << "mass_fraction: " << massFraction << std::endl;
      gridFile.close();
   }

   // computes the volume of all the particles in the cluster
   void computeTotalParticleVolume()
   {
      totalParticleVolume = 0.0;
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--) totalParticleVolume += particleHandler.getObject(i) -> getVolume();
   }

   // linearily decreases the central force modulus until a minimum value
   void dampForce()
   {
      forceModulus = maximumForceModulus - (maximumForceModulus - finalForceModulus)*(getTime() - t0)/gravityAndForceTuningDuration;
   }

   // damps the particles velocity by multiplying them by the user defined velocity damping modulus
   void dampVelocities()
   {
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         particleHandler.getObject(i) -> setVelocity(velocityDampingModulus*(particleHandler.getObject(i) -> getVelocity()));
      }
   }

   // linearily increases the central force modulus until a maximum value
   void increaseForce()
   {
      forceModulus = maximumForceModulus*(getTime() - t0)/gravityAndForceTuningDuration;
   }

   // inserts the particles in the simulation domain
   void insertParticles()
   {
      for (int i = 0; i < numberOfClusters; i++)
      {
         nParticlesInserted = 0;

         // insert particles till the desired amount is reached
         while(nParticlesInserted < numberOfParticlesPerCluster)
         {
            if (particleInsertionSuccessful(i)) nParticlesInserted++;
            else
            {
               std::cout << "CANNOT INSERT ALL THE PARTICLES" << std::endl;
               std::cout << "MAYBE TRY TO CHANGE THE clusterSizeToLatticeRatio VARIABLE"<< std::endl;
               std::cout << "QUITTING"<< std::endl;

               exit(1);
            }
         }
         std::cout << "Inserted cluster n. " << i + 1 << "/" << numberOfClusters << std::endl;
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
   }

   // loads the piston by moving it downward
   void loadPiston()
   {
      pistonHeight -= pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }

   // creates the adjacency matrix file output
   void makeAmatFile()
   {
      std::ostringstream amatName;
      std::cout.unsetf(std::ios::floatfield);
      amatName << getName() << ".amat";

      amatFile.open(amatName.str(), std::ios::out);
   }

   // creates a housing of cylindrical shape to settle the formed clusters in
   void makeBaseAndHousing()
   {
      // clearing every wall
      wallHandler.clear();

      // base of the casing
      base.setSpecies(speciesWall);
      base.set(Vec3D(0.,0.,-1.),Vec3D(0.,0.,particleHandler.getLowestPositionComponentParticle(2) -> getPosition().Z - 1.1*(1.0 + sizeDispersityParticle)*radiusParticle));
      basePointer = wallHandler.copyAndAddObject(base);

      // computes the maximum radial distance of the particles from the z axis
      double localRadialDistanceSquared;
      double maxRadialDistanceSquared = 0.0;

      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         localRadialDistanceSquared = pow(particleHandler.getObject(i) -> getPosition().X, 2.0) + pow(particleHandler.getObject(i) -> getPosition().Y, 2.0);
         if (localRadialDistanceSquared > maxRadialDistanceSquared) maxRadialDistanceSquared = localRadialDistanceSquared;
      }

      // cylindrical wall of radius R = cylinder_to_cluster_size_ratio*(mean_particle_radius + max_radial_distance)
      AxisymmetricIntersectionOfWalls(cylinder);
      cylinder.setSpecies(speciesWall);
      cylinder.setPosition(Vec3D(0.0,0.0,0.0));
      cylinder.setOrientation(Vec3D(0.0,0.0,1.0));
      cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(cylinderToClusterSizeRatio*(radiusParticle + pow(maxRadialDistanceSquared, 0.5)),0.0,0.0));
      wallHandler.copyAndAddObject(cylinder);
   }

   // creates the personalized data output file and writes the first row
   void makeCdatFile()
   {
      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "time \t Eel \t Ekin/Eel \t coord_number \t n_bounds \t meanRadius \t Vgran/Vtot \t Fcomp \t Fmean \t dTotMean \t dTotMax \t cm_X \t cm_Y \t cm_Z \t Fpiston \t h_piston" << std::endl;
   }

   // creates the piston to compress the clusters
   void makeCompressionPiston()
   {
      // sets piston initial height and piston displacement per time step
      pistonStepDisplacement = (1.0 - sizeDispersityParticle)*radiusParticle/particleRadiusToPistonDisplacementPerStepRatio;
      pistonHeight = particleHandler.getHighestPositionComponentParticle(2) -> getPosition().Z + 1.5*(1.0 + sizeDispersityParticle)*radiusParticle;

      // initializes the piston
      piston.setSpecies(speciesWall);
      piston.set(Vec3D(0.,0.,1.),Vec3D(0.,0.,pistonHeight));
      pistonPointer = wallHandler.copyAndAddObject(piston);
   }

   // makes the extended data analysis
   void makeDataAnalysis()
   {
      // resets counters and variables
      int totalInteractionCounter = 0;
      Vec3D distanceFromForceCenter;
      Vec3D localMin, localMax;
      localMin.setZero();
      localMax.setZero();
      centerOfMass.setZero();

      averageForceOnParticle = 0.0;
      meanClusterRadius = 0.0;
      meanCoordinationNumber = 0.0;
      volumeRatio = 0.0;
      maxTotalRelativeOverlap = 0.0;
      meanTotalRelativeOverlap = 0.0;
      pistonForce = 0.0;

      // loops over each particle to compute mean coordination number, center of mass, mean cluster radius and cluster volume to total particle volume ratio
      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
         centerOfMass += (particleHandler.getObject(i) -> getVolume())*(particleHandler.getObject(i) -> getPosition());

         distanceFromForceCenter = particleHandler.getObject(i) -> getPosition();
         if (distanceFromForceCenter.X > localMax.X) localMax.X = distanceFromForceCenter.X;
         if (distanceFromForceCenter.X < localMin.X) localMin.X = distanceFromForceCenter.X;
         if (distanceFromForceCenter.Y > localMax.Y) localMax.Y = distanceFromForceCenter.Y;
         if (distanceFromForceCenter.Y < localMin.Y) localMin.Y = distanceFromForceCenter.Y;
         if (distanceFromForceCenter.Z > localMax.Z) localMax.Z = distanceFromForceCenter.Z;
         if (distanceFromForceCenter.Z < localMin.Z) localMin.Z = distanceFromForceCenter.Z;
      }
      meanCoordinationNumber /= particleHandler.getNumberOfObjects();
      meanClusterRadius = (localMax.X - localMin.X + localMax.Y - localMin.Y + localMax.Z - localMin.Z)/6.0 + radiusParticle;
      volumeRatio = 4.0*constants::pi*pow(meanClusterRadius,3.0)/(3.0*totalParticleVolume);
      centerOfMass /= totalParticleVolume;

      // loops over every interaction to compute average force acting on particles, maximum and mean relative particle overlap and total piston applied force (during compression)
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

   // creates the grid file and writes the first row
   void makeGridFile()
   {
      std::ostringstream gridName;
      std::cout.unsetf(std::ios::floatfield);
      gridName << getName() << ".grid";

      gridFile.open(gridName.str(), std::ios::out);
      gridFile << "Grid_length: " << gridLength << "\t grid_resolution_length: " << gridResolutionLength << "\t fictiousGridParticleRadiusRatio: " << fictiousGridPointRadiusRatio << std::endl;
   }

   // tries to insert a particle in a spherical area centered around the initial lattice site of granule n; returns true if it succeeds, false otherwise
   bool particleInsertionSuccessful(int n)
   {
      // initialization of parameters
      int insertionFailCounter = 0;
      double rad, theta, phi;
      Vec3D particlePosition;
      SphericalParticle p0;

      // setup of particle properties and initial conditions (besides position)
      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setRadius(radiusParticle*(1.0 + sizeDispersityParticle*random.getRandomNumber(-1.0,1.0)));
      p0.setSpecies(particlePureSpeciesVector[particleHandler.getNumberOfObjects()]);

      // in this cycle a random position inside of a sphere contained in the bounding box is taken
      // if the particle is not in contact with the others then insert it, otherwise increment the fail counter
      // the maximum number of failed attempts is capped to 1000
      do
      {
         rad = random.getRandomNumber(p0.getRadius(),0.5*boxSize - 1.01*(p0.getRadius()));
         theta = constants::pi*random.getRandomNumber(-1.0,1.0);
         phi = 0.5*constants::pi*random.getRandomNumber(-1.0,1.0);
         particlePosition.X = boxCentres[n].X + rad*sin(theta)*cos(phi);
         particlePosition.Y = boxCentres[n].Y + rad*sin(theta)*sin(phi);
         particlePosition.Z = boxCentres[n].Z + rad*cos(theta);

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

   // refreshes the matrix containing the interaction information and updates the cohesion coefficients accordingly
   void refreshAdjacencyAndCohesionMatrices()
   {
      // resets the number of bounds between particles (after the cluster formation)
      numberOfIntraClusterBonds = 0;

      // updates the amat file if needed
      if (isAmatOutputON) amatFile << std::setprecision(3) << std::left << "# " << getTime() << std::endl;

      // resets the matrix content to 0
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         if (isAmatOutputON) amatFile << i << "\t";

         for (int j = 0; j < particleHandler.getNumberOfObjects(); j++)
         {
            if (isAmatOutputON && adjacencyMatrix[i][j] !=0) amatFile << j << "\t";

            adjacencyMatrix[i][j] = 0;
         }

         if (isAmatOutputON) amatFile << std::endl;
      }

      // re-allocates the contacts in the coordination matrix
      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         // the contact is preserved if it was in place at the previous check
         if (adjacencyMatrixPreviousCheck[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()])
         {
            adjacencyMatrix[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
            adjacencyMatrix[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;
            numberOfIntraClusterBonds++;
         }
      }

      // sets the old coordination matrix equal to the new one
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         for (int j = i + 1; j < particleHandler.getNumberOfObjects(); j++)
         {
            adjacencyMatrixPreviousCheck[i][j] = adjacencyMatrix[i][j];
            adjacencyMatrixPreviousCheck[j][i] = adjacencyMatrix[i][j];
         }
      }

      // updates the cohesion interactions
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         for (int j = i + 1; j < particleHandler.getNumberOfObjects(); j++)
         {
            if (adjacencyMatrix[i][j]) mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleIntercluster);
            else mixedSpeciesMatrix[i][j] -> setCohesionStiffness(kcParticleLoose);
         }
      }
   }

   // computes an overestimated radius of the cluster and the size of the box containing the cluster (for the initial cluster agglomeration)
   void setBoxSize()
   {
      // over-estimation of the cluster radius [V_c = total(V_i)/l -> R^3 ~ Np*rp^3/l -> R ~ rp*(Np/l)^(1/3), where l = underestimated packing fraction = 0.5 ]
      clusterRadiusPreFormation = radiusParticle*pow(totalNumberOfParticles/0.5,1.0/3.0);

      // the box size containing he cluster is rescaled according to the user defined clusterSizeToLatticeRatio
      boxSize = 2.0*clusterRadiusPreFormation/clusterSizeToLatticeRatio;

      std::cout << "Cluster average radius assuming packing fraction of 0.5 " << clusterRadiusPreFormation << std::endl;
      std::cout << "Cubic lattice size " << boxSize << std::endl << std::endl;
   }

   // sets the centers of the lattice sites where the clusters will be formed
   void setInitialClusterLatticeSites()
   {
      boxCentres = new Vec3D[nXLatticeSites*nYLatticeSites*nZLatticeSites];

      // offsets for the lattice sites placement
      double offsetX = -0.5*(nXLatticeSites - 1)*boxSize;
      double offsetY = -0.5*(nYLatticeSites - 1)*boxSize;
      double offsetZ = 0.0;
      // double offsetX = -0.5*nXLatticeSites*boxSize;
      // double offsetY = -0.5*nYLatticeSites*boxSize;
      // double offsetZ = -0.5*boxSize;

      for (int i = 0; i < nZLatticeSites; i++) // loop along Z
      {
         for (int j = 0; j < nYLatticeSites; j++) // loop along Y
         {
            for (int k = 0; k < nXLatticeSites; k++) // loop along X
            {
               boxCentres[k + j*nXLatticeSites + i*nXLatticeSites*nYLatticeSites].X = offsetX + k*boxSize;
               boxCentres[k + j*nXLatticeSites + i*nXLatticeSites*nYLatticeSites].Y = offsetY + j*boxSize;
               boxCentres[k + j*nXLatticeSites + i*nXLatticeSites*nYLatticeSites].Z = offsetZ + i*boxSize;
            }
         }
      }
   }

   // sets the domain boundaries of the simulation
   void setDomainLimits()
   {
      setXMin(-0.5*boxSize*nXLatticeSites);
      setYMin(-0.5*boxSize*nYLatticeSites);
      setZMin(-0.5*boxSize);

      setXMax(0.5*boxSize*nXLatticeSites);
      setYMax(0.5*boxSize*nYLatticeSites);
      setZMax((nZLatticeSites - 0.5)*boxSize);

      pistonHeight = getZMax();
   }

   // implements the species interaction matrix (for mixed interactions)
   void setMixedSpeciesCohesionMatrix()
   {
      // allocation of the species matrix
      for (int i = 0; i < totalNumberOfParticles; i++)
      {
         // temporary vector (the matrix is allocated by rows)
         std::vector<LinearPlasticViscoelasticFrictionMixedSpecies*> temporaryRowVector;

         for (int j = 0; j < totalNumberOfParticles; j++)
         {
            if (i == j) // on the main diagonal there are the mixed particle[i]-wall collision species (this is not necessary, the main diagonal could be left empty)
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
            else // on every [i,j] = [j,i] element of the matrix there is the mixed particle[i]-particle[j] collision species
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

         // allocation of the temporary vector in the matrix
         mixedSpeciesMatrix.push_back(temporaryRowVector);
      }

      // OUTPUTS
      // std::cout << "BIG-WALL stiffness and dissipation coefficients: " << mixedSpeciesMatrix[0][0] -> getLoadingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getDissipation() << "\n";
      // std::cout << "BIG-WALL friction coefficients: " << particleWallSlidingFrictionCoeff << " " << particleWallRollingFrictionCoeff << " " << particleWallTorsionFrictionCoeff << "\n";
      // std::cout << "BIG-WALL tangential stiffnesses: " << mixedSpeciesMatrix[0][0] -> getSlidingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getRollingStiffness() << " " << mixedSpeciesMatrix[0][0] -> getTorsionStiffness() << "\n";
      // std::cout << "BIG-WALL tangential dissipation coefficients: " << mixedSpeciesMatrix[0][0] -> getSlidingDissipation() << " " << mixedSpeciesMatrix[0][0] -> getRollingDissipation() << " " << mixedSpeciesMatrix[0][0] -> getTorsionDissipation() << "\n";
      // std::cout << "BIG-WALL collision time: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle) << "\n";

      // prints the ratio between collision time and time step for a mixed particle-wall interaction (should be > 50 to be super-precise)
      std::cout << "tC_BW/dt: " << std::setprecision(4) << mixedSpeciesMatrix[0][0] -> getCollisionTime(massParticle)/getTimeStep() << std::endl << std::endl;
   }

   // clears the handler and implements the species vector (for pure interactions)
   void setSpeciesVector()
   {
      speciesHandler.clear();

      // WALL-WALL
      // implementation of the wall species
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
      // every particle has its own species
      particlePureSpeciesVector.reserve(totalNumberOfParticles);

      // implementation of every pure particle species
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

      // OUTPUTS
      // std::cout << "BIG stiffness and dissipation coefficients: " << particlePureSpeciesVector[0] -> getLoadingStiffness() << " " << particlePureSpeciesVector[0] -> getDissipation() << "\n";
      // std::cout << "BIG friction coefficients: " << particleParticleSlidingFrictionCoeff << " " << particleParticleRollingFrictionCoeff << " " << particleParticleTorsionFrictionCoeff << "\n";
      // std::cout << "BIG tangential stiffnesses: " << particlePureSpeciesVector[0] -> getSlidingStiffness() << " " << particlePureSpeciesVector[0] -> getRollingStiffness() << " " << particlePureSpeciesVector[0] -> getTorsionStiffness() << "\n";
      // std::cout << "BIG tangential dissipation coefficients: " << particlePureSpeciesVector[0] -> getSlidingDissipation() << " " << particlePureSpeciesVector[0] -> getRollingDissipation() << " " << particlePureSpeciesVector[0] -> getTorsionDissipation() << "\n";
      // std::cout << "BIG collision time: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle) << "\n";

      // prints the ratio between collision time and time step for a pure particle-particle interaction (should be > 50 to be super-precise)
      std::cout << "tC_PP/dt: " << std::setprecision(4) << particlePureSpeciesVector[0] -> getCollisionTime(massParticle)/getTimeStep() << std::endl << std::endl;
   }

   // creates the matrices containing the interaction information: "1" -> in contact, "0" -> not in contact
   void setupAdjacencyMatrices()
   {
      // resets the number of bounds between particles (after the cluster formation)
      numberOfIntraClusterBonds = 0;

      // creates the matrix and fills it with zeroes
      for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
      {
         std::vector<int> temporaryRowVector;
         for (int j = 0; j < particleHandler.getNumberOfObjects(); j++) temporaryRowVector.push_back(0);
         adjacencyMatrix.push_back(temporaryRowVector);
         adjacencyMatrixPreviousCheck.push_back(temporaryRowVector);
      }

      // now checks for interactions and allocates a "1" to corresponding interacting particle indices
      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         adjacencyMatrix[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
         adjacencyMatrix[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;
         adjacencyMatrixPreviousCheck[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
         adjacencyMatrixPreviousCheck[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;

         // updates the number of bounds
         numberOfIntraClusterBonds++;
      }
   }

   // unloads the piston by moving it upwards
   void unloadPiston()
   {
      pistonHeight += pistonStepDisplacement;
      pistonPointer -> set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,pistonHeight));
   }

   // writes the personalized data output to the output file
   void writeToCdatFile()
   {
      cdatFile <<
      getTime() << "   " <<
      getElasticEnergy() << "   " <<
      getKineticEnergy()/getElasticEnergy() << "   " <<
      meanCoordinationNumber << "   " <<
      numberOfIntraClusterBonds << "   " <<
      meanClusterRadius << "   " <<
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


   //  ----- GLOBAL FUNCTIONS --------------------------------------------------
   // overrides the usual printTime() functions to display variables of interest according to the simulation stage
   void printTime() const override
   {
      if (stage < 5)
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(4) << getTime() << ", tmax = " << getTimeMax() <<
         ", E_ratio = " << getKineticEnergy()/getElasticEnergy() <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfIntraClusterBonds << ", rMean = " << meanClusterRadius << ", vF = " << volumeRatio <<
         ", Force Modulus = " << forceModulus << ", g_z = " << getGravity().Z << ", dMean = " << meanTotalRelativeOverlap << ", dMax = " << maxTotalRelativeOverlap <<
         std::endl;
      }
      else
      {
         std::cout << "t = " << std::setprecision(3) << std::left << std::setw(4) << getTime() << ", tmax = " << getTimeMax() <<
         ", E_ratio = " << getKineticEnergy()/getElasticEnergy() <<
         ", cN = " << meanCoordinationNumber << ", nB = " << numberOfIntraClusterBonds <<
         ", Fpiston = " << pistonForce << ", Hpiston = " << pistonHeight << ", dMean = " << meanTotalRelativeOverlap << ", dMax = " << maxTotalRelativeOverlap <<
         std::endl;
      }

      std::cout.flush();
   }


   //  ----- VARIABLES ---------------------------------------------------------
   // particles
   double radiusParticle;
   double sizeDispersityParticle;
   double densityParticle;
   double volumeParticle;
   double massParticle;
   double totalParticleVolume;
   int nParticlesInserted;

   // clusters
   int numberOfClusters;
   int numberOfParticlesPerCluster;
   int totalNumberOfParticles;
   double clusterRadiusPreFormation;
   double meanClusterRadius;
   int nXLatticeSites;
   int nYLatticeSites;
   int nZLatticeSites;
   Vec3D *boxCentres;
   Vec3D centerOfMass;

   // geometry
   double boxSize;
   double clusterSizeToLatticeRatio;
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
   int numberOfIntraClusterBonds;
   std::vector< std::vector<int> > adjacencyMatrix;
   std::vector< std::vector<int> > adjacencyMatrixPreviousCheck;

   // species
   LinearPlasticViscoelasticFrictionSpecies *speciesParticle, *speciesWall;
   std::vector<LinearPlasticViscoelasticFrictionSpecies*> particlePureSpeciesVector;
   std::vector< std::vector<LinearPlasticViscoelasticFrictionMixedSpecies*> > mixedSpeciesMatrix;

   // output
   double cdatOutputTimeInterval;
   std::ofstream cdatFile;
   std::ofstream amatFile;
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
   double maximumCompressiveForce;
   double cylinderToClusterSizeRatio;
   double particleRadiusToPistonDisplacementPerStepRatio;
   double pistonForce;
   double pistonStepDisplacement;
   double pistonHeight;

   // global
   int stage;
   double t0;
   double energyRatioTolerance;
   double gravityAndForceTuningDuration;
   bool isCdatOutputON;
   bool isAmatOutputON;
   bool isGridOutputON;
};


int main(int argc, char *argv[]) //  ----- MAIN --------------------------------
{
   // TIME PARAMETERS
   double timeStep = 1.0e-5;        // DEM solver integration time step
   double timeMax = 10.0;        // tMax at which the simulaton gets shut down
   double timeStepForDataOutput = 0.001;        // time between data prints

   // SETUP PARAMETERS
   double particleRadius = 5.0e-4;        // mean particle radius (uniform PSD)
   double sizeDispersionParticles = 0.1;        // extrema of the PSD (rP = random)
   double densityParticles = 1500.0;         // particle density
   double densityWall = 5000.0;        // walls density

   // INTERACTION PARAMETERS (LUDING COHESION MODEL)
   double plasticStiffness = 200.0;       // plastic stiffness (when below plasticity depth threshold)
   double elasticStiffness = 1000.0;         // elastic stiffness (when above plasticity depth threshold)
   double interclusterCohesionStiffness = 10.0;       // cohesion stiffness between particles of the same cluster (should be relatively high)
   double looseParticleCohesionStiffness = 0.0;       // cohesion stiffness between particles belonging to different clusters (or after breakage)
   double phiParticle = 0.5;        // plasticity depth
   double particleParticleRestitutionCoefficient = 0.5;        // restitution coefficient in a particle-particle collision
   double particleWallRestitutionCoefficient = 0.7;         // restitution coefficient in a particle-wall collision
   double elasticStiffnessWalls = 1000.0;       // walls elastic stiffness (usually >= than particle elastic stiffness)

   // FRICTION PARAMETERS
   double muSlidingParticleParticle = 0.5;         // sliding friction coefficient in particle-particle interaction
   double muRollingParticleParticle = 0.3;         // rolling friction coefficient in particle-particle interaction
   double muSlidingParticleWall = 0.5;       // sliding friction coefficient in particle-wall interaction
   double muRollingParticleWall = 0.3;       // rolling friction coefficient in particle-wall interaction

   // FORCE AND DAMPING PARAMETERS
   double maximumForceMultiplicationFactor = 1000.0;        // multiplication constant for the central local force responsible for cluster agglomeration
   double finalForceMultiplicationFactor = 0.0001;       // final multiplication constant of the agglomeration force before the former is turned off
   double forceTuningInterval = timeStep;       // every how much time the force is increased/decreased
   double velocityDampingModulus = 0.9;         // multiplication constant for the particle velocity when damping the particle velocity (after agglomeration)
   double velocityDampingInterval = 0.01;       // every how much time the velocity is damped

   // GRANULES PARAMETERS
   int numberOfParticlesPerCluster = 100;         // number of particles per cluster
   double clusterSizeToLatticeRatio = 0.5;         // ratio between assumed (overestimated) cluster size and its lattice size (must be < 1.0)
   int nXLatticeSites = 1;       // number of initial cluster lattice sites in x direction
   int nYLatticeSites = 1;       // number of initial cluster lattice sites in y direction
   int nZLatticeSites = 10;       // number of initial cluster lattice sites in z direction

   // INTERNAL STRUCTURE ANALYSIS PARAMETERS
   int gridLength = 300;         // length of the grid used to analyse the cluster internal structure: the higher the better the spatial resolution (300 ensures ~99% precision)
   double fictiousGridParticleRadiusRatio = 1.0e-5;         // ratio between a fictious particle used to test the structure and the mean particle size (so I can use the checkParticleForInteraction() function)

   // COMPRESSION PARAMETERS
   double maximumCompressiveForce = 10.0;        // max force to be reached during uniaxial compression
   double cylinderToClusterSizeRatio = 1.2;        // the ratio between the cylinder casing radius and the cluster radius
   double particleRadiusToPistonDisplacementPerStepRatio = 10000.0;        // the ratio between the MINIMUM particle size and the piston displacement at each time step (piston_velocity = particle_radius/ratio/time_step)

   // GLOBAL PARAMETERS
   double energyRatioTolerance = 1.0e-4;        // ratio E_kinetic/E_elastic to satisfy to consider the system as static (for settling purpouses)
   double gravityAndForceTuningDuration = 0.2;        // time interval during which force and gravity are tuned
   bool computeCdat = true;         // switches the advanced data output ON/OFF
   bool computeGrid = false;         // switches the internal granule structure analysis ON/OFF
   bool computeAmat = true;         // switches the adjacency matrix output ON/OFF

   // for (int i = 0; i < iMax; i++)   // placeholder for a loop for parameter study
   // {
         // INITIALIZATION
         dbClusters_compression problem;

         // MERCURY STANDARD FUNCTIONS
         problem.setTimeStep(timeStep);
         problem.setTimeMax(timeMax);
         problem.setGravity(Vec3D(0.00,0.00,0.00));
         problem.setSystemDimensions(3);
         problem.setSaveCount(timeStepForDataOutput/problem.getTimeStep());
         problem.setXBallsAdditionalArguments("-v0 -p 10");

         // PROBLEM SPECIFIC FUNCTIONS
         problem.setCdatOutputTimeInterval(timeStepForDataOutput);
         problem.setEnergyRatioTolerance(energyRatioTolerance);
         problem.setClustersProperties(numberOfParticlesPerCluster, clusterSizeToLatticeRatio, nXLatticeSites, nYLatticeSites, nZLatticeSites);
         problem.setParticleProperties(particleRadius, sizeDispersionParticles, densityParticles);
         problem.setWallDensity(densityWall);
         problem.setParticleParticleFrictionCoefficients(muSlidingParticleParticle, muRollingParticleParticle, 0.0);
         problem.setParticleWallFrictionCoefficients(muSlidingParticleWall, muRollingParticleWall, 0.0);
         problem.setParticleElastoPlasticProperties(plasticStiffness, elasticStiffness, interclusterCohesionStiffness, looseParticleCohesionStiffness, phiParticle);
         problem.setParticleParticleRestitutionCoefficients(particleParticleRestitutionCoefficient);
         problem.setWallStiffnessAndRestitutionCoefficients(elasticStiffnessWalls, particleWallRestitutionCoefficient);
         problem.setForceAndVelocityProperties(maximumForceMultiplicationFactor*particleRadius*elasticStiffness, finalForceMultiplicationFactor*particleRadius*elasticStiffness, forceTuningInterval, velocityDampingModulus, velocityDampingInterval);
         problem.setGravityAndForceTuningDuration(gravityAndForceTuningDuration);
         problem.setInternalStructureAnalysisParameters(gridLength, fictiousGridParticleRadiusRatio);
         problem.setCompressionParameters(maximumCompressiveForce, cylinderToClusterSizeRatio, particleRadiusToPistonDisplacementPerStepRatio);
         problem.setAdditionalDataOutput(computeCdat, computeGrid, computeAmat);

         // NAME SETTING
         std::ostringstream name;
         name.str("");
         name.clear();
         std::cout.unsetf(std::ios::floatfield);
         name << "Clusters_uniaxialCompression_" << numberOfParticlesPerCluster << "_" << nXLatticeSites << "_" << nYLatticeSites << "_" << nZLatticeSites
         << "_pR_" << particleRadius << "_rDisp_" << sizeDispersionParticles
         << "_kP_" << plasticStiffness << "_kE_" << elasticStiffness
         << "_kCi_" << interclusterCohesionStiffness << "_kCe_" << looseParticleCohesionStiffness << "_phi_" << phiParticle
         << "_muS_" << muSlidingParticleParticle << "_muR_" << muRollingParticleParticle
         << "_fRatio_" << maximumForceMultiplicationFactor << "_Fmax_" << maximumCompressiveForce;
         problem.setName(name.str());

         // RUNS THE CODE
         problem.solve();
   // }

   return 0;
}
