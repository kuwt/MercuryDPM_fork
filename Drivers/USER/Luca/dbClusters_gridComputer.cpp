/*
 *** dbClusters GRID COMPUTER ***
 An already existing cluster restart file is loaded and the .grid file fully computed and exported.
 */

#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Walls/InfiniteWall.h>
#include <fstream>

/*
 ToDo:
 */

class dbClusters_gridComputer : public Mercury3D
{
private:
   void setupInitialConditions() override
   {
      std::cout << "THIS MESSAGE SHOULD NOT BE DISPLAYED." << std::endl;
      std::cout << "SINCE YOU ARE READING IT SOMETHING WENT WRONG." << std::endl;
      std::cout << "EXITING..." << std::endl;

      exit(1);
   }

   void actionsOnRestart() override
   {
      std::cout << "RESTARTING" << std::endl;

      speciesParticle = new LinearPlasticViscoelasticFrictionSpecies;
      speciesParticle -> setDensity(2000.0);
      speciesParticle -> setStiffnessAndRestitutionCoefficient(1000.0, 1.0, 0.01);
      speciesParticle -> setUnloadingStiffnessMax(2000.0);
      speciesParticle -> setCohesionStiffness(0.0);
      speciesParticle -> setPenetrationDepthMax(0.05);

      speciesParticle -> setSlidingFrictionCoefficient(0.5);
      speciesParticle -> setSlidingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setSlidingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setRollingFrictionCoefficient(0.3);
      speciesParticle -> setRollingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setRollingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setTorsionFrictionCoefficient(0.0);
      speciesParticle -> setTorsionStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setTorsionDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesParticle);

      std::cout << particleHandler.getNumberOfObjects() << std::endl;

      makeDataAnalysis();
      makeGridFileAndComputeInternalStructure();
   }

   void actionsAfterTimeStep() override
   {

   }

   void actionsAfterSolve() override
   {

   }


public:
   // FUNCTIONS CALLED IN MAIN ----------------------------------------
   void setClustersProperties(double axesRatio)
   {
      granuleAxesRatio = axesRatio;
   }

   void setInternalStructureAnalysisParameters(int length, double ratio)
   {
      gridLength = length;
      fictiousGridPointRadiusRatio = ratio;
   }


   // FUNCTIONS CALLED IN THE CLASS -----------------------------------
   void makeDataAnalysis()
   {
      double totalParticleVolume = 0.0;
      Vec3D distanceFromForceCenter;
      Vec3D localMin, localMax;
      localMin.setZero();
      localMax.setZero();
      centerOfMass.setZero();

      radiusParticle = 0.0;
      meanClusterRadius = 0.0;

      for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      {
         totalParticleVolume += pow(particleHandler.getObject(i) -> getRadius(), 3.0);
         centerOfMass += pow(particleHandler.getObject(i) -> getRadius(), 3.0)*(particleHandler.getObject(i) -> getPosition());
         radiusParticle += particleHandler.getObject(i) -> getRadius();
         distanceFromForceCenter = particleHandler.getObject(i) -> getPosition();
         if (distanceFromForceCenter.X > localMax.X) localMax.X = distanceFromForceCenter.X;
         if (distanceFromForceCenter.X < localMin.X) localMin.X = distanceFromForceCenter.X;
         if (distanceFromForceCenter.Y > localMax.Y) localMax.Y = distanceFromForceCenter.Y;
         if (distanceFromForceCenter.Y < localMin.Y) localMin.Y = distanceFromForceCenter.Y;
         if (distanceFromForceCenter.Z > localMax.Z) localMax.Z = distanceFromForceCenter.Z;
         if (distanceFromForceCenter.Z < localMin.Z) localMin.Z = distanceFromForceCenter.Z;
      }

      centerOfMass /= totalParticleVolume;
      radiusParticle /= particleHandler.getNumberOfObjects();
      meanClusterRadius = (localMax.X - localMin.X + localMax.Y - localMin.Y + localMax.Z - localMin.Z)/6.0;
   }

   // creates the grid file and writes the first row
   void makeGridFileAndComputeInternalStructure()
   {
      Vec3D gridPoint;
      SphericalParticle p0;
      double nPointsInsideAgglomerateBoundary;
      double nPointsInsideComponents;

      gridPoint.setZero();
      gridResolutionLength = 2.0*meanClusterRadius/(gridLength - 1.0);
      nPointsInsideAgglomerateBoundary = 0;
      nPointsInsideComponents = 0;

      p0.setRadius(radiusParticle*fictiousGridPointRadiusRatio);
      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setSpecies(speciesParticle);

      std::ostringstream gridName;
      std::cout.unsetf(std::ios::floatfield);
      gridName << getName() << ".grid";
      gridFile.open(gridName.str(), std::ios::out);
      gridFile << "Grid_length: " << gridLength << "\t grid_resolution_length: " << gridResolutionLength << "\t fictiousGridParticleRadiusRatio: " << fictiousGridPointRadiusRatio << std::endl;

      // creates fictious particles at every grid node
      // if they are inside the agglomerate perimeter then save them in the .grid file
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

         std::cout << i << "/" << gridLength << std::endl;
      }

      massFraction = nPointsInsideComponents/nPointsInsideAgglomerateBoundary;

      gridFile << "n_points_inside_boundary: " << nPointsInsideAgglomerateBoundary << "\t n_points_inside_components: " << nPointsInsideComponents << "\t mass_fraction: " << massFraction << std::endl;
      gridFile.close();
   }


   // VARIABLES -------------------------------------------------------
   // particles
   double radiusParticle;

   // granules
   int totalNumberOfParticles;
   double meanClusterRadius;
   Vec3D centerOfMass;

   // granule shape
   double granuleAxesRatio;

   // species
   LinearPlasticViscoelasticFrictionSpecies *speciesParticle;

   // output
   std::ofstream gridFile;

   // data analysis
   double massFraction;
   double gridLength;
   double gridResolutionLength;
   double fictiousGridPointRadiusRatio;
};


int main(int argc, char *argv[])
{
   // SHAPE PARAMETERS
   double axesRatio = 1.0;

   // INTERNAL STRUCTURE ANALYSIS PARAMETERS
   int gridLength = 300;
   double fictiousGridParticleRadiusRatio = 1.0e-5;

   // INITIALIZATION
   dbClusters_gridComputer problem;
   problem.readArguments(argc, argv);

   // NAME SETTING
   std::ostringstream name;
   name.str("");
   name.clear();
   std::cout.unsetf(std::ios::floatfield);
   // name << argv[2];
   name << "gridComputer_1000_010";
   problem.setName(name.str());

   problem.setTimeStep(0.01);
   problem.setTimeMax(problem.getTimeStep());
   problem.setGravity(Vec3D(0.00,0.00,0.00));
   problem.setSystemDimensions(3);
   problem.setSaveCount(0.01/problem.getTimeStep());

   problem.setClustersProperties(axesRatio);
   problem.setInternalStructureAnalysisParameters(gridLength, fictiousGridParticleRadiusRatio);

   problem.solve();

   return 0;
}
