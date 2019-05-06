#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <math.h>
#include <fstream>
#include <Walls/IntersectionOfWalls.h>
#include "CompressionPistonSurface.h"
#include <vector>
#include <numeric>

/*
 * ToDo:
 * - make the box size dependent on the particle max size
 */


 class Tester : public Mercury3D
 {
private:
   // particle properties
   double particleRadiusBig = 0.001;
   double particleRadiusSmall = particleRadiusBig/10.0;
   double particleDensity = 1500.0;
   double e = 0.4;
   double k1 = 500.0;
   double k2 = 1000.0;
   double kCbig = 100.0;
   double kCsmall = 500.0;
   double phi = 0.10;

   double muPP_S = 0.50;
   double muPP_R = 0.05;
   double muPP_T = 0.00;

   double maximumForce = k2*2.0*particleRadiusBig*particleRadiusSmall/(particleRadiusBig + particleRadiusSmall)*1.0;
   double minimumForce = k2*2.0*particleRadiusBig*particleRadiusSmall/(particleRadiusBig + particleRadiusSmall)*0.0001;
   double forceUnloadingMultiplicationFactor = 1.0001;
   double forceLoadingDuration = 5.0;

   int nSmallParticles = 50;
   int stage = 0;

   void setupInitialConditions() override
   {
      setXMin(-2.0*particleRadiusBig);
      setYMin(-2.0*particleRadiusBig);
      setZMin(-2.0*particleRadiusBig);

      setXMax(2.0*particleRadiusBig);
      setYMax(2.0*particleRadiusBig);
      setZMax(2.0*particleRadiusBig);

      // -------------------------------------------------------------

      setSpecies();
      // makeWalls();
      insertParticles();

      // -------------------------------------------------------------

      // std::ostringstream cdatName;
      // std::cout.unsetf(std::ios::floatfield);
      // cdatName << getName() << ".cdat";
      //
      // cdatFile.open(cdatName.str(), std::ios::out);
   }

   void actionsAfterTimeStep() override
   {
      // if (fmod(getTime(), 0.01) < getTimeStep()) makeDataAnalysis();
      particleHandler.getObject(0) -> setPosition(Vec3D(0.,0.,0.));

      if (stage == 0)
      {

         // applyCentralForce();

         if (fmod(getTime() - t0,10.0*getTimeStep()) < getTimeStep()) increaseForce();

         if (forceModulus > maximumForce)
         {
            std::cout << "Maximum force reached." << std::endl;
            std::cout << "DAMPING CENTRAL FORCE" << std::endl << std::endl;

            t0 = getTime();
            stage++;
         }
      }

      // DAMPING OF LOCAL FORCE
      if (stage == 1)
      {
         applyCentralForce();

         if (fmod(getTime() - t0,10.0*getTimeStep()) < getTimeStep()) dampForce();

         if (forceModulus < minimumForce)
         {
            std::cout << "Minimum force reached." << std::endl;

            std::cout << "DISSIPATING ENERGY" << std::endl << std::endl;
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

   void setSpecies()
   {
      speciesHandler.clear();
      speciesParticleBig = new LinearPlasticViscoelasticFrictionSpecies;
      speciesParticleBig -> setDensity(particleDensity);
      speciesParticleBig -> setStiffnessAndRestitutionCoefficient(k1, e, 4.0*pow(particleRadiusBig,3.0)*particleDensity);
      speciesParticleBig -> setUnloadingStiffnessMax(k2);
      speciesParticleBig -> setCohesionStiffness(kCbig);
      speciesParticleBig -> setPenetrationDepthMax(phi);

      speciesParticleBig -> setSlidingFrictionCoefficient(muPP_S);
      speciesParticleBig -> setSlidingStiffness(speciesParticleBig -> getLoadingStiffness()*2.0/7.0);
      speciesParticleBig -> setSlidingDissipation(speciesParticleBig -> getDissipation()*2.0/7.0);
      speciesParticleBig -> setRollingFrictionCoefficient(muPP_R);
      speciesParticleBig -> setRollingStiffness(speciesParticleBig -> getLoadingStiffness()*2.0/7.0);
      speciesParticleBig -> setRollingDissipation(speciesParticleBig -> getDissipation()*2.0/7.0);
      speciesParticleBig -> setTorsionFrictionCoefficient(muPP_T);
      speciesParticleBig -> setTorsionStiffness(speciesParticleBig -> getLoadingStiffness()*2.0/7.0);
      speciesParticleBig -> setTorsionDissipation(speciesParticleBig -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesParticleBig);

      // -------------------------------------------------------------

      speciesParticleSmall = new LinearPlasticViscoelasticFrictionSpecies;
      speciesParticleSmall -> setDensity(particleDensity);
      speciesParticleSmall -> setStiffnessAndRestitutionCoefficient(k1, e, 4.0*pow(particleRadiusSmall,3.0)*particleDensity);
      speciesParticleSmall -> setUnloadingStiffnessMax(k2);
      speciesParticleSmall -> setCohesionStiffness(kCsmall);
      speciesParticleSmall -> setPenetrationDepthMax(phi);

      speciesParticleSmall -> setSlidingFrictionCoefficient(muPP_S);
      speciesParticleSmall -> setSlidingStiffness(speciesParticleSmall -> getLoadingStiffness()*2.0/7.0);
      speciesParticleSmall -> setSlidingDissipation(speciesParticleSmall -> getDissipation()*2.0/7.0);
      speciesParticleSmall -> setRollingFrictionCoefficient(muPP_R);
      speciesParticleSmall -> setRollingStiffness(speciesParticleSmall -> getLoadingStiffness()*2.0/7.0);
      speciesParticleSmall -> setRollingDissipation(speciesParticleSmall -> getDissipation()*2.0/7.0);
      speciesParticleSmall -> setTorsionFrictionCoefficient(muPP_T);
      speciesParticleSmall -> setTorsionStiffness(speciesParticleSmall -> getLoadingStiffness()*2.0/7.0);
      speciesParticleSmall -> setTorsionDissipation(speciesParticleSmall -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesParticleSmall);

      // -------------------------------------------------------------

      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setStiffnessAndRestitutionCoefficient(k1, e, 4.0*pow(0.5*(particleRadiusBig + particleRadiusSmall),3.0)*particleDensity);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setUnloadingStiffnessMax(k2);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setCohesionStiffness(0.5*(kCbig + kCsmall));
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setPenetrationDepthMax(phi);

      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setSlidingFrictionCoefficient(muPP_S);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setSlidingStiffness(speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setSlidingDissipation(speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setRollingFrictionCoefficient(muPP_R);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setRollingStiffness(speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setRollingDissipation(speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setTorsionFrictionCoefficient(muPP_T);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setTorsionStiffness(speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> setTorsionDissipation(speciesHandler.getMixedObject(speciesParticleBig, speciesParticleSmall) -> getDissipation()*2.0/7.0);

      // -------------------------------------------------------------
   }

   bool particleInsertionSuccessful()
   {
      int insertionFailCounter = 0;
      double rad, theta, phi;
      Vec3D particlePosition;

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setRadius(particleRadiusSmall);
      p0.setSpecies(speciesParticleSmall);

      do
      {
         rad = particleRadiusBig + 1.1*particleRadiusSmall;
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
      particleHandler.clear();

      p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
      p0.setPosition(Vec3D(0.0, 0.0, 0.0));
      p0.setRadius(particleRadiusBig);
      p0.setSpecies(speciesParticleBig);
      particleHandler.copyAndAddObject(p0);

      int nParticlesInserted = 0;

      while(nParticlesInserted < nSmallParticles)
      {
         if (particleInsertionSuccessful()) nParticlesInserted++;
         else
         {
            std::cout << "PARTICLE INSERTION FAILURE. QUITTING...";
            exit(-1);
         }
      }

      std::cout << std::endl << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
      t0 = getTime();
   }

   void makeWalls()
   {
      wallHandler.clear();
      wall.setSpecies(speciesParticleBig);

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

   void applyCentralForce()
   {
      Vec3D distanceFromForceCenter;

      for (int i=1; i<particleHandler.getNumberOfObjects(); i++)
      {
         distanceFromForceCenter = particleHandler.getObject(i) -> getPosition();

         particleHandler.getObject(i) -> addForce(-forceModulus*distanceFromForceCenter);
         particleHandler.getObject(i) -> setVelocity(0.99*(particleHandler.getObject(i) -> getVelocity()));
      }
   }

   void increaseForce()
   {
      forceModulus = maximumForce*(getTime() - t0)/forceLoadingDuration;
   }

   void dampForce()
   {
      forceModulus /= forceUnloadingMultiplicationFactor;
   }

   // // analyzes simulations data of interest
   // void makeDataAnalysis()
   // {
   //    pistonPressure = 0.0;
   //    bedHeight = (particleHandler.getHighestPositionComponentParticle(2) -> getPosition()).Z;
   //
   //    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
   //    {
   //       if ((*i) -> getI() -> getIndex() == pistonPointer -> getIndex()) pistonPressure -= ((*i) -> getForce()).Z;;
   //    }
   //
   //    pistonPressure /= constants::pi*pow(getXMax(), 2.0);
   //
   //    cdatFile << getTime() << " " << bedHeight << " " << pistonHeight << " " << pistonPressure << std::endl;
   // }

   //  ----- GLOBAL FUNCTIONS -----
   // void printTime() const override
   // {
   //    std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << getTimeMax() <<
   //    ", p_z max = " << bedHeight << ", h = " << pistonHeight << ", P = " << pistonPressure << std::endl;
   //    std::cout.flush();
   // }


   // VARIABLES -------------------------------------------------------
   double t0;
   LinearPlasticViscoelasticFrictionSpecies *speciesParticleBig;
   LinearPlasticViscoelasticFrictionSpecies *speciesParticleSmall;
   BaseParticle p0;
   InfiniteWall wall;
   Vec3D f1;
   Vec3D f2;

   double pistonPressure;
   double bedHeight;
   double forceModulus;

   AxisymmetricIntersectionOfWalls* cylPointer;
   InfiniteWall* pistonPointer;

   std::ofstream cdatFile;
};


int main(int argc, char *argv[])
{
   Tester problem;

   // INITIALIZATION
   problem.setTimeStep(1.0e-5);
   problem.setTimeMax(7.0);
   problem.setGravity(Vec3D(0.0,0.0,0.0));
   problem.setSystemDimensions(3);
   problem.setSaveCount(0.01/problem.getTimeStep());

   problem.setName("TESTER_3");
   problem.solve();

   return 0;
}
