#include <Mercury3D.h>
#include <Particles/BaseParticle.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <math.h>
#include <fstream>

 class Tester : public Mercury3D
 {
private:
   void setupInitialConditions() override
   {
      setXMin(-5.0*particleRadius);
      setYMin(-5.0*particleRadius);
      setZMin(-5.0*particleRadius);

      setXMax(5.0*particleRadius);
      setYMax(5.0*particleRadius);
      setZMax(5.0*particleRadius);

      // -------------------------------------------------------------

      setSpecies();
      makeWalls();
      makeParticles();

      // -------------------------------------------------------------

      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
      cdatFile << "time \t position_1 \t position_2 \t velocity_1 \t velocity_2 \t force_1-2" << std::endl;
   }

   void actionsAfterTimeStep() override
   {
      if (fmod(getTime(), 0.001) < getTimeStep())
      {
         cdatFile <<
         getTime() << "   " <<
         particleHandler.getObject(0) -> getPosition() << "   " <<
         particleHandler.getObject(0) -> getVelocity() << "   " <<
         particleHandler.getObject(0) -> getForce() << std::endl;
      }
   }

   void actionsAfterSolve() override
   {
      cdatFile.close();
   }

public:
   void setParticleProperties(double pR, double density, double e, double kp, double ke, double kc, double ph)
   {
      particleRadius = pR;
      particleDensity = density;
      particleRestitutionCoefficient = e;
      kP = kp;
      kE = ke;
      kC = kc;
      phi = ph;
   }

   void setSpecies()
   {
      speciesHandler.clear();
      speciesParticle = new LinearPlasticViscoelasticFrictionSpecies;
      speciesParticle -> setDensity(particleDensity);
      speciesParticle -> setStiffnessAndRestitutionCoefficient(kP, particleRestitutionCoefficient, 4.0*pow(particleRadius,3.0)*particleDensity);
      speciesParticle -> setUnloadingStiffnessMax(kE);
      speciesParticle -> setCohesionStiffness(kC);
      speciesParticle -> setPenetrationDepthMax(phi);

      speciesParticle -> setSlidingFrictionCoefficient(0.0);
      speciesParticle -> setSlidingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setSlidingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setRollingFrictionCoefficient(0.0);
      speciesParticle -> setRollingStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setRollingDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesParticle -> setTorsionFrictionCoefficient(0.0);
      speciesParticle -> setTorsionStiffness(speciesParticle -> getLoadingStiffness()*2.0/7.0);
      speciesParticle -> setTorsionDissipation(speciesParticle -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(speciesParticle);
      // -------------------------------------------------------------
   }

   void makeWalls()
   {
      wallHandler.clear();

      wall.setSpecies(speciesParticle);
      wall.set(Vec3D(1.0,0.0,0.0),Vec3D(0.0,0.0,getXMax()));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(-1.0,0.0,0.0),Vec3D(0.0,0.0,getXMin()));
      wallHandler.copyAndAddObject(wall);

      wall.setSpecies(speciesParticle);
      wall.set(Vec3D(0.0,1.0,0.0),Vec3D(0.0,0.0,getYMax()));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,0.0,getYMin()));
      wallHandler.copyAndAddObject(wall);

      wall.setSpecies(speciesParticle);
      wall.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
      wallHandler.copyAndAddObject(wall);
   }

   void makeParticles()
   {
      p0.setRadius(particleRadius);
      p0.setSpecies(speciesParticle);

      // particle X-
      p0.setPosition(Vec3D(0.0, 0.0, 0.0));
      p0.setVelocity(Vec3D(0.0000,0.0,0.0));
      particleHandler.copyAndAddObject(p0);

      // particle X+
      // p0.setPosition(Vec3D(0.035, 0.0, 0.0));
      // p0.setVelocity(Vec3D(-0.0000,0.0,0.0));
      // particleHandler.copyAndAddObject(p0);
   }

   //  ----- GLOBAL FUNCTIONS -----
   // void printTime() const override
   // {
   //    std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() << std::endl;
   // }


   // VARIABLES -------------------------------------------------------
   LinearPlasticViscoelasticFrictionSpecies *speciesParticle;
   SphericalParticle p0;
   InfiniteWall wall;
   std::ofstream cdatFile;

   double particleRadius;
   double particleDensity;
   double particleRestitutionCoefficient;
   double kP;
   double kE;
   double kC;
   double phi;
};


int main(int argc, char *argv[])
{
   double particleRadius = 0.01;
   double particleDensity = 1500.0;
   double particleRestitutionCoefficient = 1.;
   double kP = 500.0;
   double kE = 1000.0;
   double kC = 500.0;
   double phi = 0.10;

   Tester collisionTester;

   // -------------------------------------------------------------

   // INITIALIZATION
   collisionTester.setTimeStep(1.0e-6);
   collisionTester.setTimeMax(2.0);
   collisionTester.setGravity(Vec3D(0.0,0.0,-10.0));
   collisionTester.setSystemDimensions(3);
   collisionTester.setSaveCount(0.001/collisionTester.getTimeStep());

   collisionTester.setParticleProperties(particleRadius, particleDensity, particleRestitutionCoefficient, kP, kE, kC, phi);

   // NAME SETTING
   std::ostringstream name;
   name.str("");
   name.clear();
   std::cout.unsetf(std::ios::floatfield);
   name << "BinaryCollisionTester___kP_" << kP << "_kE_" << kE << "_kC_" << kC << "_phi_" << phi << "_rho_" << particleDensity << "_e_" << particleRestitutionCoefficient;
   collisionTester.setName(name.str());

   collisionTester.solve();

   return 0;
}
