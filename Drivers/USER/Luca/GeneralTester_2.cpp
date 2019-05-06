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
   double nParticles = 1000;
   double particleRadius = 0.01;
   double particleDensity = 1500.0;
   double e = 0.6;
   double k1 = 500.0;
   double k2 = 5000.0;
   double kC = 50000.0;
   double kC_2 = 0.0;
   double phi = 0.1;

   double muPP_S = 0.50;
   double muPP_R = 0.05;
   double muPP_T = 0.00;

   // geometry properties
   double cylinderRadius = 0.1;
   double cylinderHeight = 2.0*cylinderRadius;
   double omega = 0.5*constants::pi;
   int nWalls = 3;
   double angle = 0.0;

   void setupInitialConditions() override
   {
      setXMin(-cylinderRadius);
      setYMin(-cylinderRadius);
      setZMin(0.0);

      setXMax(cylinderRadius);
      setYMax(cylinderRadius);
      setZMax(cylinderHeight);

      // -------------------------------------------------------------

      setSpecies();
      makeWalls();
      makeParticles();

      // -------------------------------------------------------------

      std::ostringstream cdatName;
      std::cout.unsetf(std::ios::floatfield);
      cdatName << getName() << ".cdat";

      cdatFile.open(cdatName.str(), std::ios::out);
   }

   void actionsAfterTimeStep() override
   {
      if (fmod(getTime(), 0.01) < getTimeStep()) makeDataAnalysis();
      makeRotation();
   }

   void actionsAfterSolve() override
   {
      cdatFile.close();
   }

public:

   void setSpecies()
   {
      speciesHandler.clear();
      specieParticle1 = new LinearPlasticViscoelasticFrictionSpecies;
      specieParticle1 -> setDensity(particleDensity);
      specieParticle1 -> setStiffnessAndRestitutionCoefficient(k1, e, 4.0*pow(particleRadius,3.0)*particleDensity);
      specieParticle1 -> setUnloadingStiffnessMax(k2);
      specieParticle1 -> setCohesionStiffness(kC);
      specieParticle1 -> setPenetrationDepthMax(phi);

      specieParticle1 -> setSlidingFrictionCoefficient(muPP_S);
      specieParticle1 -> setSlidingStiffness(specieParticle1 -> getLoadingStiffness()*2.0/7.0);
      specieParticle1 -> setSlidingDissipation(specieParticle1 -> getDissipation()*2.0/7.0);
      specieParticle1 -> setRollingFrictionCoefficient(muPP_R);
      specieParticle1 -> setRollingStiffness(specieParticle1 -> getLoadingStiffness()*2.0/7.0);
      specieParticle1 -> setRollingDissipation(specieParticle1 -> getDissipation()*2.0/7.0);
      specieParticle1 -> setTorsionFrictionCoefficient(muPP_T);
      specieParticle1 -> setTorsionStiffness(specieParticle1 -> getLoadingStiffness()*2.0/7.0);
      specieParticle1 -> setTorsionDissipation(specieParticle1 -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(specieParticle1);

      // -------------------------------------------------------------

      specieParticle2 = new LinearPlasticViscoelasticFrictionSpecies;
      specieParticle2 -> setDensity(particleDensity);
      specieParticle2 -> setStiffnessAndRestitutionCoefficient(k1, e, 4.0*pow(particleRadius,3.0)*particleDensity);
      specieParticle2 -> setUnloadingStiffnessMax(k2);
      specieParticle2 -> setCohesionStiffness(kC_2);
      specieParticle2 -> setPenetrationDepthMax(phi);

      specieParticle2 -> setSlidingFrictionCoefficient(muPP_S);
      specieParticle2 -> setSlidingStiffness(specieParticle2 -> getLoadingStiffness()*2.0/7.0);
      specieParticle2 -> setSlidingDissipation(specieParticle2 -> getDissipation()*2.0/7.0);
      specieParticle2 -> setRollingFrictionCoefficient(muPP_R);
      specieParticle2 -> setRollingStiffness(specieParticle2 -> getLoadingStiffness()*2.0/7.0);
      specieParticle2 -> setRollingDissipation(specieParticle2 -> getDissipation()*2.0/7.0);
      specieParticle2 -> setTorsionFrictionCoefficient(muPP_T);
      specieParticle2 -> setTorsionStiffness(specieParticle2 -> getLoadingStiffness()*2.0/7.0);
      specieParticle2 -> setTorsionDissipation(specieParticle2 -> getDissipation()*2.0/7.0);
      speciesHandler.addObject(specieParticle2);

      // -------------------------------------------------------------

      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setStiffnessAndRestitutionCoefficient(k1, e, 4.0*pow(particleRadius,3.0)*particleDensity);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setUnloadingStiffnessMax(k2);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setCohesionStiffness(0.5*(kC + kC_2));
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setPenetrationDepthMax(phi);

      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setSlidingFrictionCoefficient(muPP_S);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setSlidingStiffness(speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setSlidingDissipation(speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setRollingFrictionCoefficient(muPP_R);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setRollingStiffness(speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setRollingDissipation(speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> getDissipation()*2.0/7.0);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setTorsionFrictionCoefficient(muPP_T);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setTorsionStiffness(speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> getLoadingStiffness()*2.0/7.0);
      speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> setTorsionDissipation(speciesHandler.getMixedObject(specieParticle1, specieParticle2) -> getDissipation()*2.0/7.0);

      // -------------------------------------------------------------
   }

   void makeWalls()
   {
      wallHandler.clear();

      wall.setSpecies(specieParticle1);
      wall.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
      wallHandler.copyAndAddObject(wall);
      wall.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
      wallHandler.copyAndAddObject(wall);

      cylinder.setSpecies(specieParticle1);
      cylinder.setPosition(Vec3D(0.0,0.0,0.0));
      cylinder.setOrientation(Vec3D(0.0,0.0,1.0));
      cylinder.addObject(Vec3D(1.0,0.0,0.0),Vec3D(cylinderRadius,0.0,0.0));
      cylPointer = wallHandler.copyAndAddObject(cylinder);

      for (int i = 0; i < nWalls; i++)
      {
         CompressionPistonSurface w;

         w.setSpecies(specieParticle1);
         w.set(1.5*cylinderHeight, 0.5*cylinderHeight, 0.01*cylinderRadius, constants::pi*i/nWalls);
         wVector.push_back(w);
         wPointerVector.push_back(wallHandler.copyAndAddObject(wVector[i]));
      }
   }

   void makeParticles()
   {
      p0.setRadius(particleRadius);
      p0.setSpecies(specieParticle1);
      p0.setVelocity(Vec3D(0.0,0.0,0.0));

      // p0.setPosition(Vec3D(0.0,0.5*cylinderRadius,0.5*cylinderHeight));
      // particleHandler.copyAndAddObject(p0);
      // p0.setPosition(Vec3D(0.0,-0.5*cylinderRadius,0.5*cylinderHeight));
      // particleHandler.copyAndAddObject(p0);

      double rad, theta;

      for (int i=0; i < nParticles; i++)
      {
         rad = random.getRandomNumber(p0.getRadius(),cylinderRadius - 1.01*(p0.getRadius()));
         theta = constants::pi*random.getRandomNumber(-1.0,1.0);
         p0.setPosition(Vec3D(rad*cos(theta),rad*sin(theta),0.5*cylinderHeight*(1 + random.getRandomNumber(-1.0,1.0))));
         particleHandler.copyAndAddObject(p0);
      }
   }

   void makeRotation()
   {
      if (getTime() < 2.0)
      {
         // lowering the piston
         for (int i = 0; i < wPointerVector.size(); i++) wPointerVector[i] -> move(-0.5*cylinderHeight*getTimeStep()/2.0);

         // wiggle
         angle = 0.5*constants::pi/1.5*std::sin(2.0*constants::pi*(getTime() - 2.0)/1.5)*getTimeStep();
         for (int i = 0; i < wPointerVector.size(); i++) wPointerVector[i] -> rotate(angle);
      }
      else
      {
         // rotation
         angle += omega*getTimeStep();
         for (int i = 0; i < wPointerVector.size(); i++) wPointerVector[i] -> rotate(omega*getTimeStep());
      }
   }

   // analyzes simulations data of interest
   void makeDataAnalysis()
   {
      // meanCoordinationNumber = 0.0;
      //
      // for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
      // {
      //    meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
      // }
      //
      // meanCoordinationNumber /= particleHandler.getNumberOfObjects();

      // pistonForce = 0.0;
      // pistonPressure = 0.0;
      // basePressure = 0.0;
      // cylinderPressure = 0.0;
      // meanTotalRelativeOverlap = 0.0;
      // maxTotalRelativeOverlap = 0.0;
      // meanPistonRelativeOverlap = 0.0;
      // maxPistonRelativeOverlap = 0.0;
      // meanBaseRelativeOverlap = 0.0;
      // maxBaseRelativeOverlap = 0.0;
      // meanCylinderRelativeOverlap = 0.0;
      // maxCylinderRelativeOverlap = 0.0;
      //
      // int totalInteractionCounter = 0;
      // int pistonInteractionCounter = 0;
      // int baseInteractionCounter = 0;
      // int cylinderInteractionCounter = 0;

      surfaceForce = 0.0;
      surfaceTorque1 = 0.0;
      surfaceTorqueAbs1 = 0.0;
      surfaceTorque2 = 0.0;
      surfaceTorqueAbs2 = 0.0;
      surfaceTorqueTot = 0.0;
      surfaceTorqueAbsTot = 0.0;

      for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
      {
         for (int j = 0; j < wPointerVector.size(); j++)
         {
            // piston interactions
            if ((*i) -> getI() -> getIndex() == wPointerVector[j] -> getIndex())
            {
               surfaceTorqueTot += ((*i) -> getContactPoint()).X*((*i) -> getForce()).Y - ((*i) -> getContactPoint()).Y*((*i) -> getForce()).X;
            }
         }
      }

      cdatFile << getTime() << "\t" << surfaceTorqueTot << "\t" << angle << std::endl;

      // pistonForce = pistonPressure;
      // pistonPressure /= constants::pi*pow(casingRadius,2.0);
      // basePressure /= constants::pi*pow(casingRadius,2.0);
      // if (pistonHeight > particleBedHeight) cylinderPressure /= 2.0*constants::pi*casingRadius*particleBedHeight;
      // else cylinderPressure /= 2.0*constants::pi*casingRadius*pistonHeight;
      // meanTotalRelativeOverlap /= totalInteractionCounter;
      // meanPistonRelativeOverlap /= pistonInteractionCounter;
      // meanBaseRelativeOverlap /= baseInteractionCounter;
      // meanCylinderRelativeOverlap /= cylinderInteractionCounter;
   }

   //  ----- GLOBAL FUNCTIONS -----
   void printTime() const override
   {
      std::cout << "t = " << std::setprecision(3) << std::left << std::setw(6) << getTime() << ", tmax = " << std::setprecision(3) << std::left << std::setw(4) << getTimeMax() <<
      ", h_1 = " << wPointerVector[0] -> getPosition() <<
      ", alpha_1 = " << wPointerVector[0] -> getAngle() <<
      ", tau = " << surfaceTorqueTot << std::endl;
      std::cout.flush();
   }


   // VARIABLES -------------------------------------------------------
   double t0;
   int numberOfParticlesExited;

   LinearPlasticViscoelasticFrictionSpecies *specieParticle1;
   LinearPlasticViscoelasticFrictionSpecies *specieParticle2;
   BaseParticle p0;
   InfiniteWall wall;
   AxisymmetricIntersectionOfWalls cylinder;
   AxisymmetricIntersectionOfWalls* cylPointer;

   std::vector<CompressionPistonSurface> wVector;
   std::vector<CompressionPistonSurface*> wPointerVector;

   double surfaceForce;
   double surfaceTorque1;
   double surfaceTorqueAbs1;
   double surfaceTorque2;
   double surfaceTorqueAbs2;
   double surfaceTorqueTot;
   double surfaceTorqueAbsTot;

   std::ofstream cdatFile;
};


int main(int argc, char *argv[])
{
   Tester problem;

   // INITIALIZATION
   problem.setTimeStep(1.0e-5);
   problem.setTimeMax(10.0);
   problem.setGravity(Vec3D(0.0,0.0,-9.81));
   problem.setSystemDimensions(3);
   problem.setSaveCount(0.01/problem.getTimeStep());

   problem.setName("TESTER_500_5000_50000_3");
   problem.solve();

   return 0;
}
