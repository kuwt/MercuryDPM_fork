/* Perry's Comb for fingering experiments. Based on Perry Harwood's experiment
 * and Binbin Jin's fingering code (Fingering.cpp). */

#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/Combtooth.h"
#include "Boundaries/DeletionBoundary.h"
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

class PerryComb : public Mercury3D { 
  public:
    PerryComb(string parsfile) {
      /* Reading config file */
      ifstream file(parsfile);
      string name;
      double var;
      while (file >> name >> var) {
        pars[name] = var;
      }
      if (pars.size()!=27) {
        logger(ERROR, ".config file not read correctly.");
        exit(-1);
      }

      /* Initialization */
      generator.setRandomSeed(int(pars["randomSeed"]));
      if (pars["fStatFile"]==0) {
        fStatFile.setFileType(FileType::NO_FILE);
      }
      if (pars["dataFile"]==0) {
        dataFile.setFileType(FileType::NO_FILE);
        logger(WARN, ".data file will not be generated.");
      } else if (pars["dataFile"]!=1) {
        dataFile.setFileType(FileType::MULTIPLE_FILES);

      }
      setTimeStep(pars["timeStep"]);
      setTimeMax(pars["timeMax"]);
      setSaveCount(int(pars["saveCount"]));
      setName(parsfile.erase(parsfile.find_last_of('.')));
      setSystemDimensions(3);
      setGravity(Vec3D(0,-pars["g"],0));
      setXMin(0);
      setXMax(pars["length"]*20);
      setYMin(-pars["length"]*20*tan(pars["alpha"]));
      setYMax(pars["height"]);
      setZMin(-pars["width"]*1.2);
      setZMax(pars["width"]*1.2);
      setXBallsAdditionalArguments("-noborder 4 -v0 -solidf");

      /* Set up species for smaller particles */
      small_ = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
      small_->setDensity(pars["rho"]);
      small_->setCollisionTimeAndRestitutionCoefficient(
        pars["collisionTime"],
        pars["restitutionCoeff"], 
        (4.0/3.0)*constants::pi*pow(pars["radius_small"],3)*pars["rho"]
      );
      small_->setSlidingFrictionCoefficient(pars["sliding_small"]);
      small_->setSlidingStiffness(2.0/7.0 * small_->getStiffness());
      small_->setSlidingDissipation(2.0/7.0 * small_->getDissipation());
      small_->setRollingFrictionCoefficient(pars["rolling_small"]);
      small_->setRollingStiffness(2.0/5.0 * small_->getStiffness());
      small_->setRollingDissipation(2.0/5.0 * small_->getDissipation());
      //small_->setTorsionFrictionCoefficient(pars["friction_small"]);
      //small_->setTorsionStiffness(2.0/5.0 * small_->getStiffness());
      //small_->setTorsionDissipation(2.0/5.0 * small_->getDissipation());

      /* Set up species for larger particles */
      large_ = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
      large_->setDensity(pars["rho"]);
      large_->setCollisionTimeAndRestitutionCoefficient(
        pars["collisionTime"],
        pars["restitutionCoeff"], 
        (4.0/3.0)*constants::pi*pow(pars["radius_large"],3)*pars["rho"]
      );
      large_->setSlidingFrictionCoefficient(pars["sliding_large"]);
      large_->setSlidingStiffness(2.0/7.0 * large_->getStiffness());
      large_->setSlidingDissipation(2.0/7.0 * large_->getDissipation());
      large_->setRollingFrictionCoefficient(pars["rolling_large"]);
      large_->setRollingStiffness(2.0/5.0 * large_->getStiffness());
      large_->setRollingDissipation(2.0/5.0 * large_->getDissipation());
      //large_->setTorsionFrictionCoefficient(pars["friction_large"]);
      //large_->setTorsionStiffness(2.0/5.0 * large_->getStiffness());
      //large_->setTorsionDissipation(2.0/5.0 * large_->getDissipation());

      /* Configure the mixed species */
      auto mixed_ = speciesHandler.getMixedObject(small_,large_);
      mixed_->setCollisionTimeAndRestitutionCoefficient(
        pars["collisionTime"],
        pars["restitutionCoeff"],
        (4.0/3.0)*constants::pi*pow(pars["radius_small"],3)*pars["rho"],
        (4.0/3.0)*constants::pi*pow(pars["radius_large"],3)*pars["rho"]
      );
      mixed_->setSlidingFrictionCoefficient(pars["sliding_large"]);
      mixed_->setSlidingStiffness(2.0/7.0 * mixed_->getStiffness());
      mixed_->setSlidingDissipation(2.0/7.0 * mixed_->getDissipation());
      mixed_->setRollingFrictionCoefficient(pars["rolling_large"]);
      mixed_->setRollingStiffness(2.0/5.0 * mixed_->getStiffness());
      mixed_->setRollingDissipation(2.0/5.0 * mixed_->getDissipation());
      //mixed_->setTorsionFrictionCoefficient(pars["friction_large"]);
      //mixed_->setTorsionStiffness(2.0/5.0 * mixed_->getStiffness());
      //mixed_->setTorsionDissipation(2.0/5.0 * mixed_->getDissipation());

      /* Set up walls */
      auto leftWall = wallHandler.copyAndAddObject(IntersectionOfWalls());
      leftWall->addObject(Vec3D(0,0,-1),Vec3D(0,0,-pars["width"]));
      leftWall->addObject(Vec3D(-1,0,0),Vec3D(pars["length"],0,0));
      leftWall->setSpecies(large_);
      auto rightWall = wallHandler.copyAndAddObject(IntersectionOfWalls());
      rightWall->addObject(Vec3D(0,0,1),Vec3D(0,0,pars["width"]));
      rightWall->addObject(Vec3D(-1,0,0),Vec3D(pars["length"],0,0));
      rightWall->setSpecies(large_);
      auto backWall = wallHandler.copyAndAddObject(InfiniteWall());
      backWall->set(Vec3D(-1,0,0),Vec3D(0,0,0));
      backWall->setSpecies(large_);
      auto bottomWall = wallHandler.copyAndAddObject(InfiniteWall());
      bottomWall->set(Vec3D(-sin(pars["alpha"]),-cos(pars["alpha"]),0),Vec3D(0,0,0));
      bottomWall->setSpecies(large_);

      /* Set up initial confinement */
      frontWall = wallHandler.copyAndAddObject(IntersectionOfWalls());
      frontWall->addObject(Vec3D(1,0,0),Vec3D(pars["length"],0,0));
      frontWall->setSpecies(large_);
      initialWall = true;

      /* Set up Perry's comb */
      Combtooth tooth;
      tooth.setSpecies(small_);
      for (double z = 0; z <= pars["width"]; z += pars["comb_spacing"])
      {
          tooth.set(Vec3D(0,1,0), Vec3D(pars["comb_pos"],0,z), pars["comb_radius"]);
          wallHandler.copyAndAddObject(tooth);
          if (z != 0)
          {
              tooth.set(Vec3D(0,1,0), Vec3D(pars["comb_pos"],0,-z), pars["comb_radius"]);
              wallHandler.copyAndAddObject(tooth);
          }
      }

      /* Set up deletion boundaries */
      /* 2016-08-22: JMFT: I think these deletion boundaries are unnecessary in
       * our problem, since particles shouldn't penetrate walls or run out too
       * far if the simulation is properly configured.
       *
       * And recently I have shown in other simulations that they may be
       * responsible for segmentation faults (by leaving the particle handler in
       * an improper state after removing a particle). Hence, removing them is a
       * good safety precaution (for now).
       */
      /*
      auto db_back = boundaryHandler.copyAndAddObject(DeletionBoundary());
      db_back->set(Vec3D(-1,0,0),pars["radius_small"]);
      auto db_left = boundaryHandler.copyAndAddObject(DeletionBoundary());
      db_left->set(Vec3D(0,0,-1),pars["width"]*1.2);
      auto db_right = boundaryHandler.copyAndAddObject(DeletionBoundary());
      db_right->set(Vec3D(0,0,1),pars["width"]*1.2);
      auto db_bottom = boundaryHandler.copyAndAddObject(DeletionBoundary());
      db_bottom->set(Vec3D(-sin(pars["alpha"]),-cos(pars["alpha"]),0),
                                pars["radius_small"]);
      auto db_top = boundaryHandler.copyAndAddObject(DeletionBoundary());
      db_top->set(Vec3D(sin(pars["alpha"]),cos(pars["alpha"]),0),pars["height"]);
      auto db_front = boundaryHandler.copyAndAddObject(DeletionBoundary());
      db_front->set(Vec3D(1,0,0),pars["length"]*20);
      */
    }

    void setupInitialConditions() {
      BaseParticle p0;

      /* Set up moving particles */
      int NSmallParticle = 0, NLargeParticle = 0;
      for(double x = pars["radius_large"];
            x <= pars["length"]-pars["radius_large"];
            x += 2*pars["radius_large"]) {
        for(double y = pars["radius_large"]*3;
              y <= pars["height"]-pars["radius_large"];
              y += 2*pars["radius_large"]) {
          for(double z = -pars["width"]+pars["radius_large"];
                  z <= pars["width"]-pars["radius_large"];
                  z += 2*pars["radius_large"]) {
            if(generator.getRandomNumber(0,1)<pars["ratio_small"]) {
              addSmallParticle(p0,Vec3D(x,y-x*tan(pars["alpha"]),z));
              NSmallParticle++;
            } else {
              addLargeParticle(p0,Vec3D(x,y-x*tan(pars["alpha"]),z));
              NLargeParticle++;
            }
          }        
        }
      }
      cout << NSmallParticle << " small particles are generated." << endl;
      cout << NLargeParticle << " large particles are generated." << endl;

      kineticEnergyThreshold_ = 0.5*(4.0/3.0)*constants::pi*pars["rho"]*0.05*0.05*
        (pow(pars["radius_small"],3)*NSmallParticle+
          pow(pars["radius_large"],3)*NLargeParticle);

      /* Set up particles for rough base */
      int NBaseParticle = 0;
      p0.fixParticle();
      for(double x = pars["radius_large"];
            x <= 20*pars["length"]-pars["radius_large"];
            x += 2*pars["radius_large"]) {
        for(double z = -pars["width"]*1.2+pars["radius_large"];
                z < pars["width"]*1.2+pars["radius_large"];
                z += 2*pars["radius_large"]) {
          addLargeParticle(p0,Vec3D(x,pars["radius_large"]-x*tan(pars["alpha"]),z));
          NBaseParticle++;
        }
      }
      cout << NBaseParticle << " base particles are generated." << endl;
    }

    void actionsAfterTimeStep() {
      /* Remove initial confinement after the partciles come to rest */
      if (initialWall && getNumberOfTimeSteps() % 100 == 0 && getTime()>=0.5 &&
            getKineticEnergy() < kineticEnergyThreshold_) {
        frontWall->addObject(Vec3D(-1,1,0),Vec3D(
          pars["length"],
          pars["radius_large"]*2+pars["gap"]-pars["length"]*tan(pars["alpha"]),
          0
        ));
        initialWall = false;
        cout << "front wall is lifted." << endl;
        setTime(0);
      }
    }

  private:
    LinearViscoelasticFrictionSpecies *small_, *large_;
    RNG generator;
    map<string,double> pars;
    bool initialWall;
    IntersectionOfWalls* frontWall;
    //DeletionBoundary *db_back, *db_left, *db_right, *db_bottom, *db_top, *db_front;
    double kineticEnergyThreshold_;

    void addSmallParticle(BaseParticle& P, Vec3D pos) {
      double r = pars["radius_small"]*generator.getRandomNumber(
                        pars["dispersity_min"],pars["dispersity_max"]);
      double dx = generator.getRandomNumber(-1,1)*(pars["radius_small"]-r);
      double dz = generator.getRandomNumber(-1,1)*(pars["radius_small"]-r);
      P.setRadius(r);
      P.setSpecies(small_);
      P.setPosition(pos+Vec3D(dx,0,dz));
      particleHandler.copyAndAddObject(P);
    }

    void addLargeParticle(BaseParticle& P, Vec3D pos) {
      double r = pars["radius_large"]*generator.getRandomNumber(
                        pars["dispersity_min"],pars["dispersity_max"]);
      double dx = generator.getRandomNumber(-1,1)*(pars["radius_large"]-r);
      double dz = generator.getRandomNumber(-1,1)*(pars["radius_large"]-r);
      P.setRadius(r);
      P.setSpecies(large_);
      P.setPosition(pos+Vec3D(dx,0,dz));
      particleHandler.copyAndAddObject(P);
    }
};

int main(int argc, char *argv[]) {
  if (argc > 1) {
    PerryComb* problem = new PerryComb(argv[1]);
    argv[1] = argv[0];
    problem->solve(argc-1, argv+1);
    delete problem;
  } else {
    fprintf(stderr, "Usage: %s config-file [options]\n", argv[0]);
    exit(-1);
  }
  return 0;
}
