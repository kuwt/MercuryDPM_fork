#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Mercury3D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Math/RNG.h"
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

class Fingering : public Mercury3D { 
  public:

    Fingering(string parsfile) {
      /* Reading config file */
      ifstream file(parsfile);
      string name;
      double var;
      while (file >> name >> var) {
        pars[name] = var;
      }
      if (pars.size()!=24) {
        logger(ERROR, ".config file not read correctly.");
      }

      /* The slope is best specified in degrees, not radians. */
//      pars["alpha"] = M_PI * pars["alpha"] / 180.;
      /* The coefficients of friction are best specified by their arctangents in degrees. */
//      pars["sliding_small"] = tan(M_PI * pars["sliding_small"] / 180.);
//      pars["sliding_large"] = tan(M_PI * pars["sliding_large"] / 180.);
//      pars["rolling_small"] = tan(M_PI * pars["rolling_small"] / 180.);
//      pars["rolling_large"] = tan(M_PI * pars["rolling_large"] / 180.);

      setName(parsfile.erase(parsfile.find_last_of('.')));

      /* Initialization */
      generator.setRandomSeed(int(pars["randomSeed"]));
      if (pars["fStatFile"]==0) {
        fStatFile.setFileType(FileType::NO_FILE);
        logger(WARN, ".fstat file will not be generated.");
      }
      else
        fStatFile.setFileType(FileType::MULTIPLE_FILES);

      if (pars["dataFile"]==0) {
        dataFile.setFileType(FileType::NO_FILE);
        logger(WARN, ".data file will not be generated.");
      }
      else
        dataFile.setFileType(FileType::MULTIPLE_FILES);

      setXMin(0);
      setXMax(pars["length"]*20);
      setYMin(-pars["width"]/2);
      setYMax(+pars["width"]/2);
      setZMin(0);
      setZMax(pars["height"]);
      setNumberOfDomains({8,2,2});   //For hopper

      /* JMFT: We must set up species in the constructor, or in any case before
       * we call setupInitialConditions(), if we want to use MPI. 
       * However, this can cause problems when restarting...
       * TODO: Talk to Marnix about this.
       */
      /* Set up species for smaller particles */
      small_ = new LinearViscoelasticFrictionSpecies();
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
      small_ = speciesHandler.copyAndAddObject(small_);
      logger(INFO, "Generated small species, maximum collision velocity %", 
              small_->getMaximumVelocity(
                  pars["radius_small"],
                  (4.0/3.0)*constants::pi*pow(pars["radius_small"],3)*pars["rho"]));

      /* Set up species for larger particles */
      large_ = new LinearViscoelasticFrictionSpecies();
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
      large_ = speciesHandler.copyAndAddObject(large_);
      logger(INFO, "Generated large species, maximum collision velocity %", 
              large_->getMaximumVelocity(
                  pars["radius_large"],
                  (4.0/3.0)*constants::pi*pow(pars["radius_large"],3)*pars["rho"]));

      /* The basal species will have the same properties as the larger species
       * (but should it?) 
       * We make it a separate species so that basal particles can be identified
       * in the .data. file (otherwise you can tell from their velocity which is
       * exactly zero. */
      basal_ = new LinearViscoelasticFrictionSpecies();
      basal_->setDensity(pars["rho"]);
      basal_->setCollisionTimeAndRestitutionCoefficient(
        pars["collisionTime"],
        pars["restitutionCoeff"], 
        (4.0/3.0)*constants::pi*pow(pars["radius_large"],3)*pars["rho"]
      );
      basal_->setSlidingFrictionCoefficient(pars["sliding_large"]);
      basal_->setSlidingStiffness(2.0/7.0 * basal_->getStiffness());
      basal_->setSlidingDissipation(2.0/7.0 * basal_->getDissipation());
      basal_->setRollingFrictionCoefficient(pars["rolling_large"]);
      basal_->setRollingStiffness(2.0/5.0 * basal_->getStiffness());
      basal_->setRollingDissipation(2.0/5.0 * basal_->getDissipation());
      //basal_->setTorsionFrictionCoefficient(pars["friction_large"]);
      //basal_->setTorsionStiffness(2.0/5.0 * basal_->getStiffness());
      //basal_->setTorsionDissipation(2.0/5.0 * basal_->getDissipation());
      basal_ = speciesHandler.copyAndAddObject(basal_);

      /* Prototypical particles, to be copied out later */
      smallPrototype = new SphericalParticle();
      smallPrototype->setSpecies(small_);
      largePrototype = new SphericalParticle();
      largePrototype->setSpecies(large_);

    }

    ~Fingering()
    {
        delete smallPrototype;
        delete largePrototype;
    }

    void setupInitialConditions() override {

      setTimeStep(pars["timeStep"]);
      setTimeMax(pars["timeMax"]);
      setSaveCount(int(pars["saveCount"]));
      setSystemDimensions(3);

      setGravity(Vec3D(0,0,0));

      setXBallsAdditionalArguments("-noborder 4 -v0 -solidf");

      /* Set up walls */
      auto leftWall = new IntersectionOfWalls();
      leftWall->addObject( Vec3D(0,+1,0), Vec3D(0,+pars["width"]/2,0));
      leftWall->addObject(Vec3D(-1,0,0),Vec3D(0,0,0));
      leftWall->setSpecies(basal_);
      leftWall = wallHandler.copyAndAddObject(leftWall);
      auto rightWall = new IntersectionOfWalls();
      rightWall->addObject(Vec3D(0,-1,0),Vec3D(0,-pars["width"]/2,0));
      rightWall->addObject(Vec3D(-1,0,0),Vec3D(0,0,0));
      rightWall->setSpecies(basal_);
      rightWall = wallHandler.copyAndAddObject(rightWall);
      auto backWall = new InfiniteWall();
      backWall->set(Vec3D(-1,0,0), Vec3D(-pars["length"],0,0));
      backWall->setSpecies(basal_);
      backWall = wallHandler.copyAndAddObject(backWall);
      auto bottomWall = new InfiniteWall();
      bottomWall->set(Vec3D(0,0,-1),Vec3D(0,0,0));
      bottomWall->setSpecies(basal_);
      bottomWall = wallHandler.copyAndAddObject(bottomWall);

      auto topWall = new InfiniteWall();
      topWall->set(Vec3D(0,0,+1), Vec3D(0,0,pars["height"]));
      topWall->setSpecies(basal_);
      topWall = wallHandler.copyAndAddObject(topWall);

      /* Set up initial confinement */
      // Note, frontWall is declared as a member of the driver class (so that it
      // may be referred to in actionsAfterTimeStep).
      frontWall = new IntersectionOfWalls();
      frontWall->addObject(Vec3D(1,0,0),Vec3D(0,0,0));
      frontWall->setSpecies(basal_);
      frontWall = wallHandler.copyAndAddObject(frontWall);
      stillFillingUp = true;

      /* Create a few InsertionBoundary for introducing particles in a layered
       * structure. 
       * We have to distinguish between ratio by number and ratio by
       * volume. What does the parameter "ratio_small" mean?
       * Binbin's original code did ratio by number, but for
       * experimental comparison it's better to use ratio by volume.
       */

      int nlayers = 5;

      insb = new CubeInsertionBoundary();
      /* The InsertionBoundary give particles an initial velocity so that
       * they don't 'block' the InsertionBoundary. */
      for (int layer = 0; layer < nlayers; layer++)
      {
          insb->set(smallPrototype, 1,
                    Vec3D(
                            -pars["length"], -pars["width"] / 2,
                            pars["height"] / nlayers * layer),
                    Vec3D(
                            0, pars["width"] / 2,
                            pars["height"] / nlayers * (layer + pars["ratio_small"])),
                    Vec3D(0.4, 0, 0), Vec3D(0.4, 0, 0));
          boundaryHandler.copyAndAddObject(insb);
          insb->set(largePrototype, 1,
                    Vec3D(
                            -pars["length"], -pars["width"] / 2,
                            pars["height"] / nlayers * (layer + pars["ratio_small"])),
                    Vec3D(
                            0, pars["width"] / 2,
                            pars["height"] / nlayers * (layer + 1)),
                    Vec3D(0.4, 0, 0), Vec3D(0.4, 0, 0));
          boundaryHandler.copyAndAddObject(insb);
      }


    }

    void actionsOnRestart() override
    {
        frontWall = (IntersectionOfWalls*)wallHandler.getObjectById(5);
        switch(frontWall->getNumberOfObjects())
        {
            case 1:
                stillFillingUp = true;
                break;
            case 3:
                stillFillingUp = false;
                break;
            default:
                logger(ERROR, "frontWall has % objects, which shouldn't happen. (Should be 1 or 3.)", 
                        frontWall->getNumberOfObjects());
                break;
        }
    }

    void printTime() 
    {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
              << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
              << std::endl;

        /*
        if (stillFillingUp)
        {
            Mdouble rmsvel = sqrt(2*getKineticEnergy() / getTotalMass());
            Mdouble volfrac = particleHandler.getVolume() 
                / ( pars["length"] * pars["width"] * pars["height"] );
            std::cout 
                      << " np = " << particleHandler.getNumberOfObjects() 
                      << ", m = " << getTotalMass()
                      << ", KE = " << getKineticEnergy() 
                      << ", rmsvel = " <<  rmsvel
                      << ", vol = " << particleHandler.getVolume()
                      << ", volfrac = " << volfrac
                      << std::endl;
        }
        */
        std::cout.flush();
    }

    void actionsAfterTimeStep() 
    {

        /* Are we still filling up? If not, no need to do anything here. */
        if (!stillFillingUp) 
            return;

        /* Have we met the condition for the system to stop filling up, and
         * start releasing? */
        /* We shall say that the system has come to rest provided that the
         * particles' r.m.s. speed is sufficiently low and the volume fraction
         * exceeds a certain amount. */
        Mdouble rmsvel = sqrt(2*getKineticEnergy() / getTotalMass());
        Mdouble volfrac = particleHandler.getVolume() 
                            / ( pars["length"] * pars["width"] * pars["height"] );

        std::cout 
            << " np = " << particleHandler.getNumberOfObjects() 
            << ", m = " << getTotalMass()
            << ", KE = " << getKineticEnergy() 
            << ", rmsvel = " <<  rmsvel
            << ", vol = " << particleHandler.getVolume()
            << ", volfrac = " << volfrac
            << std::endl;

        if ( rmsvel < 0.05 && volfrac > 0.5 ) 
        {

            /* Remove initial confinement after the particles come to rest */
            frontWall->addObject(Vec3D(0,0,1),
                    Vec3D(0,0, pars["gap"]));
            frontWall->addObject(Vec3D(-1,0,0),
                    Vec3D(2*pars["radius_large"],0,0));
            stillFillingUp = false;
            cout << "front wall is lifted." << endl;
            setTime(0);

            /* Set up particles for rough base. (Don't set them up before the
             * initial collapse, or you waste time calculating them (although they
             * are static so it's fairly cheap.) */
            int NBasalParticle = 4 * 2*pars["width"] * 20*pars["length"] 
                / (M_PI * pow(pars["radius_large"], 2));
            for (int n = 0; n < NBasalParticle; n++) 
            {
                double rand = generator.getRandomNumber(0,1);
                double r = rand*pars["radius_small"] + (1-rand)*pars["radius_large"];
                double x = generator.getRandomNumber(0, 20*pars["length"]);
                double y = generator.getRandomNumber(-1,1) * pars["width"];
                double z = 0;

                SphericalParticle P;
                P.fixParticle();
                P.setRadius(r);
                P.setSpecies(basal_);
                P.setPosition(Vec3D(x,y,z));
                particleHandler.copyAndAddObject(P);
            }
            cout << NBasalParticle << " base particles are generated." << endl;

            /* Turn on gravity */
            setGravity(Vec3D(
                        pars["g"]*sin(pars["alpha"]), 0, -pars["g"]*cos(pars["alpha"])
                        ));


            /* Get rid of all InsertionBoundary that we used to introduce particles */
            boundaryHandler.clear();

            /* Write the state */
            // Change the name of the problem first, so that the written
            // .restart file doesn't get overwritten. 
            auto name = getName();
            setName( name + ".layered" );
            writeRestartFile();
            // restore the original name
            setName(name); 

        }

    }

  private:
    LinearViscoelasticFrictionSpecies *small_, *large_, *basal_;
    RNG generator;
    map<string,double> pars;
    bool stillFillingUp;
    IntersectionOfWalls* frontWall;
    CubeInsertionBoundary* insb;
    SphericalParticle *smallPrototype, *largePrototype;

};

int main(int argc, char *argv[]) {
  if (argc > 1) {
    Fingering* problem = new Fingering(argv[1]);
    argv[1] = argv[0];
    problem->solve(argc-1, argv+1);
    delete problem;
  } else {
    fprintf(stderr, "Usage: %s config-file [options]\n", argv[0]);
    exit(-1);
  }
  return 0;
}
