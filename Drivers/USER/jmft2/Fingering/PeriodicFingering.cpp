/* PeriodicFingering - fingering experiments in a periodic domain. 
 * Particles are released by a reservoir and the depth is controlled by a gate. 
 * For best results, take quite a low coefficient of restitution and quite a
 * high coefficient of friction. This helps impose a no-slip condition.
 */
#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/PolydisperseInsertionBoundary.h"
#include "Math/RNG.h"
#include <iostream>
#include <fstream>
#include <map>

/* Specify angles, coefficients of friction in terms of degrees? */
#define DEGREES

/* Do we write .data. and .fstat. files? while filling up? */
#define FILESWHILEFILLING

using namespace std;

void terminationUncaughtException()
{
    cout << "Uncaught exception!" << endl;
    cout << "Perhaps your configuration file doesn't specify all the parameters?" << endl;
    exit(1);
}

class PeriodicFingering : public Mercury3D { 
  public:

    PeriodicFingering(string parsfile) {
      /* Reading config file */
      ifstream file(parsfile);
      string name;
      double var;
      while (file >> name >> var) {
        pars[name] = var;
      }
      /*
      if (pars.size()!=24) {
        logger(ERROR, ".config file not read correctly.");
      }
      */

#ifdef DEGREES
      /* The slope is specified in degrees, not radians. */
      pars.at("alpha") = M_PI * pars.at("alpha") / 180.;
      /* The coefficients of friction are best specified by their arctangents in degrees. */
      pars.at("sliding_small") = tan(M_PI * pars.at("sliding_small") / 180.);
      pars.at("sliding_large") = tan(M_PI * pars.at("sliding_large") / 180.);
      pars.at("rolling_small") = tan(M_PI * pars.at("rolling_small") / 180.);
      pars.at("rolling_large") = tan(M_PI * pars.at("rolling_large") / 180.);
#endif

      /* Initialization */
      setName(parsfile.erase(parsfile.find_last_of('.')));

      if (pars.find("randomSeed") != pars.end())
      {
          random.setRandomSeed(int(pars.at("randomSeed")));
          logger(INFO, "Set random seed to %", pars.at("randomSeed"));
      }
      else
      {
          random.randomise();
          logger(INFO, "randomSeed not specified, randomising.");
      }


      /* Do we write .data. and .fstat. files? while filling up? */
#ifdef FILESWHILEFILLING
      dataFile.setFileType(FileType::MULTIPLE_FILES);
      fStatFile.setFileType(FileType::MULTIPLE_FILES);
#else
      dataFile.setFileType(FileType::NO_FILE);
      fStatFile.setFileType(FileType::NO_FILE);
      eneFile.setFileType(FileType::NO_FILE);
#endif


      setXMin(pars.at("reservoirLength"));
      setXMax(pars.at("domainLength"));
      setYMin(0);
      setYMax(+pars.at("domainWidth"));
      setZMin(0);
      setZMax(pars.at("reservoirHeight")); // This is about to be overwritten...
      setNumberOfDomains({8,2,2});   //For hopper

      /* JMFT: We must set up species in the constructor, or in any case before
       * we call setupInitialConditions(), if we want to use MPI. 
       * However, this can cause problems when restarting...
       * TODO: Talk to Marnix about this.
       */
      /* Set up species for smaller particles */
      small_ = new LinearViscoelasticFrictionSpecies();
      small_->setDensity(pars.at("rho"));
      small_->setCollisionTimeAndRestitutionCoefficient(
        pars.at("collisionTime"),
        pars.at("restitutionCoeff"), 
        (4.0/3.0)*constants::pi*pow(pars.at("radius_small"),3)*pars.at("rho")
      );
      small_->setSlidingFrictionCoefficient(pars.at("sliding_small"));
      small_->setSlidingStiffness(2.0/7.0 * small_->getStiffness());
      small_->setSlidingDissipation(2.0/7.0 * small_->getDissipation());
      small_->setRollingFrictionCoefficient(pars.at("rolling_small"));
      small_->setRollingStiffness(2.0/5.0 * small_->getStiffness());
      small_->setRollingDissipation(2.0/5.0 * small_->getDissipation());
      //small_->setTorsionFrictionCoefficient(pars.at("friction_small"));
      //small_->setTorsionStiffness(2.0/5.0 * small_->getStiffness());
      //small_->setTorsionDissipation(2.0/5.0 * small_->getDissipation());
      small_ = speciesHandler.copyAndAddObject(small_);
      logger(INFO, "Generated small species, maximum collision velocity %", 
              small_->getMaximumVelocity(
                  pars.at("radius_small"),
                  (4.0/3.0)*constants::pi*pow(pars.at("radius_small"),3)*pars.at("rho")));

      /* Set up species for larger particles */
      large_ = new LinearViscoelasticFrictionSpecies();
      large_->setDensity(pars.at("rho"));
      large_->setCollisionTimeAndRestitutionCoefficient(
        pars.at("collisionTime"),
        pars.at("restitutionCoeff"), 
        (4.0/3.0)*constants::pi*pow(pars.at("radius_large"),3)*pars.at("rho")
      );
      large_->setSlidingFrictionCoefficient(pars.at("sliding_large"));
      large_->setSlidingStiffness(2.0/7.0 * large_->getStiffness());
      large_->setSlidingDissipation(2.0/7.0 * large_->getDissipation());
      large_->setRollingFrictionCoefficient(pars.at("rolling_large"));
      large_->setRollingStiffness(2.0/5.0 * large_->getStiffness());
      large_->setRollingDissipation(2.0/5.0 * large_->getDissipation());
      //large_->setTorsionFrictionCoefficient(pars.at("friction_large"));
      //large_->setTorsionStiffness(2.0/5.0 * large_->getStiffness());
      //large_->setTorsionDissipation(2.0/5.0 * large_->getDissipation());
      large_ = speciesHandler.copyAndAddObject(large_);
      logger(INFO, "Generated large species, maximum collision velocity %", 
              large_->getMaximumVelocity(
                  pars.at("radius_large"),
                  (4.0/3.0)*constants::pi*pow(pars.at("radius_large"),3)*pars.at("rho")));

      /* The basal species will have the same properties as the larger species
       * (but should it?) 
       * We make it a separate species so that basal particles can be identified
       * in the .data. file (otherwise you can tell from their velocity which is
       * exactly zero. */
      basal_ = new LinearViscoelasticFrictionSpecies();
      basal_->setDensity(pars.at("rho"));
      basal_->setCollisionTimeAndRestitutionCoefficient(
        pars.at("collisionTime"),
        pars.at("restitutionCoeff"), 
        (4.0/3.0)*constants::pi*pow(pars.at("radius_large"),3)*pars.at("rho")
      );
      basal_->setSlidingFrictionCoefficient(pars.at("sliding_large"));
      basal_->setSlidingStiffness(2.0/7.0 * basal_->getStiffness());
      basal_->setSlidingDissipation(2.0/7.0 * basal_->getDissipation());
      basal_->setRollingFrictionCoefficient(pars.at("rolling_large"));
      basal_->setRollingStiffness(2.0/5.0 * basal_->getStiffness());
      basal_->setRollingDissipation(2.0/5.0 * basal_->getDissipation());
      //basal_->setTorsionFrictionCoefficient(pars.at("friction_large"));
      //basal_->setTorsionStiffness(2.0/5.0 * basal_->getStiffness());
      //basal_->setTorsionDissipation(2.0/5.0 * basal_->getDissipation());
      basal_ = speciesHandler.copyAndAddObject(basal_);

      /* Prototypical particles, to be copied out later */
      smallPrototype = new BaseParticle();
      smallPrototype->setSpecies(small_);
      smallPrototype->setRadius(pars.at("radius_small"));
      largePrototype = new BaseParticle();
      largePrototype->setSpecies(large_);
      largePrototype->setRadius(pars.at("radius_large"));

    }

    ~PeriodicFingering()
    {
        delete smallPrototype;
        delete largePrototype;
    }

    void setupInitialConditions() {

      setTimeStep(pars.at("timeStep"));
      setTimeMax(pars.at("timeMax"));
      setSaveCount(int(pars.at("saveCount")));
      setSystemDimensions(3);

      // setGravity(Vec3D(0,0,0));
      /*
      setGravity(Vec3D(
                  pars.at("g")*sin(pars.at("alpha")), 
                  0, 
                  -pars.at("g")*cos(pars.at("alpha")) 
                  ));
                  */
      setGravity(Vec3D(0,0,-pars.at("g")));

      setXBallsAdditionalArguments("-noborder 4 -v0 -solidf");
      /* Set up walls */
      auto backWall = new InfiniteWall();
      backWall->set(Vec3D(-1,0,0), Vec3D(-pars.at("reservoirLength"),0,0));
      backWall->setSpecies(basal_);
      backWall = wallHandler.copyAndAddObject(backWall);
      auto bottomWall = new InfiniteWall();
      bottomWall->set(Vec3D(0,0,-1),Vec3D(0,0,0));
      bottomWall->setSpecies(basal_);
      bottomWall = wallHandler.copyAndAddObject(bottomWall);

      /* Set up initial confinement */
      // Note, liftableGate is declared as a member of the driver class (so that it
      // may be referred to in actionsAfterTimeStep).
      liftableGate = new IntersectionOfWalls();
      liftableGate->addObject(Vec3D(1,0,0),Vec3D(0,0,0));
      liftableGate->setSpecies(basal_);
      liftableGate = wallHandler.copyAndAddObject(liftableGate);

      auto endWall = new InfiniteWall();
      endWall->set(Vec3D(+1,0,0), Vec3D(pars.at("domainLength"), 0, 0));
      endWall->setSpecies(basal_);
      endWall = wallHandler.copyAndAddObject(endWall);

      /* Sides */
      auto sides = new PeriodicBoundary();
      sides->set(Vec3D(0,1,0), 0, +pars.at("domainWidth"));
      sides = boundaryHandler.copyAndAddObject(sides);

      auto topWall = new InfiniteWall();
      topWall->set(Vec3D(0,0,+1), Vec3D(0,0,pars.at("reservoirHeight")));
      topWall->setSpecies(basal_);
      topWall = wallHandler.copyAndAddObject(topWall);

      /* A PolydisperseInsertionBoundary */
      // JMFT: TODO: You have to add the PolydisperseInsertionBoundary to the
      // boundaryHandler before setting its properties. The other way round
      // doesn't work. That's because the copy constructor of the
      // PolydisperseInsertionBoundary currently (as of r2978) doesn't copy
      // species over.
      auto insb = boundaryHandler.copyAndAddObject(new
              PolydisperseInsertionBoundary());
      insb->setGeometry(pars.at("insbMaxFailed"),
              Vec3D( -pars.at("reservoirLength"), 0, 0.9*pars.at("reservoirHeight") ),
              Vec3D( 0, pars.at("domainWidth"), pars.at("reservoirHeight") ),
              Vec3D(0, 0, - sqrt(pars.at("g") * pars.at("reservoirHeight"))),
              Vec3D(0, 0, - sqrt(pars.at("g") * pars.at("reservoirHeight")))
              );

      double probsmall = pars.at("ratio_small"); // TODO
      double problarge = 1 - probsmall;
      std::cout << "probsmall = " << probsmall << std::endl;
      insb->addGenerandum(smallPrototype, probsmall, pars.at("dispersitySmall") );
      insb->addGenerandum(largePrototype, problarge, pars.at("dispersityLarge") );
      insb->checkBoundaryBeforeTimeStep(this);
      stillFillingUp = true;
    }

    void actionsOnRestart() 
    {
        /* Set the pointer for liftableGate */
        liftableGate = (IntersectionOfWalls*)wallHandler.getObjectById(2);
        /* Check whether we have already lifted the gate or not, by examining
         * the number of objects in the IntersectionOfWalls. */
        switch(liftableGate->getNumberOfObjects())
        {
            case 1:
                stillFillingUp = true;
                break;
            case 3:
                stillFillingUp = false;
                break;
            default:
                logger(ERROR, "liftableGate has % objects, which shouldn't happen. (Should be 1, if still filling up, or 3 if already filled.) Terminating.", 
                        liftableGate->getNumberOfObjects());
                break;
        }

        if (stillFillingUp)
        {
#ifdef FILESWHILEFILLING
            eneFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::MULTIPLE_FILES);
            dataFile.setFileType(FileType::MULTIPLE_FILES);
#else
            eneFile.setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::NO_FILE);
            dataFile.setFileType(FileType::NO_FILE);
#endif
        }
        else
        {
            eneFile.setFileType(FileType::ONE_FILE);

            if (pars.at("fStatFile")==0) {
                fStatFile.setFileType(FileType::NO_FILE);
                logger(WARN, ".fstat file will not be generated.");
            }
            else
                fStatFile.setFileType(FileType::MULTIPLE_FILES);

            if (pars.at("dataFile")==0) {
                dataFile.setFileType(FileType::NO_FILE);
                logger(WARN, ".data file will not be generated.");
            }
            else
                dataFile.setFileType(FileType::MULTIPLE_FILES);
        }
    }

    void printTime() const override
    {
        Mdouble rmsvel = sqrt(2*getKineticEnergy() / getTotalMass());
        Mdouble volfrac = particleHandler.getVolume() 
                            / ( pars.at("reservoirLength") * pars.at("domainWidth") * pars.at("reservoirHeight") );

        logger(INFO, 
                "t = %, tmax = %, Np = %, m = %, KE = %, rmsvel = %, vol = %, volfrac = %",
                getTime(), getTimeMax(), 
                 particleHandler.getNumberOfObjects(),
                 getTotalMass(),
                 getKineticEnergy(),
                 rmsvel,
                 particleHandler.getVolume(),
                 volfrac);

        /* Warnings about particles overlapping (try cranking up stiffness?) */
        int warningsSoFar = 0;
        Mdouble worstOverlapRatio = 0;
        Vec3D   worstOverlapPosition = Vec3D(0,0,0);
        for (auto i : interactionHandler)
        {
            Mdouble overlap = i->getOverlap();
            Mdouble overlapRatio = overlap / pars.at("radius_small");
            if (overlapRatio > worstOverlapRatio) 
            {
                worstOverlapRatio = overlapRatio;
                worstOverlapPosition = i->getContactPoint();
            }
            if (overlapRatio > 0.05)
            {
                if (warningsSoFar < 5)
                {
                    logger(WARN, "interaction at % has overlap %, overlap/radius = %, relative vel %",
                            i->getContactPoint(), overlap, overlapRatio,
                            i->getNormalRelativeVelocity());
                    warningsSoFar++;
                    if (warningsSoFar >= 5)
                        logger(WARN, "further overlap warnings this timestep will not be printed");
                }
            }
        }
        if (worstOverlapRatio > 0.05 && warningsSoFar > 5)
            logger(WARN, "worst overlap/radius = % at position %",
                    worstOverlapRatio, worstOverlapPosition);

        std::cout.flush();
    }

    void actionsAfterTimeStep() 
    {

        /* If we're done filling up, no need to do anything here, just let
         * things run. */
        if (!stillFillingUp) 
            return;

        /* Otherwise... Have we met the conditions for the system to stop
         * filling up, and start releasing? */

        /* We shall say that the system has come to rest provided that the
         * particles' r.m.s. speed is sufficiently low and the volume fraction
         * exceeds a certain amount. */

        /* There are two separate conditions:
         *   After the volume fraction is high enough, remove the
         *   InsertionBoundary.
         *
         *   After the particles' r.m.s. speed is low enough, open the gate.
         * These are tested for separately.
         */

        Mdouble rmsvelThreshold = pars.at("rmsvelThreshold");
        Mdouble volfracThreshold = pars.at("volfracThreshold");

        // Volume fraction in the reservoir
        Mdouble volfrac = particleHandler.getVolume() 
                            / ( pars.at("reservoirLength") * pars.at("domainWidth") * pars.at("reservoirHeight") );

        if ( volfrac > volfracThreshold && boundaryHandler.getNumberOfObjects() > 1 )
        {
            /* We have introduced enough particles. */

            /* Get rid of all the InsertionBoundary that we used to introduce
             * particles. Be careful --- we don't want to delete the
             * PeriodicBoundary! */
            boundaryHandler.clear();
            /* Put the periodic sides back in */
            auto sides = new PeriodicBoundary();
            sides->set(Vec3D(0,1,0), 0, +pars.at("domainWidth"));
            sides = boundaryHandler.copyAndAddObject(sides);

            /* Turn on gravity */
            setGravity(Vec3D(
                        pars.at("g")*sin(pars.at("alpha")), 
                        0, 
                        -pars.at("g")*cos(pars.at("alpha")) 
                        ));
            /* Give particles some downwards impulse so that the rmsvelThreshold
             * condition isn't broken first */
            for (auto p : particleHandler)
                p->setVelocity(Vec3D(0,0,-rmsvelThreshold * 1.1));

            logger(INFO, 
                    "t = %, removed all the InsertionBoundary. volume = %, volfrac = %",
                    getTime(), particleHandler.getVolume(), volfrac);

            writeRestartFile();
            setTime(0);
        }
        
        Mdouble rmsvel = sqrt(2*getKineticEnergy() / getTotalMass());

        if ( rmsvel < rmsvelThreshold && volfrac > volfracThreshold 
//                && boundaryHandler.getNumberOfObjects() == 1
           )
        {
            /* We have allowed the particles to settle. */

            /* Get rid of all residual velocity. (Not really that important? But
             * useful for comparison between experiments?) */
            for (auto p : particleHandler)
                p->setVelocity(Vec3D(0,0,0));

            /* Remove initial confinement after the particles come to rest */
            liftableGate->addObject(Vec3D(0,0,1),
                    Vec3D(0,0, pars.at("gateHeight")));
            liftableGate->addObject(Vec3D(-1,0,0),
                    Vec3D(2*pars.at("radius_large"),0,0));
            stillFillingUp = false;
            cout << "front wall is lifted." << endl;


            /* Set up particles for rough base. (Don't set them up before the
             * initial collapse, or you waste time calculating them (although they
             * are static so it's fairly cheap.) */
            double basalConcentration;
            if (pars.find("basalConcentration") != pars.end())
                basalConcentration = pars.at("basalConcentration");
            else
            {
                basalConcentration = 1;
                std::cout << "Set basalConcentration to default value 1" << std::endl;
            }
            int NBasalParticle = basalConcentration * pars.at("domainWidth") * pars.at("domainLength") 
                / (4 * pow(pars.at("radius_large"), 2));


            /* The rough base should be generated according to a seed, so that
             * we can make the same rough base across different experiments
             * (assuming same domain length and width) */
            if (pars.find("basalSeed") != pars.end())
            {
                baseGenerator.setRandomSeed(pars.at("basalSeed"));
                std::cout << "Set the random seed for the rough base to be " << pars.at("basalSeed") << std::endl;
            }
            else
                std::cout << "basalSeed not set: didn't change the random seed for the rough base." << std::endl;

            for (int n = 0; n < NBasalParticle; n++) 
            {
                // double rand = baseGenerator.getRandomNumber(0,1);
                // double r = rand*pars.at("radius_small") + (1-rand)*pars.at("radius_large");
                double r = pars.at("radius_large");
                double x = baseGenerator.getRandomNumber(0,1) * pars.at("domainLength");
                double y = baseGenerator.getRandomNumber(0,1) * pars.at("domainWidth");
                double z = 0;

                BaseParticle P;
                P.fixParticle();
                P.setRadius(r);
                P.setSpecies(basal_);
                P.setPosition(Vec3D(x,y,z));
                particleHandler.copyAndAddObject(P);
            }
            cout << NBasalParticle << " base particles are generated." << endl;


            /* Write the state */
            // Change the name of the problem first, so that the written
            // .restart file doesn't get overwritten. 
            auto name = getName();
            setName( name + ".layered" );
            writeRestartFile();
            // restore the original name
            setName(name); 

            /* The experiment has properly started! Start writing to .data. and
             * .fstat. files, if the user wants them. */
            setTime(0);
            eneFile.setFileType(FileType::ONE_FILE);

            if (pars.at("fStatFile")==0) {
                fStatFile.setFileType(FileType::NO_FILE);
                logger(WARN, ".fstat file will not be generated.");
            }
            else
                fStatFile.setFileType(FileType::MULTIPLE_FILES);

            if (pars.at("dataFile")==0) {
                dataFile.setFileType(FileType::NO_FILE);
                logger(WARN, ".data file will not be generated.");
            }
            else
                dataFile.setFileType(FileType::MULTIPLE_FILES);

#ifdef FILESWHILEFILLING
            dataFile.setCounter(0);
            fStatFile.setCounter(0);
#endif


            forceWriteOutputFiles(); 
        }

    }

  private:
    LinearViscoelasticFrictionSpecies *small_, *large_, *basal_;
    RNG baseGenerator;
    map<string,double> pars;
    bool stillFillingUp;
    IntersectionOfWalls* liftableGate;
    BaseParticle *smallPrototype, *largePrototype;

};

int main(int argc, char *argv[]) {
    // set_terminate (terminationUncaughtException);
    if (argc > 1) {
        PeriodicFingering* problem = new PeriodicFingering(argv[1]);
        argv[1] = argv[0];
        problem->solve(argc-1, argv+1);
        delete problem;
    } else {
        fprintf(stderr, "Usage: %s config-file [options]\n", argv[0]);
        exit(-1);
    }
    return 0;
}
