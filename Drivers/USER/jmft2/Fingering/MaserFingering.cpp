/* MaserFingering - fingering experiments in a periodic domain. 
 * Particles are released from a Maser, which is initially closed to produce a
 * segregated inflow. A finite quantity is released. 
 */
#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Mercury3D.h"
#include "Particles/SphericalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/PolydisperseInsertionBoundary.h"
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
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

class MaserFingering : public Mercury3D { 
  public:

    MaserFingering(string parsfile) {
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
      stillFillingUp = false;
      blocker = new InfiniteWall();

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
      smallPrototype = new SphericalParticle();
      smallPrototype->setSpecies(small_);
      smallPrototype->setRadius(pars.at("radius_small"));
      largePrototype = new SphericalParticle();
      largePrototype->setSpecies(large_);
      largePrototype->setRadius(pars.at("radius_large"));

    }

    ~MaserFingering()
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
      backWall = new InfiniteWall();
      backWall->set(Vec3D(-1,0,0), Vec3D(-pars.at("reservoirLength"),0,0));
      backWall->setSpecies(basal_);
      backWall = wallHandler.copyAndAddObject(backWall);
      blocker->set(Vec3D(1,0,0), Vec3D(0, 0, 0));
      blocker->setSpecies(basal_);
      blocker = wallHandler.copyAndAddObject(blocker);

      auto bottomWall = new InfiniteWall();
      bottomWall->set(Vec3D(0,0,-1),Vec3D(0,0,0));
      bottomWall->setSpecies(basal_);
      bottomWall = wallHandler.copyAndAddObject(bottomWall);

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
              Vec3D( -pars.at("reservoirLength"), 0, 0 ),
              Vec3D( -2*pars.at("radius_large"), pars.at("domainWidth"), pars.at("reservoirHeight") ),
              Vec3D(0, 0, - sqrt(pars.at("g") * pars.at("reservoirHeight"))),
              Vec3D(0, 0, - sqrt(pars.at("g") * pars.at("reservoirHeight")))
              );


      double probsmall = pars.at("ratio_small"); // TODO by volume or by number?
      double problarge = 1 - probsmall;
      std::cout << "probsmall = " << probsmall << std::endl;
      insb->addGenerandum(smallPrototype, probsmall, pars.at("dispersitySmall") );
      insb->addGenerandum(largePrototype, problarge, pars.at("dispersityLarge") );
      insb->checkBoundaryBeforeTimeStep(this);
      stillFillingUp = true;

    }

    void actionsOnRestart() override
    {
        logger(WARN, "In actionsOnRestart()");

        stillFillingUp = false;

        Mdouble volfracThreshold = pars.at("volfracThreshold");

        // Volume fraction in the reservoir
        Mdouble volfrac = particleHandler.getVolume() 
            / ( pars.at("reservoirLength") * pars.at("domainWidth") * pars.at("reservoirHeight") );

        if ( volfrac <= volfracThreshold && boundaryHandler.getNumberOfObjects() > 1 )
        {
            stillFillingUp = true;
            backWall = (InfiniteWall*) wallHandler.getObjectById(0);
            blocker = (InfiniteWall*) wallHandler.getObjectById(1);

#ifdef FILESWHILEFILLING
            eneFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::MULTIPLE_FILES);
            dataFile.setFileType(FileType::MULTIPLE_FILES);
#else
            eneFile.setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::NO_FILE);
            dataFile.setFileType(FileType::NO_FILE);
#endif

            logger(WARN, "Read restart file. Still filling up; volfrac = %", volfrac);
            return;
        }
        else
        {
            stillFillingUp = false;
            masb = (ConstantMassFlowMaserBoundary*) boundaryHandler.getObjectById(1);
            masb->activateMaser();
            masb->turnOffCopying();
            logger(WARN, "Read restart file. Done filling up; volfrac = %", volfrac);
        }

        if (!stillFillingUp 
                && particleHandler.getVolume() < pars.at("totalVolume") 
                && getTime() >= pars.at("timeToOpen"))
        {
            masb->turnOnCopying();
            logger(WARN, "Reading restart file. Maser is copying. t = %", getTime());
        }

        if (!stillFillingUp 
                && particleHandler.getVolume() >= pars.at("totalVolume") )
        {
            masb->turnOffCopying();
            logger(WARN, "Reading restart file. Maser is not copying. t = %", getTime());
        }

        if (!stillFillingUp)
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



        /* There are two separate conditions:
         *   After the volume fraction is high enough, remove the
         *   InsertionBoundary.
         *
         *   After the Maser has run for long enough, activate it.
         * These are tested for separately.
         */

        if (stillFillingUp)
        {
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
                logger(INFO, 
                        "t = %, removed all the InsertionBoundary. volume = %, volfrac = %",
                        getTime(), particleHandler.getVolume(), volfrac);


                wallHandler.removeObject(backWall->getIndex());
                wallHandler.removeObject(blocker->getIndex());

                /* Put in the Maser (initially closed) */
                /* We need to put the Maser in *after* inserting particles using the
                 * PeriodicInsertionBoundary. Otherwise, the in-Maser species is not
                 * set properly. */
                masb = new ConstantMassFlowMaserBoundary();
                masb->set(Vec3D(1,0,0),  - pars.at("reservoirLength"), 0);
                masb = boundaryHandler.copyAndAddObject(masb);
                masb->activateMaser();
                masb->turnOffCopying();
                logger(INFO, "Maser inserted and activated, but not copying yet.");


                /* Turn on gravity */
                setGravity(Vec3D(
                            pars.at("g")*sin(pars.at("alpha")), 
                            0, 
                            -pars.at("g")*cos(pars.at("alpha")) 
                            ));
                /* Give particles some downwards impulse so that the rmsvelThreshold
                 * condition isn't broken first */
                /*
                   for (auto p : particleHandler)
                   p->setVelocity(Vec3D(0,0,-rmsvelThreshold * 1.1));
                */

                /* Set up particles for rough base. (Don't set them up during
                 * the initial filling, or you waste time calculating them
                 * (although they are static so it's fairly cheap.) 
                 *
                 * We need to put in the rough base after activating the Maser,
                 * so that the basal particles don't get converted. (TODO would
                 * that be important?) */
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
                    double x = baseGenerator.getRandomNumber(-pars.at("reservoirLength")-masb->getGapSize(), pars.at("domainLength"));
                    double y = baseGenerator.getRandomNumber(0,1) * pars.at("domainWidth");
                    double z = 0;

                    SphericalParticle P;
                    P.fixParticle();
                    P.setRadius(r);
                    P.setSpecies(basal_);
                    P.setPosition(Vec3D(x,y,z));
                    particleHandler.copyAndAddObject(P);
                }
                cout << NBasalParticle << " base particles are generated." << endl;

                setTime(0);
                forceWriteOutputFiles();

                stillFillingUp = false;
            }
        }

        /* If we're no longer filling up... */
        
        /* Open up the Maser after some time */
        if (!stillFillingUp && (masb->isActivated() && !masb->isCopying()) && particleHandler.getVolume() < pars.at("totalVolume") && getTime() > pars.at("timeToOpen"))
        {
            masb->turnOnCopying();
            logger(INFO, "t=%, maser is now copying", getTime());



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

        /* Deactivate the Maser after enough particles have been released */
        if (!stillFillingUp && (masb->isActivated() && masb->isCopying()) && particleHandler.getVolume() >= pars.at("totalVolume") )
        {
            masb->closeMaser();

        }

    }


  private:
    LinearViscoelasticFrictionSpecies *small_, *large_, *basal_;
    RNG baseGenerator;
    map<string,double> pars;
    bool stillFillingUp;
    InfiniteWall *backWall, *blocker;
    ConstantMassFlowMaserBoundary *masb;
    SphericalParticle *smallPrototype, *largePrototype;

};

int main(int argc, char *argv[]) {
    // set_terminate (terminationUncaughtException);
    if (argc > 1) {
        auto problem = new MaserFingering(argv[1]);
        argv[1] = argv[0];
        problem->solve(argc-1, argv+1);
        delete problem;
    } else {
        fprintf(stderr, "Usage: %s config-file [options]\n", argv[0]);
        exit(-1);
    }
    return 0;
}

