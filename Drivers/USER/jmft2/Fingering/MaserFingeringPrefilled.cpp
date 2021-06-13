/* MaserFingeringPrefilled - fingering experiments in a periodic domain. 
 * Particles are released from a Maser, which is initially closed to produce a
 * segregated inflow. A finite quantity is released.
 * Unlike MaserFinering, 
 */
#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Mercury3D.h"
#include "Particles/SphericalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "Math/RNG.h"

#include <fstream>
#include <iostream>
#include <map>

/* Do we write .data. and .fstat. files? while filling up? */
#define FILESWHILEFILLING

using namespace std;

class MaserFingeringPrefilled : public Mercury3D { 
  public:

      MaserFingeringPrefilled(string parsfile) {
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

          /* The slope is specified in degrees, not radians. */
          pars.at("alpha") = M_PI * pars.at("alpha") / 180.;
          /* The coefficients of friction are actually specified by their
           * arctangents in degrees. */
          pars.at("sliding_small") = tan(M_PI * pars.at("sliding_small") / 180.);
          pars.at("sliding_large") = tan(M_PI * pars.at("sliding_large") / 180.);
          pars.at("rolling_small") = tan(M_PI * pars.at("rolling_small") / 180.);
          pars.at("rolling_large") = tan(M_PI * pars.at("rolling_large") / 180.);

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

          setXMin(-pars.at("reservoirLength"));
          setXMax(pars.at("domainLength"));
          setYMin(0);
          setYMax(+pars.at("domainWidth"));
          setZMin(0);
          setZMax(pars.at("reservoirHeight"));
          setNumberOfDomains({8,2,2});

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
                  (4.0/3.0) * constants::pi * pow(pars.at("radius_large"), 3) * pars.at("rho")
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

    ~MaserFingeringPrefilled()
    {
        delete smallPrototype;
        delete largePrototype;
    }


    /* Only count the z component of gravity in calculating the
     * potential energy, so that we can meaningfully compare this to the
     * kinetic energy and use this to calculate a Froude number.
     */
    Mdouble getGravitationalEnergy() const override {
        Mdouble gravitationalEnergy = 0;
        for (const BaseParticle* const p : particleHandler)
        {
            if (!(p->isFixed()))
            {
                gravitationalEnergy += p->getMass() * Vec3D::dot(
                    Vec3D(0, 0, -getGravity().getZ()), 
                    p->getPosition()
                );
            }
        }
        return gravitationalEnergy;
    }

    void setupInitialConditions() override {
        startedProper = false;

        setTimeStep(pars.at("timeStep"));
        setTimeMax(pars.at("timeMax"));
        setSaveCount(int(pars.at("saveCount")));
        setSystemDimensions(3);

        setXBallsAdditionalArguments("-noborder 4 -v0 -solidf");

        /* Set up walls */

        auto bottomWall = new InfiniteWall();
        bottomWall->set(Vec3D(0,0,-1), Vec3D(0,0,0));
        bottomWall->setSpecies(basal_);
        bottomWall = wallHandler.copyAndAddObject(bottomWall);

        /* Sides */
        auto sides = new PeriodicBoundary();
        sides->set(Vec3D(0,1,0), 0, pars.at("domainWidth"));
        sides = boundaryHandler.copyAndAddObject(sides);

        /* Probability of generating a small particle. Since ratio_small
         * specifies a proportion by volume, we need to convert this to
         * a ratio by number. */
        double r = pars.at("ratio_small") / (1 - pars.at("ratio_small")) 
                    * pow(pars.at("radius_large") / pars.at("radius_small"), 3);
        double prob = r / (1 + r);
        std::cout << "Probability of generating a small particle is " << prob << std::endl;
        /* Space between particles */
        double space = 2.1 * pars.at("radius_large") * (1 + pars.at("dispersityLarge"));
        for (double z = 2 + space; z < pars.at("reservoirHeight") ; z += space) 
            for (double x = -pars.at("reservoirLength") + space; x < -space; x += space)
                for (double y = 0; y < pars.at("domainWidth") - space; y += space) 
                {
                    // std::cout << x << " " << y  << " " << z << std::endl;
                    // std::cout.flush();

                    SphericalParticle* particle;
                    if (baseGenerator.getRandomNumber(0, 1) < prob) {
                        particle = smallPrototype;
                        // TODO dispersity 
                        particle->setRadius(pars.at("radius_small"));  
                    } else {
                        particle = largePrototype;
                        particle->setRadius(pars.at("radius_large"));
                    }

                    particle->setPosition(Vec3D(x, y, z));
                    particle->setVelocity(Vec3D(
                        2 + baseGenerator.getRandomNumber(-0.5, 0.5), 
                        baseGenerator.getRandomNumber(-0.5, 0.5), 
                        -2 + baseGenerator.getRandomNumber(-0.5, 0.5)
                    ));
                    particleHandler.copyAndAddObject(particle);
                }

        logger(INFO, "Done filling");

        /* Put in the Maser */
        logger(INFO, "Putting in Maser boundary");
        masb = new ConstantMassFlowMaserBoundary();
        masb->set(Vec3D(1,0,0), -pars.at("reservoirLength"), 0);
        masb = boundaryHandler.copyAndAddObject(masb);
        masb->activateMaser();
        masb->turnOffCopying();
        logger(INFO, "Maser inserted and activated, but not copying yet.");


        /* Turn on gravity */
        setGravity(Vec3D(
                    pars.at("g") * sin(pars.at("alpha")), 
                    0, 
                    -pars.at("g") * cos(pars.at("alpha")) 
                    ));

        makeRoughBase();
    }

    /* Set up particles for rough base. (Don't set them up during
     * the initial filling, or you waste time calculating them
     * (although they are static so it's fairly cheap.) 
     *
     * We need to put in the rough base after activating the Maser,
     * so that the basal particles don't get converted. (TODO would
     * that be important?) */
    void makeRoughBase() {
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
            double x = baseGenerator.getRandomNumber(
                -pars.at("reservoirLength")-masb->getGapSize(), 
                pars.at("domainLength")
            );
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
    }

    void actionsOnRestart() override
    {
        logger(WARN, "In actionsOnRestart()");

        Mdouble volfracThreshold = pars.at("volfracThreshold");

        // Volume fraction in the reservoir
        Mdouble volfrac = particleHandler.getVolume() 
            / ( pars.at("reservoirLength") * pars.at("domainWidth") * pars.at("reservoirHeight") );

        if ( volfrac <= volfracThreshold && boundaryHandler.getNumberOfObjects() > 1 )
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

            logger(WARN, "Read restart file. Still filling up; volfrac = %", volfrac);
            return;
        }
        else
        {
            stillFillingUp = false;
            masb = (ConstantMassFlowMaserBoundary*) boundaryHandler.getObjectById(1);
            masb->activateMaser();
            logger(WARN, "Read restart file. Done filling up; volfrac = %", volfrac);
        }

        if (
            particleHandler.getVolume() < pars.at("totalVolume") 
            && getTime() >= pars.at("timeToOpen")
        ) {
            masb->turnOnCopying();
            logger(WARN, "Reading restart file. Maser is copying. t = %", getTime());
        }

        if (particleHandler.getVolume() >= pars.at("totalVolume"))
        {
            masb->turnOffCopying();
            logger(WARN, "Reading restart file. Maser is not copying. t = %", getTime());
        }

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

    /* Turn on the Maser's copying, and reset time and indices to 0. */
    void startSimulationProper()
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

        startedProper = true;
    }

    void actionsAfterTimeStep() override
    {
        /* After allowing the Maser to settle for a long enough time, 
         * open it. */

        if (masb->isActivated() 
                && !masb->isCopying()
                && !startedProper
                && getTime() > pars.at("timeToOpen")
        ) {
            startSimulationProper();
        }

        /* Deactivate the Maser after enough particles have been released */
        if (masb->isActivated() 
            && masb->isCopying()
            && startedProper
            && particleHandler.getVolume() >= pars.at("totalVolume")
        ) {
            logger(INFO, "Released enough particles, so closing the Maser");
            masb->closeMaser();
        }
    }

  private:
    LinearViscoelasticFrictionSpecies *small_, *large_, *basal_;
    RNG baseGenerator;
    map<string,double> pars;
    bool stillFillingUp;
    bool startedProper;
    ConstantMassFlowMaserBoundary *masb;
    SphericalParticle *smallPrototype, *largePrototype;

};

int main(int argc, char *argv[]) {
    if (argc > 1) {
        auto problem = new MaserFingeringPrefilled(argv[1]);
        argv[1] = argv[0];
        problem->solve(argc-1, argv+1);
        delete problem;
    } else {
        fprintf(stderr, "Usage: %s config-file [options]\n", argv[0]);
        exit(-1);
    }
    return 0;
}
