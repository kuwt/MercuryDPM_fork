/* TimeDependentBlasius - Periodic domain; initially smooth base, rough base
 * suddenly turned on */
#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <map>
#include <valarray>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"

class TimeDependentBlasius : public Mercury2D {
    public:

        /* Constructor */
        TimeDependentBlasius(std::string parsfile)
        {
            /* Reading config file */
            std::ifstream file(parsfile);
            std::string name;
            double var;
            while (file >> name >> var) 
                pars[name] = var;

            /* The slope angle is specified in degrees, not radians. */
            pars["theta"] = M_PI * pars.at("theta") / 180.;

            /* Initialisation */
            setName(parsfile.erase(parsfile.find_last_of('.')));
            if (pars.find("randomSeed") != pars.end())
            {
                generator.setRandomSeed(int(pars.at("randomSeed")));
                logger(INFO, "Set random seed to %", pars.at("randomSeed"));
            }
            else
            {
                generator.randomise();
                logger(INFO, "Random seed not specified, randomising.");
            }

            /* We don't need any .data files when filling up. */
            eneFile.setFileType(FileType::NO_FILE);
            dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);

            setXMin(0);
            setXMax(pars.at("length"));
            setYMin(0);
            setYMax(pars.at("height"));
            setZMin(0);
            setZMax(1);

            /* Set defaults */
            if (pars.find("rho") == pars.end())
                pars["rho"] = 1;
            if (pars.find("base_rho") == pars.end())
                pars["base_rho"] = pars.at("rho");
            if (pars.find("g") == pars.end())
            {
                pars["g"] = 1;
                logger(WARN, "parameter g not specified, taking g = 1");
            }
            if (pars.find("baseRadius") == pars.end()) 
                pars["baseRadius"] = pars.at("particleRadius");
            if (pars.find("base_collisionTime") == pars.end())
                pars["base_collisionTime"] = pars.at("collisionTime");
            if (pars.find("base_restitutionCoefficient") == pars.end())
                pars["base_restitutionCoefficient"] = pars.at("restitutionCoefficient");
            if (pars.find("base_beta") == pars.end()) 
                pars["base_beta"] = pars.at("beta");
            if (pars.find("base_betaroll") == pars.end())
                pars["base_betaroll"] = pars.at("betaroll");
            if (pars.find("baseConc") == pars.end())
                pars["baseConc"] = 1;
            if (pars.find("baseDispersity") == pars.end())
                pars["baseDispersity"] = 0;

            /* Species */
            /* JMFT: We must set up species in the constructor, or in any case before
             * we call setupInitialConditions(), if we want to use MPI. 
             * However, this can cause problems when restarting...
             * TODO: Talk to Marnix about this.
             */

            auto spec_particles = new LinearViscoelasticFrictionSpecies();
            spec_particles->setDensity(pars.at("rho"));
            spec_particles->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"), 
                    pars.at("restitutionCoefficient"),
                    constants::pi*pow(pars.at("particleRadius"),2)*pars.at("rho") // note - mass per unit _area_
                    );

            spec_particles->setSlidingFrictionCoefficient(tan( pars.at("beta") * M_PI / 180. ));
            spec_particles->setSlidingStiffness(2.0/7.0 * spec_particles->getStiffness());
            spec_particles->setSlidingDissipation(2.0/7.0 * spec_particles->getDissipation());
            spec_particles->setRollingFrictionCoefficient(tan( pars.at("betaroll") * M_PI / 180. ));
            spec_particles->setRollingStiffness(2.0/5.0 * spec_particles->getStiffness());
            spec_particles->setRollingDissipation(2.0/5.0 * spec_particles->getDissipation());
            spec_particles = speciesHandler.copyAndAddObject(spec_particles);

            auto spec_base = new LinearViscoelasticFrictionSpecies();
            spec_base->setDensity(pars.find("base_rho") == pars.end() ? pars.at("base_rho") : pars.at("rho"));
            spec_base->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("base_collisionTime"), 
                    pars.at("base_restitutionCoefficient"), 
                    constants::pi*pow(pars.at("baseRadius"),2)*pars.at("rho") // note - mass per unit _area_
                    );

            spec_base->setSlidingFrictionCoefficient(tan( 
                        pars.at("base_beta") * M_PI / 180. ));
            spec_base->setSlidingStiffness(2.0/7.0 * spec_base->getStiffness());
            spec_base->setSlidingDissipation(2.0/7.0 * spec_base->getDissipation());
            spec_base->setRollingFrictionCoefficient(tan( 
                        pars.at("base_betaroll") * M_PI / 180. ));
            spec_base->setRollingStiffness(2.0/5.0 * spec_base->getStiffness());
            spec_base->setRollingDissipation(2.0/5.0 * spec_base->getDissipation());
            spec_base = speciesHandler.copyAndAddObject(spec_base);

            /* Prototypical particles */
            particlePrototype = new SphericalParticle();
            particlePrototype->setSpecies(spec_particles);
            particlePrototype->setRadius(pars.at("particleRadius"));
            basePrototype = new SphericalParticle();
            basePrototype->setSpecies(spec_base);
            // basePrototype->setSpecies(spec_particles);
            basePrototype->setRadius(pars.at("baseRadius"));

            logger(INFO, "Maximum collision speed %",
                    spec_particles->getMaximumVelocity(
                        particlePrototype->getRadius(), particlePrototype->getMass()
                    ));

            logger(INFO, "Constructor completed.");
        }

        ~TimeDependentBlasius(void) {
            delete particlePrototype;
            delete basePrototype;
        }

        void setupInitialConditions() override
        {
            setTimeStep(pars.at("timeStep"));
            setTimeMax(pars.at("timeMax"));
            setSaveCount(pars.at("saveEvery"));

            /* Gravity. While filling up, let gravity point straight down. */
            setGravity(Vec3D(0, -pars.at("g"), 0));

            /* Periodic boundaries */
            PeriodicBoundary bounds;
            bounds.set(Vec3D(+1,0,0), 0, pars.at("length"));
            boundaryHandler.copyAndAddObject(bounds);

            /* Walls */

            auto spec_particles = speciesHandler.getObject(0);
            auto spec_base = speciesHandler.getObject(1);

            // The base
            auto base = wallHandler.copyAndAddObject(new InfiniteWall);
            base->setSpecies(spec_base);
            base->set(Vec3D(0, -1, 0), Vec3D(0,0,0));

            // A lid that will be removed after the domain is filled
            auto lid = wallHandler.copyAndAddObject(new InfiniteWall);
            lid->setSpecies(spec_base);
            lid->set(Vec3D(0, +1, 0), Vec3D(0, pars.at("height"), 0));

            // A dam that will be removed after the domain is filled
            auto dam = wallHandler.copyAndAddObject(new IntersectionOfWalls);
            dam->setSpecies(spec_base);
            dam->addObject(Vec3D(1, 0, 0), Vec3D(pars.at("length")*0.49, 0, 0));
            dam->addObject(Vec3D(-1, 0, 0), Vec3D(pars.at("length")*0.51, 0, 0));

            /* CubeInsertionBoundary for introducing new particles */
            auto generandum = new SphericalParticle;
            generandum->setSpecies(spec_particles);
            generandum->setRadius(pars.at("particleRadius"));
            auto insb = new CubeInsertionBoundary;
            insb->set(
                    generandum, 40,
                    Vec3D(0, 0, 0),
                    Vec3D(pars.at("length"), pars.at("height"), 0),
                    Vec3D(0, -sqrt(pars.at("g") * pars.at("height")), 0),
                    Vec3D(0, -sqrt(pars.at("g") * pars.at("height")), 0));
            insb = boundaryHandler.copyAndAddObject(insb);
            insb->checkBoundaryBeforeTimeStep(this);
            stillFillingUp = true;
            stillSmooth = true;
            
            restartFile.setFileType(FileType::ONE_FILE);

        }

        void actionsOnRestart() override
        {
        }

        void printTime() const override
        {
            Mercury2D::printTime();
            if (stillFillingUp)
                logger(INFO, "volfrac %, volfracThreshold %, np %",
                        particleHandler.getVolume() / (pars.at("length") * pars.at("height")),
                        pars.at("volfracThreshold"),
                        particleHandler.getNumberOfObjects()
                      );
        }

        void actionsAfterTimeStep() 
        {

            /* After a certain time: 
             *   + stop introducing new particles
             *   + remove the dam
             *   + incline the plane
             */

            if (stillFillingUp && 
                    particleHandler.getVolume() / (pars.at("length") * pars.at("height")) > pars.at("volfracThreshold") )
            {
                /* Get rid of InsertionBoundary and dam and lid */

                boundaryHandler.removeObject(1);
                wallHandler.removeObject(2);
                wallHandler.removeObject(1);
                logger(INFO, "Removed InsertionBoundary and dam and lid");

                /* Turn on gravity */
                setGravity(Vec3D(
                            pars.at("g") * sin(pars.at("theta")),
                            -pars.at("g") * cos(pars.at("theta")), 
                            0));

                setTime(-pars.at("triggerTime"));
                eneFile.setFileType(FileType::ONE_FILE);
                dataFile.setFileType(FileType::MULTIPLE_FILES);
                logger(INFO, "Set gravity");

                stillFillingUp = false;
            }

            if (!stillFillingUp && stillSmooth && getTime() >= 0)
            {
                /* Rough base */
                if (pars.at("baseConc") > 0)
                    for (double xpos = 0; 
                            xpos <= pars.at("length"); 
                            xpos += 4*pars.at("baseRadius") / pars.at("baseConc"))
                    {
                        double ypos = 0;
                        SphericalParticle rbParticle;
                        auto spec_base = speciesHandler.getObject(1);
                        rbParticle.setSpecies(spec_base);
                        rbParticle.setRadius(pars.at("baseRadius") *
                                (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1)));
                        rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                        rbParticle.setVelocity(Vec3D(0,0,0));
                        rbParticle.fixParticle();

                        particleHandler.copyAndAddObject(rbParticle);
                    }

                stillSmooth = false;
                logger(INFO, "Included rough base");
            }

        }


    private:
        RNG generator;
        std::map<std::string, double> pars;
        bool stillFillingUp;
        bool stillSmooth;

        BaseParticle *particlePrototype, *basePrototype;

};

int main(const int argc, char* argv[]) {
    // std::set_terminate(terminationUncaughtException);
    if (argc > 1)
    {
        auto problem = new TimeDependentBlasius(argv[1]);
        argv[1] = argv[0];
        problem->solve(argc-1, argv+1);
        delete problem;
        return 0;
    }
    else
    {
        logger(ERROR, "Usage: % config-file [options]",
                argv[0]);
        exit(-1);
    }
}
