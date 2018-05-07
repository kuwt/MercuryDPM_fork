/* BrokenMaserMinimalExampleBlasius - A simplified form of JMFT's Blasius
 * drivers that doesn't work in parallel.
 *
 * Example usage: mpirun -n 2 ./BrokenMaserMinimalExampleBlasius db-minimal.pars 
 */
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/SubcriticalMaserBoundaryTEST.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <map>

class BrokenMaserMinimalExample : public Mercury2D {
    public:

        /* Constructor */
        BrokenMaserMinimalExample(std::string parsfile)
        {
            /* Reading config file */
            std::ifstream file(parsfile);
            std::string name;
            double var;
            while (file >> name >> var) 
                pars[name] = var;


            /* The slope angle is specified in degrees, not radians. */
            if (pars.find("theta") != pars.end())
                pars.at("theta") = M_PI * pars.at("theta") / 180.;
            else
            {
                pars["theta"] = 0;
                logger(INFO, "slope theta not specified, setting to zero");
            }

            /* Initialisation */
            setName(parsfile.erase(parsfile.find_last_of('.')));
            if (pars.find("randomSeed") != pars.end())
            {
                random.setRandomSeed(int(pars.at("randomSeed")));
                logger(INFO, "Set random seed to %", pars.at("randomSeed"));
            }
            else
            {
                random.randomise();
                logger(INFO, "Random seed not specified, randomising.");
            }

            dataFile.setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::NO_FILE);
            eneFile.setFileType(FileType::NO_FILE);

            setXMin(pars.at("xmin"));
            setXMax(pars.at("xmax"));
            setYMin(0);
            setYMax(pars.at("reservoirHeight"));
            setZMin(0);
            setZMax(1);
            setNumberOfDomains({2, 1, 1});

            /* Species */
            /* JMFT: We must set up species in the constructor, or in any case before
             * we call setupInitialConditions(), if we want to use MPI. 
             */

            spec_particles = new LinearViscoelasticFrictionSpecies();
            spec_particles->setDensity(1);
            spec_particles->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"), 
                    pars.at("restitutionCoefficient"),
                    constants::pi*pow(pars.at("particleRadius"),2) // note - mass per unit _area_
                    );

            spec_particles->setSlidingFrictionCoefficient(tan( pars.at("beta") * M_PI / 180. ));
            spec_particles->setSlidingStiffness(2.0/7.0 * spec_particles->getStiffness());
            spec_particles->setSlidingDissipation(2.0/7.0 * spec_particles->getDissipation());
            spec_particles->setRollingFrictionCoefficient(tan( pars.at("betaroll") * M_PI / 180. ));
            spec_particles->setRollingStiffness(2.0/5.0 * spec_particles->getStiffness());
            spec_particles->setRollingDissipation(2.0/5.0 * spec_particles->getDissipation());
            spec_particles = speciesHandler.copyAndAddObject(spec_particles);

            /* Prototypical particles */
            auto particlePrototype = new BaseParticle();
            particlePrototype->setSpecies(spec_particles);
            particlePrototype->setRadius(pars.at("particleRadius"));

            logger(INFO, "Maximum collision speed %",
                    spec_particles->getMaximumVelocity(
                        particlePrototype->getRadius(), particlePrototype->getMass()
                    ));

            logger(INFO, "Constructor completed.");

        }

        ~BrokenMaserMinimalExample(void) {
        }

        void setupInitialConditions() 
        {
            setTimeStep(pars.at("timeStep"));
            setTimeMax(pars.at("timeMax"));
            setSaveCount(pars.at("saveEvery"));

            /* Gravity. While filling up, let gravity point straight down. */
            setGravity(Vec3D(0, -1, 0));

            /* Walls */

            // The base
            auto base = new InfiniteWall;
            base->setSpecies(spec_particles);
            // base->setSpecies(spec_base);
            base->set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            base = wallHandler.copyAndAddObject(base);

            lid = new InfiniteWall();
            // lid->setSpecies(spec_base);
            lid->setSpecies(spec_particles);
            lid->set(Vec3D(0, +1, 0), Vec3D(0, pars.at("reservoirHeight"), 0));
            lid = wallHandler.copyAndAddObject(lid);

            auto back = new InfiniteWall();
            // back->setSpecies(spec_base);
            back->setSpecies(spec_particles);
            back->set(Vec3D(-1, 0, 0), 
                      Vec3D(pars.at("xmin") - 2*pars.at("particleRadius"), 0, 0)
                     );
            back = wallHandler.copyAndAddObject(back);
                for (double xpos = pars.at("xmin"); 
                        xpos <= pars.at("xmax"); 
                        xpos += 4*pars.at("particleRadius"))
                {
                    if (xpos < 0)
                        continue;

                    double ypos = 0;
                    BaseParticle rbParticle;
                    // rbParticle.setSpecies(spec_base);
                    rbParticle.setSpecies(spec_particles);
                    rbParticle.setRadius(pars.at("particleRadius"));

                    rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                    rbParticle.setVelocity(Vec3D(0,0,0));
                    rbParticle.fixParticle();

                    particleHandler.copyAndAddObject(rbParticle);
                }

            /* CubeInsertionBoundary for introducing new particles */
            auto generandum = new BaseParticle;
            generandum->setSpecies(spec_particles);
            generandum->setRadius(pars.at("particleRadius"));
            double velvar = 0;
            insb = new CubeInsertionBoundary();
            insb->set(
                generandum, 0, 
                    Vec3D(pars.at("xmin") - pars.at("reservoirLength"), 
                          0, 0),
                    Vec3D(pars.at("xmin"), 
                          pars.at("reservoirHeight"), 0),
                    Vec3D(0, - sqrt(pars.at("reservoirHeight")), 0),
                    Vec3D(0, - sqrt(pars.at("reservoirHeight")), 0),
                    pars.at("particleRadius"),
                    pars.at("particleRadius")
                );
            insb = boundaryHandler.copyAndAddObject(insb);
            insb->checkBoundaryBeforeTimeStep(this);

            stillFillingUp = true;

            restartFile.setFileType(FileType::ONE_FILE);
        }

        void actionsOnRestart()
        {
        }

        void printTime() const override
        {
            Mercury2D::printTime();
            logger(INFO, "t %, np local %", 
                    getTime(), particleHandler.getNumberOfRealObjectsLocal());
        }

        void computeExternalForces(BaseParticle* CI) override
        {
            if (!CI->isFixed())
            {
                // Applying the force due to gravity
                CI->addForce(getGravity() * CI->getMass());

                // Wall forces
                computeForcesDueToWalls(CI);

                // Force controller if inside Maser
                // if (masb != nullptr && masb->isMaserParticle(CI))
                if (masb != nullptr && CI->isMaserParticle())
                    CI->addForce(Vec3D(
                            - CI->getMass() * getGravity().X * CI->getVelocity().X / pars.at("reservoirVel"), 
                        0, 0));
            }
        }

        void actionsAfterTimeStep() 
        {
            // logger(INFO, "In actionsAfterTimeStep()");

            // logger(INFO, "t = %, np = %", getTime(), particleHandler.getNumberOfRealObjectsLocal());
            /* Are we still filling up? If not, no need to do anything here. */
            if (!stillFillingUp) 
                return;

            /* Have we filled the reservoir up enough already? */

            // Volume fraction in the reservoir
            Mdouble volfrac = particleHandler.getVolume()
                / ( pars.at("reservoirLength") * pars.at("reservoirHeight") );

            if (volfrac > 0.8)
            {
                /* We have introduced enough particles. Get rid of the
                 * InsertionBoundary, and put in all the other boundaries
                 * (maser and deletion) */
                boundaryHandler.removeObject(insb->getId());

                /* Deletion boundary */
                auto delb = new DeletionBoundary;
                delb->set(Vec3D(1,0,0), pars.at("xmax"));
                delb = boundaryHandler.copyAndAddObject(delb);

                /* MaserBoundary */
                masb = new SubcriticalMaserBoundaryTEST();
                masb->set(Vec3D(1.0, 0.0, 0.0), 
                        pars.at("xmin") - pars.at("reservoirLength") + pars.at("particleRadius"), 
                        pars.at("xmin"));
                logger(INFO, "About to put in masb");
                masb = boundaryHandler.copyAndAddObject(masb);
                logger(INFO, "Have put in the masb");
                logger(INFO, "About to activate masb");
                masb->activateMaser();
                logger(INFO, "Have activated masb");

                /* Get rid of the lid */
                wallHandler.removeObject(lid->getIndex());

                /* Give an impulse to all the particles in the Maser. */
                for (BaseParticle* const p : particleHandler)
                    if (! p->isFixed())
                        p->setVelocity(p->getVelocity() + Vec3D(pars.at("reservoirVel"), 0, 0));

                /* Turn on gravity */
                setGravity(Vec3D(
                            sin(pars.at("theta")),
                            -cos(pars.at("theta")), 
                            0));

                stillFillingUp = false;

                dataFile.setFileType(FileType::MULTIPLE_FILES);
                forceWriteOutputFiles();
            }
            // logger(INFO, "Have completed actionsAfterTimeStep");
        }

    private:

        std::map<std::string,double> pars;
        bool stillFillingUp;

        CubeInsertionBoundary* insb;
        LinearViscoelasticFrictionSpecies *spec_particles;
        SubcriticalMaserBoundaryTEST* masb;
        InfiniteWall* lid;

};

int main(const int argc, char* argv[]) {
    if (argc > 1)
    {
        auto problem = new BrokenMaserMinimalExample(argv[1]);
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
