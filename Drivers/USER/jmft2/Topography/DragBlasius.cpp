/* DragBlasius - As with MaserBlasius, but now the particles inside the Maser
 * experience an additional drag force that is proportional to their velocity,
 * hopefully letting them move at a controlled velocity. */
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/SubcriticalMaserBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <map>

class DragBlasius : public Mercury2D {
    public:

        /* Constructor */
        DragBlasius(std::string parsfile)
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
                logger(INFO, "randomSeed not specified, randomising.");
            }

            /* roughBaseType:   0 for sudden 'step function', 
             *                  1 for mollified 'rising tanh' (default), 
             *                  2 for bubble 'expanding tanh'
             * transitionLengthscale
             */
            if (pars.find("roughBaseType") == pars.end())
            {
                pars["roughBaseType"] = 1; 
                logger(INFO, "roughBaseType not specified, setting to 1 (for a rising tanh)");
            }
            else
                logger(INFO, "roughBaseType is set to %", pars.at("roughBaseType"));

            /* We want to use a switch, so we need to switch over an int,, not a
             * double. */
            roughBaseType = (int) pars.at("roughBaseType");

            if (pars.find("transitionLengthscale") == pars.end())
            {
                pars["transitionLengthscale"] = 0;
                logger(INFO, "transitionLengthscale not specified, setting to 0");
            }

            /* Initially, we don't want to write .data. or .fstat. files, while
             * we are filling up the Maser. Later we will. */
            dataFile.setFileType(FileType::NO_FILE);
            // dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);
            eneFile.setFileType(FileType::NO_FILE);

            setXMin(pars.at("xmin"));
            setXMax(pars.at("xmax"));
            setYMin(0);
            setYMax(pars.at("reservoirHeight"));
            setZMin(0);
            setZMax(1);

            /* Species */
            /* JMFT: We must set up species in the constructor, or in any case before
             * we call setupInitialConditions(), if we want to use MPI. 
             */

            spec_particles = new LinearViscoelasticFrictionSpecies();
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

        ~DragBlasius(void) {
        }

        void setupInitialConditions() 
        {
            setTimeStep(pars.at("timeStep"));
            setTimeMax(pars.at("timeMax"));
            setSaveCount(pars.at("saveEvery"));

            /* Gravity. While filling up, let gravity point straight down. */
            setGravity(Vec3D(0, -pars.at("g"), 0));

            /* Walls */

            // The base
            auto base = new InfiniteWall;
            base->setSpecies(spec_particles);
            base->set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            base = wallHandler.copyAndAddObject(base);

            lid = new InfiniteWall();
            lid->setSpecies(spec_particles);
            lid->set(Vec3D(0, +1, 0), Vec3D(0, pars.at("reservoirHeight"), 0));
            lid = wallHandler.copyAndAddObject(lid);

            auto back = new InfiniteWall();
            back->setSpecies(spec_particles);
            back->set(Vec3D(-1, 0, 0), 
                      Vec3D(pars.at("xmin") - pars.at("reservoirLength"), 0, 0)
                     );
            back = wallHandler.copyAndAddObject(back);

            /* Rough base. The type of rough base is determined by the parameter
             * roughBaseType: 
             *      0 - sudden step transition 
             *      1 - mollified 'rising tanh' transition (default)
             *      2 - 'growing bubble' transition
             */
            BaseParticle rbParticle;
            rbParticle.setSpecies(spec_particles);
            if (pars.at("baseConc") > 0)
                for (double xpos = pars.at("xmin"); 
                        xpos <= pars.at("xmax"); 
                        xpos += 4*pars.at("baseRadius") / pars.at("baseConc"))
                {
                    if (pars.at("transitionLengthscale") == 0 && xpos < 0)
                        continue;

                    double ypos;
                    double radius;
                    bool toInsert = false;
                    switch (roughBaseType)
                    {
                        case 0:
                            ypos = 0;
                            radius = pars.at("baseRadius")
                                        * (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1));
                            toInsert = (xpos >= 0);
                            break;

                        case 1:
                            ypos = ((pars.at("transitionLengthscale") == 0) 
                                    ? 0 
                                    : 0.5 * pars.at("baseRadius") * (
                                            (tanh(xpos / pars.at("transitionLengthscale")) - 1) 
                                          )
                                   );
                            radius = pars.at("baseRadius")
                                        * (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1));
                            toInsert = true;
                            break;

                        case 2:
                            ypos = 0;
                            double scaleFactor;
                            scaleFactor = ( (pars.at("transitionLengthscale") == 0) 
                                            ? 1 : tanh(xpos / pars.at("transitionLengthscale"))
                                          );
                            radius = scaleFactor * pars.at("baseRadius")
                                        * (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1));
                            toInsert = (xpos > 0);
                            break;

                        default:
                            logger(ERROR, "roughBaseType should be one of 0, 1, 2");
                            break;
                    }

                    rbParticle.setRadius(pars.at("baseRadius")  *
                            (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1)));

                    rbParticle.setRadius(radius);
                    rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                    rbParticle.setVelocity(Vec3D(0,0,0));
                    rbParticle.fixParticle();

                    if (toInsert)
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
                          // pars.at("reservoirHeight") - 8*pars.at("particleRadius"), 
                          0, 0),
                    Vec3D(pars.at("xmin"), 
                          pars.at("reservoirHeight"), 0),
//                    Vec3D(pars.at("reservoirVel") - velvar, -velvar, 0),
//                    Vec3D(pars.at("reservoirVel") + velvar, +velvar, 0),
                    Vec3D(0, - sqrt(pars.at("g") * pars.at("reservoirHeight")), 0),
                    Vec3D(0, - sqrt(pars.at("g") * pars.at("reservoirHeight")), 0),
                    pars.at("particleRadius") * (1 - pars.at("dispersity")),
                    pars.at("particleRadius") * (1 + pars.at("dispersity"))
                );
            insb = boundaryHandler.copyAndAddObject(insb);
            insb->checkBoundaryBeforeTimeStep(this);

            stillFillingUp = true;

            // dataFile.setFileType(FileType::MULTIPLE_FILES);
            // fStatFile.setFileType(FileType::MULTIPLE_FILES);
            restartFile.setFileType(FileType::ONE_FILE);

        }

        /* If restarting, we need to assign the pointers properly. 
         * This is a little messy but it must be done. */
        void actionsOnRestart()
        {
            if (wallHandler.getNumberOfObjects() == 3)
            {
                lid = (InfiniteWall*) wallHandler.getObjectById(1);

                insb = (CubeInsertionBoundary*) boundaryHandler.getObjectById(0);
                stillFillingUp = true;
            }
            else
            {
                masb = (SubcriticalMaserBoundary*) boundaryHandler.getObjectById(1);
                stillFillingUp = false;
            }

        }

        void printTime() const override
        {
            Mercury2D::printTime();

            Mdouble volfracThreshold = pars.at("volfracThreshold");
            // Volume fraction in the reservoir
            Mdouble volfrac = particleHandler.getVolume() 
                / ( pars.at("reservoirLength") * pars.at("reservoirHeight") );

            logger(INFO, "volfrac %, volfracThreshold %", 
                    volfrac, volfracThreshold);


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
                // if (masb != nullptr && CI->isMaserParticle())
                if (masb != nullptr && masb->isMaserParticle(CI))
                    CI->addForce(Vec3D(
                            - CI->getMass() * getGravity().X * CI->getVelocity().X / pars.at("reservoirVel"), 
                        0, 0));
            }
        }

        void actionsAfterTimeStep() 
        {

            /* Are we still filling up? If not, no need to do anything here. */
            if (!stillFillingUp) 
                return;

            /* Have we filled the reservoir up enough already? */

            // Volume fraction in the reservoir
            Mdouble volfrac = particleHandler.getVolume()
                / ( pars.at("reservoirLength") * pars.at("reservoirHeight") );

            if (volfrac > pars.at("volfracThreshold"))
            {
                /* We have introduced enough particles. Get rid of the
                 * InsertionBoundary, and put in all the other boundaries
                 * (maser and deletion) */
                if (boundaryHandler.getNumberOfObjects() != 1)
                    logger(ERROR, "boundaryHandler should have 1 object, but has %",
                            boundaryHandler.getNumberOfObjects());
                else 
                    boundaryHandler.removeObject(insb->getId());

                /* Deletion boundary */
                auto delb = new DeletionBoundary;
                delb->set(Vec3D(1,0,0), pars.at("xmax"));
                delb = boundaryHandler.copyAndAddObject(delb);

                /* MaserBoundary */
                /* If restarting from a system that already has a Maser, we will
                 * never actually reach this. */
                masb = new SubcriticalMaserBoundary();
                masb->set(Vec3D(1.0, 0.0, 0.0), 
                        pars.at("xmin") - pars.at("reservoirLength") + pars.at("particleRadius"), 
                        pars.at("xmin"));
                masb = boundaryHandler.copyAndAddObject(masb);
                masb->activateMaser();

                /* Get rid of the lid */
                wallHandler.removeObject(lid->getIndex());

                /* Give an impulse to all the particles in the Maser. */
                for (BaseParticle* const p : particleHandler)
                    if (! p->isFixed())
                        p->setVelocity(p->getVelocity() + Vec3D(pars.at("reservoirVel"), 0, 0));

                /* Give a label to their initial y-position (I dunno why, might
                 * be useful for studying mixing) */
                for (BaseParticle* const p : particleHandler)
                    p->setInfo(p->getPosition().Y);

                /* Turn on gravity */
                setGravity(Vec3D(
                            pars.at("g") * sin(pars.at("theta")),
                            -pars.at("g") * cos(pars.at("theta")), 
                            0));

                stillFillingUp = false;

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

                forceWriteOutputFiles();
            }
        }

    private:

        std::map<std::string,double> pars;
        bool stillFillingUp;

        int roughBaseType; 

        CubeInsertionBoundary* insb;
        LinearViscoelasticFrictionSpecies *spec_particles;
        // LinearViscoelasticFrictionSpecies *spec_base;
        SubcriticalMaserBoundary* masb;
        InfiniteWall* lid;


};

int main(const int argc, char* argv[]) {
    if (argc > 1)
    {
        DragBlasius* problem = new DragBlasius(argv[1]);
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
