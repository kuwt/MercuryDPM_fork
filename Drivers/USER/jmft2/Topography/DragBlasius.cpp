/* DragBlasius - As with MaserBlasius, but now the particles inside the Maser
 * experience an additional drag force that is proportional to their velocity,
 * hopefully letting them move at a controlled velocity. */
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <map>

#include "sys/times.h"
//#include "sys/vtimes.h"
#include<ctime>
#include <valarray>

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

            /* roughBaseType:   0 for sudden 'step function', 
             *                  1 for mollified 'rising tanh' 
             *                  2 for bubble 'expanding tanh'
             *                          (Very slow, don't use this!)
             *                  3 for 'mollified half' (same as 1 but starting only
             *                          after x = 0) (default) 
             *                  4 for Amalia's 'alternating'  base (abusing
             *                          baseDispersity for the larger particles)
             * transitionLengthscale: For roughBaseType = 1, 2 or 3, you need to
             *                  specify this. (default to 0)
             */
            if (pars.find("roughBaseType") == pars.end())
            {
                pars["roughBaseType"] = 1; 
                logger(INFO, "roughBaseType not specified, setting to 3 (for a rising half-tanh)");
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

            if (pars.find("beta") == pars.end())
            {
                pars["beta"] = 0.0;
                logger(INFO, "beta not specified, setting to 0");
            }

            if (pars.find("betaroll") == pars.end())
            {
                pars["betaroll"] = 0.0;
                logger(INFO, "betaroll not specified, setting to 0");
            }

            if (pars.find("collisionTime") == pars.end())
                logger(ERROR, "collisionTime not specified, needed for LinearViscoelastic contact model");

            if (pars.find("elasticModulus") != pars.end())
                logger(WARN, "elasticModulus specified, but this is not used for the LinearViscoelastic contact model");

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

            speciesP = new LinearViscoelasticFrictionSpecies();
            speciesP->setDensity(pars.at("rho"));
            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"), 
                    pars.at("restitutionCoefficient"),
                    constants::pi*pow(pars.at("particleRadius"),2)*pars.at("rho") // note - mass per unit _area_
                    );

            speciesP->setSlidingFrictionCoefficient(tan( pars.at("beta") * M_PI / 180. ));
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(tan( pars.at("betaroll") * M_PI / 180. ));
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP = speciesHandler.copyAndAddObject(speciesP);

            /* Prototypical particles */
            auto particlePrototype = new SphericalParticle();
            particlePrototype->setSpecies(speciesP);
            particlePrototype->setRadius(pars.at("particleRadius"));

            logger(INFO, "Maximum collision speed %",
                    speciesP->getMaximumVelocity(
                        particlePrototype->getRadius(), particlePrototype->getMass()
                    ));

            /* Profiling */
            /* TODO The filename size limit of 1024 is arbitrary and unnecessary. */
            char profilingFileName[1024];
            if (snprintf(profilingFileName, 1024, "%s.prof", getName().c_str()) > 1024)
                logger(ERROR, "profilingFileName exceeds 1024 bytes");

            profilingFile = fopen(profilingFileName, "a+");
            setbuf(profilingFile, NULL);
            fprintf(profilingFile, "tsim Np treal proj\n");
            startTime_ = std::time(nullptr);

            logger(INFO, "Constructor completed.");

        }

        ~DragBlasius(void) {
            fclose(profilingFile);
        }

        void setupInitialConditions() 
        {
            setTimeStep(pars.at("timeStep"));
            setTimeMax(pars.at("timeMax"));
            setSaveCount(pars.at("saveEvery"));

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

            /* Gravity. While filling up, let gravity point straight down. */
            setGravity(Vec3D(0, -pars.at("g"), 0));

            /* Walls */

            // The base
            auto base = new InfiniteWall;
            base->setSpecies(speciesP);
            base->set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            base = wallHandler.copyAndAddObject(base);

            lid = new InfiniteWall();
            lid->setSpecies(speciesP);
            lid->set(Vec3D(0, +1, 0), Vec3D(0, pars.at("reservoirHeight"), 0));
            lid = wallHandler.copyAndAddObject(lid);

            auto backWall = new InfiniteWall();
            backWall->setSpecies(speciesP);
            backWall->set(Vec3D(-1, 0, 0), 
                    Vec3D(pars.at("xmin") - pars.at("reservoirLength") - 7*pars.at("particleRadius"), 
                        0, 0) );
            backWall = wallHandler.copyAndAddObject(backWall);

            /* Rough base. The type of rough base is determined by the parameter
             * roughBaseType: 
             *      0 - sudden step transition 
             *      1 - mollified 'rising tanh' transition 
             *      2 - 'growing bubble' transition
             *      3 - rising half-tanh
             *      4 - Amalia's alternating base (baseRadius for smaller)
             */
            SphericalParticle rbParticle;
            rbParticle.setSpecies(speciesP);
            if (pars.at("baseConc") > 0)
            {
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

                        case 3:
                            ypos = ((pars.at("transitionLengthscale") == 0) 
                                    ? 0 
                                    : pars.at("baseRadius") * (
                                            (tanh(xpos / pars.at("transitionLengthscale")) - 1) 
                                          )
                                   );
                            radius = pars.at("baseRadius")
                                        * (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1));
                            toInsert = (xpos >= 0);
                            break;

                        case 4:
                            logger(ERROR, "Sorry, roughBaseType = 4 (Amalia's alternating) is not implemented in DragBlasius. Try PEChute!");

                        default:
                            logger(ERROR, "roughBaseType should be one of 0, 1, 2, 3 or 4 (but not 4)");
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
            }

            /* CubeInsertionBoundary for introducing particles when initialising
             * the system (will use Maser during the simulation proper) */
            auto generandum = new SphericalParticle;
            generandum->setSpecies(speciesP);
            generandum->setRadius(pars.at("particleRadius"));
            insb = new CubeInsertionBoundary();
            insb->set(
                    generandum, 0,
                    Vec3D(pars.at("xmin") - pars.at("reservoirLength"),
                          0, 0),
                    Vec3D(pars.at("xmin"),
                          pars.at("reservoirHeight"), 0),
                    Vec3D(0, -sqrt(pars.at("g") * pars.at("reservoirHeight")), 0),
                    Vec3D(0, -sqrt(pars.at("g") * pars.at("reservoirHeight")), 0));
            insb = boundaryHandler.copyAndAddObject(insb);
            insb->checkBoundaryBeforeTimeStep(this);

            stillFillingUp = true;
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
                masb = (ConstantMassFlowMaserBoundary*) boundaryHandler.getObjectById(1);
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

            logger(INFO, "t %, n %", 
                    getTime(), particleHandler.getNumberOfObjects());
            logger(INFO, "volfrac %, volfracThreshold %", 
                    volfrac, volfracThreshold);

            /* Warnings about particles overlapping (try cranking up stiffness?) */
            int warningsSoFar = 0;
            Mdouble worstOverlapRatio = 0;
            Vec3D   worstOverlapPosition = Vec3D(0,0,0);
            for (auto i : interactionHandler)
            {
                Mdouble overlap = i->getOverlap();
                Mdouble overlapRatio = overlap / pars.at("particleRadius");
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

            /* Profiling */
            std::time_t timeSinceStart = std::time(nullptr) - startTime_;
            std::time_t projectedTimeLeft = timeSinceStart * (getTimeMax() - getTime())/getTime();
            fprintf(profilingFile, "%f %d %d %d\n",
                getTime(), particleHandler.getNumberOfObjects(),
                timeSinceStart, 
                projectedTimeLeft
            );

        }

        void computeExternalForces(BaseParticle* CI) override
        {
            if (!CI->isFixed())
            {
                DPMBase::computeExternalForces(CI);

                // Force controller if inside Maser
                if (masb != nullptr && masb->isMaserParticle(CI))
                {
                    Mdouble dragForce = - CI->getMass() * getGravity().X * CI->getVelocity().X / pars.at("reservoirVel"); 
                    // logger(INFO, "applying a drag force to particle id % at position % velocity %",
                    //        CI->getId(), CI->getPosition(), CI->getVelocity());
                    CI->addForce(Vec3D(dragForce, 0, 0));
                }
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
                boundaryHandler.removeObject(insb->getId());

                /* Deletion boundary. Put this in after the initialisation. */
                auto delb = new DeletionBoundary;
                delb->set(Vec3D(1,0,0), pars.at("xmax"));
                delb = boundaryHandler.copyAndAddObject(delb);

                /* MaserBoundary */
                /* If restarting from a system that already has a Maser, we will
                 * never actually reach this. */
                masb = new ConstantMassFlowMaserBoundary();
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
                std::cout << "The experiment has properly started!" << std::endl;
                setTime(0);
                eneFile.setFileType(FileType::ONE_FILE);

                /* default is to save the dataFile */
                if (pars.find("dataFile") != pars.end() && pars.at("dataFile") == 0) {
                    dataFile.setFileType(FileType::NO_FILE);
                    logger(INFO, ".data files will not be saved.");
                }
                else
                {
                    dataFile.setFileType(FileType::MULTIPLE_FILES);
                    logger(INFO, ".data files will be saved.");
                }

                /* default is not to save the fStatFile */
                fStatFile.setFileType(FileType::NO_FILE);
                if (pars.find("fStatFile") != pars.end() && pars.at("fStatFile") == 1) 
                {
                    fStatFile.setFileType(FileType::MULTIPLE_FILES);
                    logger(INFO, ".fstat files will be saved");
                }
                else
                {
                    fStatFile.setFileType(FileType::NO_FILE);
                    logger(INFO, ".fstat file will not be saved.");
                }

                forceWriteOutputFiles();
            }
        }

    private:

        std::map<std::string,double> pars;
        bool stillFillingUp;

        int roughBaseType; 

        CubeInsertionBoundary* insb;
        LinearViscoelasticFrictionSpecies *speciesP;
        ConstantMassFlowMaserBoundary* masb;
        InfiniteWall* lid;

        std::time_t startTime_;
        FILE* profilingFile;
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
