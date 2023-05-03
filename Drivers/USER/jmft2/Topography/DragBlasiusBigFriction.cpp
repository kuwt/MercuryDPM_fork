/* DragBlasiusBigFriction - As with DragBlasius, but with a big basal friction
 * instead of a bumpy base. The parameters beta and betaroll fix interparticular
 * friction. 
 * There is no basal friction for x < 0.
 * The basal friction for x > 0 is set by baseBeta and baseBetaRoll
 */
#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <map>
#include <valarray>

class DragBlasiusBigFriction : public Mercury2D {
    public:

        /* Constructor */
        DragBlasiusBigFriction(std::string parsfile)
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

            if (pars.find("roughBaseType") != pars.end())
                logger(WARN, "roughBaseType is specified, but shouldn't be for DragBlasiusBigFriction");

            if (pars.find("transitionLengthscale") != pars.end())
                logger(WARN, "transitionLengthscale is specified, but shouldn't be for DragBlasiusBigFriction");

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

            spec_baseLeft = new LinearViscoelasticFrictionSpecies();
            spec_baseLeft->setDensity(pars.at("rho"));
            spec_baseLeft->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"), 
                    pars.at("restitutionCoefficient"),
                    constants::pi*pow(pars.at("particleRadius"),2)*pars.at("rho") // note - mass per unit _area_
                    );

            // spec_baseLeft->setSlidingFrictionCoefficient( std::numeric_limits<double>::infinity() );
            spec_baseLeft->setSlidingFrictionCoefficient( 0 );
            spec_baseLeft->setSlidingStiffness(2.0/7.0 * spec_baseLeft->getStiffness());
            spec_baseLeft->setSlidingDissipation(2.0/7.0 * spec_baseLeft->getDissipation());
            // spec_baseLeft->setRollingFrictionCoefficient( std::numeric_limits<double>::infinity() );
            spec_baseLeft->setRollingFrictionCoefficient( 0 );
            spec_baseLeft->setRollingStiffness(2.0/5.0 * spec_baseLeft->getStiffness());
            spec_baseLeft->setRollingDissipation(2.0/5.0 * spec_baseLeft->getDissipation());
            spec_baseLeft = speciesHandler.copyAndAddObject(spec_baseLeft);

            spec_baseRight = new LinearViscoelasticFrictionSpecies();
            spec_baseRight->setDensity(pars.at("rho"));
            spec_baseRight->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"), 
                    pars.at("restitutionCoefficient"),
                    constants::pi*pow(pars.at("particleRadius"),2)*pars.at("rho") // note - mass per unit _area_
                    );

            // spec_baseRight->setSlidingFrictionCoefficient( std::numeric_limits<double>::infinity() );
            spec_baseRight->setSlidingFrictionCoefficient( tan( pars.at("baseBeta") * M_PI / 180. ) );
            spec_baseRight->setSlidingStiffness(2.0/7.0 * spec_baseRight->getStiffness());
            spec_baseRight->setSlidingDissipation(2.0/7.0 * spec_baseRight->getDissipation());
            // spec_baseRight->setRollingFrictionCoefficient( std::numeric_limits<double>::infinity() );
            spec_baseRight->setRollingFrictionCoefficient( tan( pars.at("baseBetaRoll") * M_PI / 180. ) );
            spec_baseRight->setRollingStiffness(2.0/5.0 * spec_baseRight->getStiffness());
            spec_baseRight->setRollingDissipation(2.0/5.0 * spec_baseRight->getDissipation());
            spec_baseRight = speciesHandler.copyAndAddObject(spec_baseRight);

            /* Prototypical particles */
            auto particlePrototype = new SphericalParticle();
            particlePrototype->setSpecies(spec_particles);
            particlePrototype->setRadius(pars.at("particleRadius"));

            logger(INFO, "Maximum collision speed %",
                    spec_particles->getMaximumVelocity(
                        particlePrototype->getRadius(), particlePrototype->getMass()
                    ));

            logger(INFO, "Constructor completed.");

        }

        ~DragBlasiusBigFriction(void) {
        }

        void setupInitialConditions() override
        {
            setTimeStep(pars.at("timeStep"));
            setTimeMax(pars.at("timeMax"));
            setSaveCount(pars.at("saveEvery"));

            /* Gravity. While filling up, let gravity point straight down. */
            setGravity(Vec3D(0, -pars.at("g"), 0));

            /* Walls */

            // The base
            auto baseLeft = new IntersectionOfWalls;
            baseLeft->setSpecies(spec_baseLeft);
            baseLeft->addObject(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            baseLeft->addObject(Vec3D(-1, 0, 0), Vec3D(0, 0, 0));
            baseLeft = wallHandler.copyAndAddObject(baseLeft);

            auto baseRight = new IntersectionOfWalls;
            baseRight->setSpecies(spec_baseRight);
            baseRight->addObject(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            baseRight->addObject(Vec3D(+1, 0, 0), Vec3D(0, 0, 0));
            baseRight = wallHandler.copyAndAddObject(baseRight);

            lid = new InfiniteWall();
            lid->setSpecies(spec_particles);
            lid->set(Vec3D(0, +1, 0), Vec3D(0, pars.at("reservoirHeight"), 0));
            lid = wallHandler.copyAndAddObject(lid);

            auto back = new InfiniteWall();
            back->setSpecies(spec_particles);
            back->set(Vec3D(-1, 0, 0), 
                    Vec3D(pars.at("xmin") - pars.at("reservoirLength") - 7*pars.at("particleRadius"), 0, 0)
                    );
            back = wallHandler.copyAndAddObject(back);

            /* CubeInsertionBoundary for introducing new particles */
            auto generandum = new SphericalParticle;
            generandum->setSpecies(spec_particles);
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
        void actionsOnRestart() override
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

                /* Deletion boundary */
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

        CubeInsertionBoundary* insb;
        LinearViscoelasticFrictionSpecies *spec_particles;
        LinearViscoelasticFrictionSpecies *spec_baseLeft;
        LinearViscoelasticFrictionSpecies *spec_baseRight;
        ConstantMassFlowMaserBoundary* masb;
        InfiniteWall* lid;


};

int main(const int argc, char* argv[]) {
    if (argc > 1)
    {
        DragBlasiusBigFriction* problem = new DragBlasiusBigFriction(argv[1]);
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
