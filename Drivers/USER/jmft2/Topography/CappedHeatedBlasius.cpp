/* CappedHeatedBlasius - As with MaserBlasius but using a ConstantMassFlowMaserBoundary for
 * upstream conditions, instead of a SubcriticalMaserBoundary. */
#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "Boundaries/HeaterBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <valarray>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"

void terminationUncaughtException()
{
    std::cout << "Uncaught exception!" << std::endl;
    std::cout << "Perhaps your configuration file doesn't specify all the parameters?" << std::endl;
    exit(1);
}

class CappedHeatedBlasius : public Mercury2D {
    public:

        /* Constructor */
        CappedHeatedBlasius(std::string parsfile)
        {
            /* Reading config file */
            std::ifstream file(parsfile);
            std::string name;
            double var;
            while (file >> name >> var) 
                pars[name] = var;


            /* The slope angle is specified in degrees, not radians. */
            pars.at("theta") = M_PI * pars.at("theta") / 180.;

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

            /* Initially, we don't want to write .data. or .fstat. files, while
             * we are filling up the Maser. Later we will. */
            // dataFile.setFileType(FileType::NO_FILE);
            dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);
            // eneFile.setFileType(FileType::NO_FILE);

            setXMin(pars.at("xmin"));
            setXMax(pars.at("xmax"));
            setYMin(0);
            setYMax(pars.at("reservoirHeight"));
            setZMin(0);
            setZMax(1);
            // setNumberOfDomains(2,2,2);

            /* Species */
            /* JMFT: We must set up species in the constructor, or in any case before
             * we call setupInitialConditions(), if we want to use MPI. 
             * However, this can cause problems when restarting...
             * TODO: Talk to Marnix about this.
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

            /*
            spec_base = new LinearViscoelasticFrictionSpecies();
            spec_base->setDensity(pars.at("rho"));
            spec_base->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"), 
                    pars.at("restitutionCoefficient"),
                    constants::pi*pow(pars.at("baseRadius"),2)*pars.at("rho") // note - mass per unit _area_
                    );

            spec_base->setSlidingFrictionCoefficient(tan( pars.at("beta") * M_PI / 180. ));
            spec_base->setSlidingStiffness(2.0/7.0 * spec_base->getStiffness());
            spec_base->setSlidingDissipation(2.0/7.0 * spec_base->getDissipation());
            spec_base->setRollingFrictionCoefficient(tan( pars.at("betaroll") * M_PI / 180. ));
            spec_base->setRollingStiffness(2.0/5.0 * spec_base->getStiffness());
            spec_base->setRollingDissipation(2.0/5.0 * spec_base->getDissipation());
            spec_base = speciesHandler.copyAndAddObject(spec_base);
            */

            /* Prototypical particles */
            particlePrototype = new SphericalParticle();
            particlePrototype->setSpecies(spec_particles);
            particlePrototype->setRadius(pars.at("particleRadius"));
            basePrototype = new SphericalParticle();
            // basePrototype->setSpecies(spec_base);
            basePrototype->setSpecies(spec_particles);
            basePrototype->setRadius(pars.at("baseRadius"));

            logger(INFO, "Maximum collision speed %",
                    spec_particles->getMaximumVelocity(
                        particlePrototype->getRadius(), particlePrototype->getMass()
                    ));

            logger(INFO, "Constructor completed.");


        }

        ~CappedHeatedBlasius(void) {
            delete particlePrototype;
            delete basePrototype;
        }

        void setupInitialConditions() override
        {
            setTimeStep(pars.at("timeStep"));
            setTimeMax(pars.at("timeMax"));
            setSaveCount(pars.at("saveEvery"));

            /* Gravity. While filling up, let gravity point straight down. */
            setGravity(Vec3D(0,-pars.at("g"), 0));

            /* Walls */
            // We should not have, or need, a back wall if we are using a Maser.

            // The base
            IntersectionOfWalls base;
            base.setSpecies(spec_particles);
            // base.setSpecies(spec_base);
            base.addObject(Vec3D(0, -1, 0), Vec3D(0,0,0));
            base.addObject(Vec3D(-1, 0, 0), Vec3D(pars.at("xmax"), 0,0));
            wallHandler.copyAndAddObject(base);

            // A dam that will be lifted after the Maser is filled
            /*
            dam = new IntersectionOfWalls();
            // dam->setSpecies(spec_base);
            dam->setSpecies(spec_particles);
            // dam->addObject(Vec3D(+1, 0, 0), Vec3D(pars.at("xmin") - 7*pars.at("particleRadius"), 0, 0));
            dam->addObject(Vec3D(+1, 0, 0), 
                    Vec3D(pars.at("xmin") - 0.5*pars.at("reservoirLength") - pars.at("particleRadius"), 0, 0));
            // dam->addObject(Vec3D(-1, 0, 0), Vec3D(pars.at("xmin") - 9*pars.at("particleRadius"), 0, 0));
            dam->addObject(Vec3D(-1, 0, 0), 
                    Vec3D(pars.at("xmin") - 0.5*pars.at("reservoirLength") + pars.at("particleRadius"), 0, 0));
            dam = wallHandler.copyAndAddObject(dam);
            */
            lid = new IntersectionOfWalls();
            // lid->setSpecies(spec_base);
            lid->setSpecies(spec_particles);
            lid->addObject(Vec3D(0, +1, 0), Vec3D(0, pars.at("reservoirHeight"), 0));
            lid->addObject(Vec3D(-1, 0, 0), Vec3D(pars.at("xmin"), 0, 0));
            lid = wallHandler.copyAndAddObject(lid);
            back = new InfiniteWall();
            // back->setSpecies(spec_base);
            back->setSpecies(spec_particles);
            back->set(Vec3D(-1, 0, 0), 
                      Vec3D(pars.at("xmin") - pars.at("reservoirLength"), 0, 0)
                     );
            back = wallHandler.copyAndAddObject(back);

            /* Lining the top and bottom of the Maser */
            if (pars.at("liningConc") > 0 && pars.at("liningRadius") > 0)
                for (double xpos = pars.at("xmin"); 
                        xpos >= pars.at("xmin") - pars.at("reservoirLength"); 
                        xpos -= 4*pars.at("liningRadius") / pars.at("liningConc"))
                {
                    double ypos = 0;
                    SphericalParticle rbParticle;
                    // rbParticle.setSpecies(spec_base);
                    rbParticle.setSpecies(spec_particles);
                    rbParticle.setRadius(pars.at("liningRadius")  *
                            (1 + pars.at("liningDispersity") * random.getRandomNumber(-1,1)));
                    rbParticle.setPosition(Vec3D(xpos, 0, 0));
                    rbParticle.setVelocity(Vec3D(0,0,0));
                    rbParticle.fixParticle();
                    particleHandler.copyAndAddObject(rbParticle);
                    rbParticle.setPosition(Vec3D(xpos, pars.at("reservoirHeight"), 0));
                    particleHandler.copyAndAddObject(rbParticle);
                }

            /* Rough base for x > 0 */
            if (pars.at("baseConc") > 0 && pars.at("baseRadius"))
                for (double xpos = 0; 
                        xpos <= pars.at("xmax"); 
                        xpos += 4*pars.at("baseRadius") / pars.at("baseConc"))
                {
                    double ypos = 0;
                    SphericalParticle rbParticle;
                    // rbParticle.setSpecies(spec_base);
                    rbParticle.setSpecies(spec_particles);
                    rbParticle.setRadius(pars.at("baseRadius")  *
                            (1 + pars.at("baseDispersity") * random.getRandomNumber(-1,1)));
                    rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                    rbParticle.setVelocity(Vec3D(0,0,0));
                    rbParticle.fixParticle();

                    particleHandler.copyAndAddObject(rbParticle);
                }

            /* CubeInsertionBoundary for introducing new particles */
            BaseParticle* generandum = new SphericalParticle;
            generandum->setSpecies(spec_particles);
            generandum->setRadius(pars.at("particleRadius"));
            // double velvar = pars.at("reservoirTemperature") * sqrt(pars.at("g") * pars.at("particleRadius"));
            double velvar = 0;
            insb = new CubeInsertionBoundary();
            insb->set(
                    generandum, 0,
                    Vec3D(pars.at("xmin") - pars.at("reservoirLength"), // + 6*pars.at("particleRadius"), 
                            // pars.at("reservoirHeight") - 8*pars.at("particleRadius"),
                          0, 0),
                    Vec3D(pars.at("xmin"), // - 6*pars.at("particleRadius"), 
                          pars.at("reservoirHeight"), 0),
                    Vec3D(0, -sqrt(pars.at("g") * pars.at("reservoirHeight")), 0),
                    Vec3D(0, -sqrt(pars.at("g") * pars.at("reservoirHeight")), 0));
            insb = boundaryHandler.copyAndAddObject(insb);
            insb->checkBoundaryBeforeTimeStep(this);
            stillFillingUp = true;
            


            // dataFile.setFileType(FileType::MULTIPLE_FILES);
            // fStatFile.setFileType(FileType::MULTIPLE_FILES);
            restartFile.setFileType(FileType::ONE_FILE);

        }

        /* If restarting, we need to assign the pointers properly. 
         * This is a little messy but it must be done. */
        void actionsOnRestart() override
        {
            if (wallHandler.getNumberOfObjects() == 4)
            {
                // dam = (IntersectionOfWalls*) wallHandler.getObjectById(1);
                lid = (IntersectionOfWalls*) wallHandler.getObjectById(1);
                back = (InfiniteWall*) wallHandler.getObjectById(2);

                insb = (CubeInsertionBoundary*) boundaryHandler.getObjectById(0);
                // masb = (MaserBoundary*) boundaryHandler.getObjectById(6);
                stillFillingUp = true;
            }
            else
            {
                delb = (DeletionBoundary*) boundaryHandler.getObjectById(2);
                stillFillingUp = false;
            }

        }

        void printTime() const override
        {
            Mercury2D::printTime();

            /* Are we still filling up? If not, no need to do anything here. */
//            if (!stillFillingUp) 
//                return;

            Mdouble volfracThreshold = pars.at("volfracThreshold");
            // Volume fraction in the reservoir
            Mdouble volfrac = particleHandler.getVolume() 
                / ( pars.at("reservoirLength") * pars.at("reservoirHeight") );

            logger(INFO, "volfrac %, volfracThreshold %", 
                    volfrac, volfracThreshold);


        }

        void actionsAfterTimeStep() 
        {

            /* Are we still filling up? If not, no need to do anything here. */
            if (!stillFillingUp) 
                return;

            /* Have we filled the reservoir up enough already? */

            Mdouble volfracThreshold = pars.at("volfracThreshold");
            // Volume fraction in the reservoir
            Mdouble volfrac = particleHandler.getVolume() 
                / ( pars.at("reservoirLength") * pars.at("reservoirHeight") );

            
            if (volfrac > volfracThreshold)
            {
                /* We have introduced enough particles. Get rid of the
                 * InsertionBoundary, and put in all the other boundaries
                 * (heater, maser and deletion) */
                if (boundaryHandler.getNumberOfObjects() != 1)
                    logger(ERROR, "boundaryHandler should have 1 object, but has %",
                            boundaryHandler.getNumberOfObjects());
                else 
                    boundaryHandler.clear();

                /* HeaterBoundary */
                auto heater = new HeaterBoundary();
                heater->set2D(
                          Vec3D(pars.at("xmin") - pars.at("reservoirLength"), 
                                0.0, 0.0), 
                          Vec3D(pars.at("xmin"),
                                pars.at("reservoirHeight") * 1.2, 0.0), 
                          pars.at("reservoirTemperature") 
                        );
                heater = boundaryHandler.copyAndAddObject(heater);

                /* The deletion boundary below the chute takes away any particles
                 * that have left the chute. */
                delb = new DeletionBoundary;
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

                /* Move the back wall back to make way for the Maser's shifting! */
                back->set(Vec3D(-1, 0, 0), 
                        Vec3D(pars.at("xmin") - pars.at("reservoirLength") - 6.2*pars.at("particleRadius")*(1+pars.at("dispersity")), 0, 0)
                        );

                /* Lift the gate. We need to shift all three parts of the dam
                 * leftwards in order to follow the
                 * ConstantMassFlowMaserBoundary. */
                /*
                dam->clear();
                dam->addObject(Vec3D(+1, 0, 0), 
                        Vec3D(pars.at("xmin") - 0.5*pars.at("reservoirLength") - pars.at("particleRadius") - masb->getGapSize(), 0, 0));
                // dam->addObject(Vec3D(-1, 0, 0), Vec3D(pars.at("xmin") - 9*pars.at("particleRadius"), 0, 0));
                dam->addObject(Vec3D(-1, 0, 0), 
                    Vec3D(pars.at("xmin") - 0.5*pars.at("reservoirLength") + pars.at("particleRadius") - masb->getGapSize(), 0, 0));
                dam->addObject(Vec3D(0, +1, 0), 
                        Vec3D(pars.at("xmin") + pars.at("particleRadius"), 
                            pars.at("gateHeight"), 0));
                wallHandler.removeObject(lid->getIndex());
                // wallHandler.removeObject(back->getIndex());
                */

                /* Give an impulse to all the particles in the Maser. */
                for (BaseParticle* const p : particleHandler)
                {
                    if (! p->isFixed())
                    {
                        p->setVelocity(p->getVelocity() + Vec3D(pars.at("reservoirVel"), 0, 0));
                        // fprintf(stderr, "impulsed i = %d, p = %p\n", i, p);
                    }
                }

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
            }
        }

    private:

        std::map<std::string,double> pars;
        bool stillFillingUp;


        CubeInsertionBoundary* insb;
        IntersectionOfWalls* dam;
        IntersectionOfWalls* lid;
        InfiniteWall* back;
        ConstantMassFlowMaserBoundary* masb;
        DeletionBoundary* delb;
        LinearViscoelasticFrictionSpecies *spec_particles, *spec_base;
        BaseParticle *particlePrototype, *basePrototype;

};

int main(const int argc, char* argv[]) {
    // std::set_terminate(terminationUncaughtException);
    if (argc > 1)
    {
        CappedHeatedBlasius* problem = new CappedHeatedBlasius(argv[1]);
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

