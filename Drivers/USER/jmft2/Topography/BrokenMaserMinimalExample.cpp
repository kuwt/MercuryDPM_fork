/* BrokenMaserMinimalExampleBlasius - A simplified form of JMFT's Blasius
 * drivers that doesn't work in parallel.
 *
 * Example usage: mpirun -n 2 ./BrokenMaserMinimalExampleBlasius db-minimal.pars 
 */
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/CubeInsertionBoundary.h"
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

            dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);

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
            auto base = new InfiniteWall();
            base->setSpecies(spec_particles);
            // base->setSpecies(spec_base);
            base->set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            base = wallHandler.copyAndAddObject(base);

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
                    Vec3D(1, - sqrt(pars.at("reservoirHeight")), 0),
                    Vec3D(1, - sqrt(pars.at("reservoirHeight")), 0),
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

        void actionsAfterTimeStep() 
        {
            // logger(INFO, "In actionsAfterTimeStep()");

            // logger(INFO, "t = %, np = %", getTime(), particleHandler.getNumberOfRealObjectsLocal());
            /* Are we still filling up? If not, no need to do anything here. */
            if (!stillFillingUp) 
                return;

            /* Have we filled the reservoir up enough already? */

            if (getTime() > 20*getTimeStep())
            {
                boundaryHandler.removeObject(insb->getId());

                /* MaserBoundary */
                auto masb = new SubcriticalMaserBoundaryTEST();
                masb->set(Vec3D(1.0, 0.0, 0.0), 
                        pars.at("xmin") - pars.at("reservoirLength") + pars.at("particleRadius"), 
                        pars.at("xmin"));
                logger(INFO, "About to put in masb");
                masb = boundaryHandler.copyAndAddObject(masb);
                logger(INFO, "Have put in the masb");
                logger(INFO, "About to activate masb");
                masb->activateMaser();
                logger(INFO, "Have activated masb");


                stillFillingUp = false;

            }
            // logger(INFO, "Have completed actionsAfterTimeStep");
        }

    private:

        std::map<std::string,double> pars;
        bool stillFillingUp;

        CubeInsertionBoundary* insb;
        LinearViscoelasticFrictionSpecies *spec_particles;

};

int main(const int argc, char* argv[]) {
    if (argc > 1)
    {
        BrokenMaserMinimalExample* problem = new BrokenMaserMinimalExample(argv[1]);
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
