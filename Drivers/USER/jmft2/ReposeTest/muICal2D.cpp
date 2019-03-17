/* muICal2D - Drop some particles onto a base y=0, sloped at an angle
 * theta(t) to the horizontal. Let theta(t) vary slowly with time, so that the
 * flow may adjust to an equilibrium. Measure I (by measuring u and h) at
 * different times, and plot I against theta. Then obtain the properties mu_1,
 * mu_2 and I_0 by fitting an appropriate curve. 
 *
 * This is similar to the periodic current test (CurrentTest2D), but this one is
 * meant to tell you about steady flows, not starting and stopping.
 * */

#include "Mercury2D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#include <map>

#include "CG/CG.h"
#include "CG/TimeSmoothedCG.h"
#include "CG/TimeAveragedCG.h"

#define MAX_STRLEN 1024

void terminationUncaughtException()
{
    std::cout << "Uncaught exception!" << std::endl;
    std::cout << "Perhaps your configuration file doesn't specify all the parameters?" << std::endl;
    exit(1);
}

class muICal2D : public Mercury2D
{
    public:
        muICal2D(std::string parsfile)
        {

            std::ifstream file(parsfile);
            std::string name;
            double var;
            while (file >> name >> var)
            {
                pars[name] = var;
            }

            setName(parsfile.erase(parsfile.find_last_of('.')));

            std::cout << pars.at("theta") << std::endl;

            pars.at("theta")            = M_PI * pars.at("theta") / 180.;
            pars.at("betaslide")        = M_PI * pars.at("betaslide") / 180.;
            pars.at("betaroll")         = M_PI * pars.at("betaroll") / 180.;
            pars.at("base_betaslide")   = M_PI * pars.at("base_betaslide") / 180.;
            pars.at("base_betaslide")   = M_PI * pars.at("base_betaslide") / 180.;

            setTimeMax(pars.at("timeMax"));
            setTimeStep(pars.at("timeStep"));
            setSaveCount(pars.at("saveEvery"));

            setXMin(-pars.at("length")/2);
            setXMax(+pars.at("length")/2);
            setYMin(0);
            setYMax(pars.at("height"));
            setZMin(0);
            setZMax(1);

            dataFile.setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::NO_FILE);

            /* Define species */
            speciesP = new LinearViscoelasticFrictionSpecies();
            speciesP->setDensity(pars.at("rho"));
            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"),
                    pars.at("restitutionCoefficient"), 
                    constants::pi*pow(pars.at("particleRadius"), 2)*pars.at("rho")
            );
            speciesP->setSlidingFrictionCoefficient(tan(pars.at("betaslide"))); 
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(tan(pars.at("betaroll")));
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(0);
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP = speciesHandler.copyAndAddObject(speciesP);

            speciesB = new LinearViscoelasticFrictionSpecies();
            speciesB->setDensity(pars.at("rho"));
            speciesB->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"),
                    pars.at("restitutionCoefficient"), 
                    constants::pi*pow(pars.at("baseRadius"), 2)*pars.at("rho")
            );
            speciesB->setSlidingFrictionCoefficient(tan(pars.at("base_betaslide"))); 
            speciesB->setSlidingStiffness(2.0/7.0 * speciesB->getStiffness());
            speciesB->setSlidingDissipation(2.0/7.0 * speciesB->getDissipation());
            speciesB->setRollingFrictionCoefficient(tan(pars.at("base_betaroll")));
            speciesB->setRollingStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setRollingDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB->setTorsionFrictionCoefficient(0);
            speciesB->setTorsionStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setTorsionDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB = speciesHandler.copyAndAddObject(speciesB);

            /* Walls */
            InfiniteWall bottomwall;
            bottomwall.setSpecies(speciesB);
            bottomwall.set(Vec3D(0,-1,0), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(bottomwall);

            /* Periodic boundaries */
            PeriodicBoundary bounds;
            bounds.set(Vec3D(1,0,0), -pars.at("length")/2, +pars.at("length")/2);
            boundaryHandler.copyAndAddObject(bounds);

            /* Rough base */
            if (pars.at("baseConc") > 0)
                for (double xpos = -pars.at("length")/2; 
                        xpos <= pars.at("length")/2;
                        xpos += 4*pars.at("baseRadius") / pars.at("baseConc"))
                {
                    double ypos = 0;

                    BaseParticle rbParticle;
                    rbParticle.setSpecies(speciesB);
                    rbParticle.setRadius(pars.at("baseRadius") * 
                            (1 + pars.at("baseDispersity") * generator.getRandomNumber(-1,1)));
                    rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                    rbParticle.setVelocity(Vec3D(0,0,0));
                    rbParticle.fixParticle();
                    particleHandler.copyAndAddObject(rbParticle);
                }

            /* Gravity initially points downwards */
            setGravity(Vec3D(0, -1, 0));
        }

        muICal2D(std::string parsfile, double theta) : muICal2D(parsfile) 
        {
            pars["theta"] = theta;

            char* name = (char*)calloc(sizeof(char), 1024);
            snprintf(name, 1024, "%s-%.1f",
                    parsfile.erase(parsfile.find_last_of('.')).c_str(), 
                    pars.at("theta") * 180. / M_PI);
            setName(name);
            free(name);

            // setName(parsfile.erase(parsfile.find_last_of('.')) + "-" + std::to_string(pars.at("theta") * 180. / M_PI));
        }

        void setupInitialConditions()
        {
            /* A CubeInsertionBoundary for introducing the particles. We will
             * remove this after a few (arbitrary number of) timesteps. If the
             * InsertionBoundary is doing its job properly then it will stop
             * introducing particles after a while anyway. */
            BaseParticle* p0 = new BaseParticle;
            p0->setSpecies(speciesP);
            insb = new CubeInsertionBoundary;
            insb->set( p0, 1,
                    Vec3D( -pars.at("length")/2 + 1*pars.at("particleRadius"),
                            0, 0),
                    Vec3D( +pars.at("length")/2 - 2*pars.at("particleRadius"),
                            pars.at("height"), 0),
                    /*
                    -1*sqrt(pars.at("g") * pars.at("particleRadius")) * Vec3D(1,1,0),
                    1*sqrt(pars.at("g") * pars.at("particleRadius")) * Vec3D(1,1,0),
                    */
                    Vec3D(0,0,0), Vec3D(0,0,0),
                    pars.at("particleRadius") * (1-pars.at("dispersity")),
                    pars.at("particleRadius") * (1+pars.at("dispersity"))
                    );
            insb = boundaryHandler.copyAndAddObject(insb);

            /* A dam for blocking any initial flow. */
            /*
            dam = new InfiniteWall;
            dam->setSpecies(speciesB);
            dam->set( Vec3D(1,0,0), Vec3D(pars.at("length")/2 - 0.5*pars.at("particleRadius"), 0, 0) );
            dam = wallHandler.copyAndAddObject(dam);
            */
            lid = new InfiniteWall;
            lid->setSpecies(speciesB);
            lid->set( Vec3D(0,1,0), Vec3D(0, pars.at("height"), 0) );
            lid = wallHandler.copyAndAddObject(lid);

            not_yet_removed_insb = true;

            /* Create a CG object for live statistics */
            auto cg = cgHandler.copyAndAddObject(
                    new CG<CGCoordinates::Y>);
            cg->statFile.setSaveCount(pars.at("saveEvery"));
            cg->setWidth(pars.at("particleRadius") * 1.5);
            cg->setY(0, pars.at("height")*1.5);
            cg->setN(64);
        }

        void actionsOnRestart()
        {
            if (boundaryHandler.getNumberOfObjects() == 1)
            {
                not_yet_removed_insb = false;

                /* Continue writing to the .muI file. */
                char muICal2D_fn[MAX_STRLEN];
                snprintf(muICal2D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal2D_f = fopen(muICal2D_fn, "a");
                setbuf(muICal2D_f, NULL);
                logger(INFO, "Continuing writing to .muI file\n");

            }
        }

        void actionsAfterTimeStep()
        {
            /* We remove the CubeInsertionBoundary so that it doesn't keep
             * giving new particles. After a few timesteps, it should have
             * saturated the system. */
            if (getNumberOfTimeSteps() >= pars.at("saveEvery")/2 && not_yet_removed_insb
                    && getKineticEnergy() < getTotalMass()*1e-4) 
            {
                boundaryHandler.removeObject(insb->getIndex());
                // wallHandler.removeObject(dam->getIndex());
                wallHandler.removeObject(lid->getIndex());


                /* Gravity: Only start the sloping after the container has been
                 * filled.*/
                setGravity(Vec3D(pars.at("g")*sin(pars.at("theta")), -pars.at("g")*cos(pars.at("theta")), 0.0));
                not_yet_removed_insb = false;
                logger(INFO, "time %, removed insb, dam and lid", getTime());

                /* Start writing to output files. */
                setTime(0);
                // dataFile.setFileType(FileType::MULTIPLE_FILES);
                // fStatFile.setFileType(FileType::MULTIPLE_FILES);
                forceWriteOutputFiles();

                /* Start writing to the .muI file. */
                char muICal2D_fn[MAX_STRLEN];
                snprintf(muICal2D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal2D_f = fopen(muICal2D_fn, "w");
                setbuf(muICal2D_f, NULL);
                fprintf(muICal2D_f, "time theta depth mass xmom ke basfx basfy\n"); 
                logger(INFO, "Started writing to .muI file\n");
            }

            /* If the experiment has started, then write to the .muI file. */
            if (!not_yet_removed_insb 
                    && getNumberOfTimeSteps() % dataFile.getSaveCount() == 0)
            {
                /* Calculate the forces on the basal particles 
                 * and the basal wall. */
                Vec3D basalForce;
                for (auto p : particleHandler)
                {
                    if (p->isFixed())
                    {
                        basalForce += p->getForce();
                    }
                }
                basalForce += wallHandler.getObject(0)->getForce();
                
                /* Write all these details to the .muI file. */
                fprintf(muICal2D_f, "%g %g %g %g %g %g %g %g\n",
                        getTime(), pars.at("theta"),
                        2*getCentreOfMass().Y,
                        getTotalMass(), 
                        getTotalMomentum().X, getKineticEnergy(),
                        basalForce.X, basalForce.Y);
            }
        }

        void actionsAfterSolve()
        {
            dataFile.setFileType(FileType::MULTIPLE_FILES);
            writeDataFile();
            fStatFile.setFileType(FileType::MULTIPLE_FILES);
            writeFStatFile();
        }

    private:
        std::map<std::string, double> pars;
        RNG generator;
        LinearViscoelasticFrictionSpecies* speciesP;
        LinearViscoelasticFrictionSpecies* speciesB;
        CubeInsertionBoundary* insb;
        InfiniteWall*          dam;
        InfiniteWall*          lid;
        bool not_yet_removed_insb;
        FILE * muICal2D_f;
};

int main(int argc, char ** argv) 
{
    std::set_terminate(terminationUncaughtException);

    auto problem = new muICal2D(argv[1], M_PI * atof(argv[2]) / 180.);
    argv[2] = argv[0];
    problem->solve(argc-2, argv+2);
    delete problem;
    return 0;
}

#if 0
    switch (argc)
    {
        case 2:
            {
                /* No angle has been specified, 
                 * so choose from the pars file. */
                auto problem = new muICal2D(argv[1]);
                argv[1] = argv[0];
                problem->solve(argc-1, argv+1);
                delete problem;
                break;
            }
        case 3:
            {
                /* An angle has been specified */
                auto problem 
                    = new muICal2D(argv[1], M_PI * atof(argv[2]) / 180.);
                argv[2] = argv[0];
                problem->solve(argc-2, argv+2);
                delete problem;
                break;
            }
        default:
            {
                fprintf(stderr, 
                    "Usage: %s pars-file [angle]\n", argv[0]);
                exit(-1);
                break;
            }
    }
#endif
