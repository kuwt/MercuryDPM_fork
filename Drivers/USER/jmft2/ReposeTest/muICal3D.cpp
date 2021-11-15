/* muICal3D - Drop some particles onto a base z=0, sloped at an angle
 * theta(t) to the horizontal. Let theta(t) vary slowly with time, so that the
 * flow may adjust to an equilibrium. Measure I (by measuring u and h) at
 * different times, and plot I against theta. 
 * Then obtain the properties mu_1, mu_2 and I_0 by fitting an appropriate
 * curve. (This will be done externally, using Matlab.) 
 * */

#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <map>

#define MAX_STRLEN 1024

void terminationUncaughtException()
{
    std::cout << "Uncaught exception!" << std::endl;
    std::cout << "Perhaps your configuration file doesn't specify all the parameters?" << std::endl;
    exit(1);
}

class muICal3D : public Mercury3D
{
    public:
        muICal3D(std::string parsfile)
        {
            std::ifstream file(parsfile);
            std::string name;
            double var;
            while (file >> name >> var)
            {
                pars[name] = var;
            }

            setName(parsfile.erase(parsfile.find_last_of('.')));

            theta = M_PI * pars.at("theta") / 180.;
            std::cout << theta << std::endl;

            double betaslide        = M_PI * pars.at("betaslide") / 180.;
            double betaroll         = M_PI * pars.at("betaroll") / 180.;
            double betators        = M_PI * pars.at("betators") / 180.;
            double base_betaslide   = M_PI * pars.at("base_betaslide") / 180.;
            double base_betaroll    = M_PI * pars.at("base_betaroll") / 180.;
            double base_betators   = M_PI * pars.at("base_betators") / 180.;

            g = pars.at("g");

            particleRadius = pars.at("particleRadius");
            baseRadius = pars.at("baseRadius");
            double rho = pars.at("rho");

            setTimeMax(pars.at("timeMax"));
            setTimeStep(pars.at("timeStep"));
            setSaveCount(pars.at("saveCount"));

            setXMin(0);
            setXMax(pars.at("length"));
            setYMin(0);
            setYMax(pars.at("width"));
            setZMin(0);
            setZMax(+pars.at("height"));

            dataFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            /* Define species */
            speciesP = new LinearViscoelasticFrictionSpecies();
            speciesP->setDensity(rho);
            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"),
                    pars.at("restitutionCoefficient"), 
                    constants::pi*pow(particleRadius, 2) * rho
            );
            speciesP->setSlidingFrictionCoefficient(tan(betaslide)); 
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(tan(betaroll));
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(tan(betators));
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP = speciesHandler.copyAndAddObject(speciesP);

            speciesB = new LinearViscoelasticFrictionSpecies();
            speciesB->setDensity(rho);
            speciesB->setCollisionTimeAndRestitutionCoefficient(
                    pars.at("collisionTime"),
                    pars.at("restitutionCoefficient"), 
                    constants::pi*pow(baseRadius, 2) * rho
            );
            speciesB->setSlidingFrictionCoefficient(tan(base_betaslide)); 
            speciesB->setSlidingStiffness(2.0/7.0 * speciesB->getStiffness());
            speciesB->setSlidingDissipation(2.0/7.0 * speciesB->getDissipation());
            speciesB->setRollingFrictionCoefficient(tan(base_betaroll));
            speciesB->setRollingStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setRollingDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB->setTorsionFrictionCoefficient(tan(base_betators));
            speciesB->setTorsionStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setTorsionDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB = speciesHandler.copyAndAddObject(speciesB);

            /* Walls */
            InfiniteWall bottomwall;
            bottomwall.setSpecies(speciesB);
            bottomwall.set(Vec3D(0 ,0, -1), Vec3D(0, 0, 0));
            wallHandler.copyAndAddObject(bottomwall);

            /* Periodic boundaries */
            PeriodicBoundary xbounds, ybounds;
            xbounds.set(Vec3D(1, 0, 0), getXMin(), getXMax());
            ybounds.set(Vec3D(0, 1, 0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(ybounds);
            boundaryHandler.copyAndAddObject(ybounds);

            /* Rough base */
            if (pars.at("baseConc") > 0)
                for (double xpos = getXMin(); 
                    xpos <= getXMax();
                    xpos += 4*baseRadius / sqrt(pars.at("baseConc"))
                ) {
                    for (double ypos = getYMin();
                            ypos <= getYMax();
                            ypos += 4*baseRadius / sqrt(pars.at("baseConc")))
                    {
                        double zpos = 0;

                        SphericalParticle rbParticle;
                        rbParticle.setSpecies(speciesB);
                        rbParticle.setRadius(baseRadius *
                                (1 + pars.at("baseDispersity") * generator.getRandomNumber(-1,1)));
                        rbParticle.setPosition(Vec3D(xpos, ypos, zpos));
                        rbParticle.fixParticle();
                        particleHandler.copyAndAddObject(rbParticle);
                    }
                }
        }


        /*
         * Create a setup with the slope at the angle specified in radians.
         */
        muICal3D(std::string parsfile, double newTheta) : muICal3D(parsfile) 
        {
            theta = newTheta;
            const double degrees = theta * 180. / M_PI;
            setName(parsfile.erase(parsfile.find_last_of('.')) + "-" + std::to_string(degrees));
        }

        void setupInitialConditions()
        {
            /* Start writing to the .muI file. */
            char muICal3D_fn[MAX_STRLEN];
            snprintf(muICal3D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
            muICal3D_f = fopen(muICal3D_fn, "w");
            setbuf(muICal3D_f, NULL);
            fprintf(muICal3D_f, "time theta n depth mass xmom ke basfx basfy\n"); 
            fprintf(stderr, "Started writing to .muI file\n");

            setGravity(Vec3D(0, 0, -g));

            /* A CubeInsertionBoundary for introducing the particles. We will
             * remove this after a few (arbitrary number of) timesteps. If the
             * InsertionBoundary is doing its job properly then it will stop
             * introducing particles after a while anyway. */
            BaseParticle* p0 = new SphericalParticle;
            p0->setSpecies(speciesP);
            insb = new CubeInsertionBoundary;
            insb->set(p0, 1,
                    Vec3D(
                        getXMin() + particleRadius,
                        getYMin() + particleRadius,
                        getZMin()
                    ),
                    Vec3D(
                        getXMax() - 2*particleRadius,
                        getYMax() - particleRadius,
                        getZMax()
                    ),
                    Vec3D(0,0,0), Vec3D(0,0,0),
                    particleRadius * (1-pars.at("dispersity")),
                    particleRadius * (1+pars.at("dispersity"))
                    );
            insb = boundaryHandler.copyAndAddObject(insb);
            insb->insertParticles(this);

            /* Dam and lid to block any initial flow */
            lid = new InfiniteWall;
            lid->setSpecies(speciesB);
            lid->set(Vec3D(0,0,1), Vec3D(0, 0, getZMax()));
            lid = wallHandler.copyAndAddObject(lid);

            not_yet_removed_insb = true;
        }

        void actionsOnRestart()
        {
            if (boundaryHandler.getNumberOfObjects() == 1)
            {
                not_yet_removed_insb = false;

                /* Continue writing to the .muI file. */
                char muICal3D_fn[MAX_STRLEN];
                snprintf(muICal3D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal3D_f = fopen(muICal3D_fn, "a");
                setbuf(muICal3D_f, NULL);
                logger(INFO, "Continuing writing to .muI file\n");

            }
        }

        void actionsAfterTimeStep()
        {
            if (!not_yet_removed_insb)
                return;

            /* We remove the CubeInsertionBoundary so that it doesn't keep
             * giving new particles. After a few timesteps, it should have
             * saturated the system. */
            if (getNumberOfTimeSteps() >= dataFile.getSaveCount() / 2
                    && getKineticEnergy() < getTotalMass()*1e-4)
            {
                boundaryHandler.removeObject(insb->getIndex());
                wallHandler.removeObject(lid->getIndex());


                setGravity(Vec3D(g*sin(theta), 0, -g*cos(theta)));
                not_yet_removed_insb = false;
                logger(INFO, "time %, removed insb and lid", getTime());

                /* Start writing to output files. */
                setTime(0);
                forceWriteOutputFiles();
            }
        }

        Vec3D calculateBasalForce()
        {
            Vec3D basalForce;
            for (auto p : particleHandler)
            {
                if (p->isFixed())
                {
                    basalForce += p->getForce();
                }
            }
            basalForce += wallHandler.getObject(0)->getForce();

            return basalForce;
        }

        void writeOutputFiles()
        {
            Mercury3D::writeOutputFiles();
            // if (not_yet_removed_insb)
            //     return;

            // Calculate the forces on the basal particles and wall.
            Vec3D basalForce = calculateBasalForce();

            /* Write all these details to the .muI file. */
            fprintf(
                muICal3D_f, "%g %g %d %g %g %g %g %g %g\n",
                getTime(),
                theta,
                particleHandler.getNumberOfObjects(),
                2*getCentreOfMass().Z,
                getTotalMass(), 
                getTotalMomentum().X,
                getKineticEnergy(),
                basalForce.X,
                basalForce.Y
            );
        }

        void printTime()
        {
            Mercury3D::printTime();
            Vec3D basalForce = calculateBasalForce();
            printf(
                "%g %g %d %g %g %g %g %g %g\n",
                getTime(),
                theta,
                particleHandler.getNumberOfObjects(),
                2*getCentreOfMass().Z,
                getTotalMass(), 
                getTotalMomentum().X,
                getKineticEnergy(),
                basalForce.X,
                basalForce.Y
            );

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
        FILE * muICal3D_f;

        double g;
        double theta;  // slope angle in radians
        double particleRadius;
        double baseRadius;
};


int main(int argc, char ** argv) 
{
    std::set_terminate(terminationUncaughtException);

    auto problem = new muICal3D(argv[1], M_PI * atof(argv[2]) / 180.);
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
                auto problem = new muICal3D(argv[1]);
                argv[1] = argv[0];
                problem->solve(argc-1, argv+1);
                delete problem;
                break;
            }
        case 3:
            {
                /* An angle has been specified */
                auto problem = new muICal3D(argv[1], M_PI * atof(argv[2]) / 180.);
                argv[2] = argv[0];
                problem->solve(argc-2, argv+2);
                delete problem;
                break;
            }
        default:
            {
                fprintf(stderr, "Usage: %s pars-file [angle]\n", argv[0]);
                exit(-1);
                break;
            }
    }
#endif
