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
#include "Particles/BaseParticle.h"
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
            setYMin(-pars.at("width")/2);
            setYMax(+pars.at("width")/2);
            setZMin(0);
            setZMax(+pars.at("height"));

            dataFile.setFileType(FileType::NO_FILE);
            // dataFile.setFileType(FileType::MULTIPLE_FILES);
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
            speciesP->setTorsionFrictionCoefficient(tan(pars.at("betators")));
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
            speciesB->setTorsionFrictionCoefficient(tan(pars.at("base_betators")));
            speciesB->setTorsionStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setTorsionDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB = speciesHandler.copyAndAddObject(speciesB);

            /* Walls */
            InfiniteWall bottomwall;
            bottomwall.setSpecies(speciesB);
            bottomwall.set(Vec3D(0,0,-1), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(bottomwall);

            /* Periodic boundaries */
            PeriodicBoundary bounds, sides;
            bounds.set(Vec3D(1,0,0), -pars.at("length")/2, +pars.at("length")/2);
            boundaryHandler.copyAndAddObject(bounds);
            sides.set(Vec3D(0,1,0), -pars.at("width")/2, +pars.at("width")/2);
            boundaryHandler.copyAndAddObject(sides);


            /* Rough base */
            int NBasalParticle = pars.at("baseConc") * pars.at("width") * pars.at("length")
                / (4 * pow(pars.at("baseRadius"), 2));
            for (int n = 0; n < NBasalParticle; n++)
            {
                double r = pars.at("baseRadius") * (1 + pars.at("baseDispersity") * generator.getRandomNumber(-1,1));
                double x = generator.getRandomNumber(-0.5, 0.5) * pars.at("length");
                double y = generator.getRandomNumber(-0.5, 0.5) * pars.at("width");
                double z = 0;

                BaseParticle rbParticle;
                rbParticle.setSpecies(speciesB);
                rbParticle.setRadius(r);
                rbParticle.setPosition(Vec3D(x, y, z));
                rbParticle.fixParticle();

                particleHandler.copyAndAddObject(rbParticle);
            }
            std::cout << NBasalParticle << " base particles are generated." << std::endl;

            setGravity(Vec3D(0, 0, -1));
        }



        muICal3D(std::string parsfile, double theta) : muICal3D(parsfile) 
        {
            pars.at("theta") = theta;
            setName(parsfile.erase(parsfile.find_last_of('.')) + "-" + std::to_string(pars.at("theta") * 180. / M_PI));
        }

        void setupInitialConditions()
        {
            /* A CubeInsertionBoundary for introducing the particles. We will
             * remove this after a few (arbitrary number of) time steps. If the
             * InsertionBoundary is doing its job properly then it will stop
             * introducing particles after a while anyway. */
            BaseParticle* p0 = new BaseParticle;
            p0->setSpecies(speciesP);
            insb = new CubeInsertionBoundary;
            insb->set( p0, 6, 
                    Vec3D(  -pars.at("length")/2  + 1*pars.at("particleRadius"), 
                            -pars.at("width")/2   + 1*pars.at("particleRadius"), 
                            0*pars.at("particleRadius")
                        ),
                    Vec3D(  +pars.at("length")/2  - 2*pars.at("particleRadius"),
                            +pars.at("width")/2   - 1*pars.at("particleRadius"), 
                             pars.at("height")    + 0*pars.at("particleRadius")
                        ),
                    Vec3D(0,0,0), Vec3D(0,0,0),
                    pars.at("particleRadius") * (1-pars.at("dispersity")),
                    pars.at("particleRadius") * (1+pars.at("dispersity"))
                    );
            insb = boundaryHandler.copyAndAddObject(insb);

            /* Dam and lid to block any initial flow */
            lid = new InfiniteWall;
            lid->setSpecies(speciesB);
            lid->set(Vec3D(0,0,1), Vec3D(0, 0, pars.at("height")));
            lid = wallHandler.copyAndAddObject(lid);

            not_yet_removed_insb = true;
        }

        void actionsAfterTimeStep()
        {
            /* We remove the CubeInsertionBoundary so that it doesn't keep
             * giving new particles. After a few time steps, it should have
             * saturated the system. */
            if (getNumberOfTimeSteps() >= pars.at("saveEvery")/2 && not_yet_removed_insb
                    && getKineticEnergy() < getTotalMass()*1e-4) 
            {
                boundaryHandler.removeObject(insb->getIndex());
                wallHandler.removeObject(lid->getIndex());
                setGravity(Vec3D(pars.at("g")*sin(pars.at("theta")), 0,
                                 -pars.at("g")*cos(pars.at("theta"))));

                not_yet_removed_insb = false;
                fprintf(stderr, "time %f, removed insb, dam and lid\n", getTime());
                setTime(0);

                size_t MAX_STRLEN = 1024;
                char fn[MAX_STRLEN];
                snprintf(fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal3D_f = fopen(fn, "w");
                setbuf(muICal3D_f, NULL);
                fprintf(muICal3D_f, "time theta depth mass xmom ke\n"); 
                fprintf(stderr, "Started writing to .muI file\n");

            }

            if (!not_yet_removed_insb 
                    && getNumberOfTimeSteps() % dataFile.getSaveCount() == 0)
                fprintf(muICal3D_f, "%g %g %g %g %g %g\n",
                        getTime(), pars.at("theta"),
                        2*getCentreOfMass().Z,
                        getTotalMass(), 
                        getTotalMomentum().X, getKineticEnergy());
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
};

int main(int argc, char ** argv) 
{
    std::set_terminate(terminationUncaughtException);

    switch (argc)
    {
        case 2:
            {
                auto problem = new muICal3D(argv[1]);
                argv[1] = argv[0];
                problem->solve(argc-1, argv+1);
                delete problem;
                break;
            }
        case 3:
            {
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

    return 0;
}
