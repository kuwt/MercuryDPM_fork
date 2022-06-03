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
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#define MAX_STRLEN 1024
#define DEGREES (M_PI / 180.)


class muICal2D : public Mercury2D
{
    public:
        muICal2D(std::string parsfile, double thetaInDegrees)
        {

            std::ifstream file(parsfile);
            std::string name;
            double var;
            while (file >> name >> var)
            {
                pars[name] = var;
            }

            theta = thetaInDegrees * DEGREES;
            setName(parsfile.erase(parsfile.find_last_of('.')) + "-" + std::to_string(thetaInDegrees));

            std::cout << theta << std::endl;

            pars.at("betaslide")        = pars.at("betaslide") * DEGREES;
            pars.at("betaroll")         = pars.at("betaroll") * DEGREES;
            pars.at("base_betaslide")   = pars.at("base_betaslide") * DEGREES;
            pars.at("base_betaslide")   = pars.at("base_betaslide") * DEGREES;

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
            bottomwall.set(Vec3D(0, -1, 0), Vec3D(0, 0, 0));
            wallHandler.copyAndAddObject(bottomwall);

            /* Periodic boundaries */
            PeriodicBoundary bounds;
            bounds.set(Vec3D(1, 0, 0), getXMin(), getXMax());
            boundaryHandler.copyAndAddObject(bounds);

            /* Rough base */
            if (pars.at("baseConc") > 0)
                for (double xpos = getXMin();
                     xpos <= getXMax();
                     xpos += 4*pars.at("baseRadius") / pars.at("baseConc"))
                {
                    double ypos = 0;

                    SphericalParticle rbParticle;
                    rbParticle.setSpecies(speciesB);
                    rbParticle.setRadius(pars.at("baseRadius") *
                            (1 + pars.at("baseDispersity") * generator.getRandomNumber(-1,1)));
                    rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                    rbParticle.setVelocity(Vec3D(0, 0, 0));
                    rbParticle.fixParticle();
                    particleHandler.copyAndAddObject(rbParticle);
                }
        }

        void setupInitialConditions() override
        {
            /* Gravity initially points downwards */
            setGravity(Vec3D(0, -1, 0));

            /* A CubeInsertionBoundary for introducing the particles. We will
             * remove this after a few (arbitrary number of) timesteps. If the
             * InsertionBoundary is doing its job properly then it will stop
             * introducing particles after a while anyway. */
            BaseParticle* p0 = new SphericalParticle;
            p0->setSpecies(speciesP);
            insb = new CubeInsertionBoundary;
            PSD myPSD;
            myPSD.setDistributionUniform( pars.at("particleRadius") * (1-pars.at("dispersity")),  pars.at("particleRadius") * (1+pars.at("dispersity")),1000);
            insb->set(p0, 1,
                    Vec3D(getXMin() + 1*pars.at("particleRadius"),
                          getYMin(), 0),
                    Vec3D(getXMax() - 2*pars.at("particleRadius"),
                          getYMax(), 0),
                    Vec3D(0,0,0), Vec3D(0,0,0)
                    );
            insb->setPSD(myPSD);
            insb = boundaryHandler.copyAndAddObject(insb);

            lid = new InfiniteWall;
            lid->setSpecies(speciesB);
            lid->set(Vec3D(0,1,0), Vec3D(0, pars.at("height"), 0));
            lid = wallHandler.copyAndAddObject(lid);

            notYetRemovedInsb = true;
        }

        void actionsOnRestart() override
        {
            if (boundaryHandler.getNumberOfObjects() == 1)
            {
                notYetRemovedInsb = false;

                /* Continue writing to the .muI file. */
                char muICal2D_fn[MAX_STRLEN];
                snprintf(muICal2D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal2D_f = fopen(muICal2D_fn, "a");
                setbuf(muICal2D_f, nullptr);
                logger(INFO, "Continuing writing to .muI file\n");
            }
        }

    Vec3D calculateBasalForce()
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

        return basalForce;
    }

    void actionsAfterTimeStep() override
    {
        if (!notYetRemovedInsb)
            return;

        /* We remove the CubeInsertionBoundary so that it doesn't keep
         * giving new particles. After a few timesteps, it should have
         * saturated the system. */
        if (getNumberOfTimeSteps() >= dataFile.getSaveCount() / 2
            && getKineticEnergy() < getTotalMass()*1e-4)
            {
                boundaryHandler.removeObject(insb->getIndex());
                wallHandler.removeObject(lid->getIndex());

                /* Gravity: Only start the sloping after the container has been
                 * filled.*/
                setGravity(Vec3D(pars.at("g")*sin(theta), -pars.at("g")*cos(theta), 0.0));
                notYetRemovedInsb = false;
                logger(INFO, "time %, removed insb and lid", getTime());

                /* Start writing to output files. */
                setTime(0);
                forceWriteOutputFiles();

                /* Start writing to the .muI file. */
                char muICal2D_fn[MAX_STRLEN];
                snprintf(muICal2D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal2D_f = fopen(muICal2D_fn, "w");
                setbuf(muICal2D_f, nullptr);
                fprintf(muICal2D_f, "time theta depth mass xmom ke basfx basfy\n");
                logger(INFO, "Started writing to .muI file\n");
            }

            /* If the experiment has started, then write to the .muI file. */
            if (notYetRemovedInsb)
                return;

            if (getNumberOfTimeSteps() % dataFile.getSaveCount() != 0)
                return;

            // Calculate the forces on the basal particles and wall.
            Vec3D basalForce = calculateBasalForce();

            // Write all these details to the .muI file.
            fprintf(
                muICal2D_f, "%g %g %g %g %g %g %g %g\n",
                getTime(),
                theta / DEGREES,
                2*getCentreOfMass().Y,
                getTotalMass(),
                getTotalMomentum().X,
                getKineticEnergy(),
                basalForce.X,
                basalForce.Y
            );
        }

        void actionsAfterSolve() override
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
        InfiniteWall*          lid;
        bool notYetRemovedInsb;
        FILE * muICal2D_f;

        double theta;  // slope angle in radians
};


int main(int argc, char ** argv)
{
    auto problem = new muICal2D(argv[1], atof(argv[2]));
    argv[2] = argv[0];
    problem->solve(argc-2, argv+2);
    delete problem;
    return 0;
}
