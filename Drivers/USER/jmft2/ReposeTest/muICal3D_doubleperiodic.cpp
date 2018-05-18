/* muICal3D_doubleperiodic - Drop some particles onto a base z=0, sloped at an angle
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
#define MAX_STRLEN 1024

#include "muICal3D.h"

class muICal3D_doubleperiodic : public Mercury3D
{
    public:
        muICal3D_doubleperiodic(muICal3D_pars_t pars_in)
        {
            fprintf(stderr, "constructor for muICal3D_doubleperiodic\n");
            pars = pars_in;

            setTimeMax(pars.timeMax);
            setTimeStep(pars.timeStep);
            setSaveCount(pars.saveEvery);

            setXMin(-pars.length/2);
            setXMax(+pars.length/2);
            setYMin(-pars.width/2);
            setYMax(+pars.width/2);
            setZMin(0);
            setZMax(+pars.height);

            dataFile.setFileType(FileType::NO_FILE);
            // dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::NO_FILE);

            /* Define species */
            speciesP = new LinearViscoelasticFrictionSpecies();
            speciesP->setDensity(pars.rho);
            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    pars.collisionTime,
                    pars.restitutionCoefficient, 
                    constants::pi*pow(pars.particleRadius, 2)*pars.rho
            );
            speciesP->setSlidingFrictionCoefficient(tan(pars.betaslide)); 
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(tan(pars.betaroll));
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(tan(pars.betators));
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP = speciesHandler.copyAndAddObject(speciesP);

            speciesB = new LinearViscoelasticFrictionSpecies();
            speciesB->setDensity(pars.rho);
            speciesB->setCollisionTimeAndRestitutionCoefficient(
                    pars.collisionTime,
                    pars.restitutionCoefficient, 
                    constants::pi*pow(pars.baseRadius, 2)*pars.rho
            );
            speciesB->setSlidingFrictionCoefficient(tan(pars.base_betaslide)); 
            speciesB->setSlidingStiffness(2.0/7.0 * speciesB->getStiffness());
            speciesB->setSlidingDissipation(2.0/7.0 * speciesB->getDissipation());
            speciesB->setRollingFrictionCoefficient(tan(pars.base_betaroll));
            speciesB->setRollingStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setRollingDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB->setTorsionFrictionCoefficient(tan(pars.base_betators));
            speciesB->setTorsionStiffness(2.0/5.0 * speciesB->getStiffness());
            speciesB->setTorsionDissipation(2.0/5.0 * speciesB->getDissipation());
            speciesB = speciesHandler.copyAndAddObject(speciesB);

            /* Walls */
            InfiniteWall bottomwall, leftwall, rightwall;
            bottomwall.setSpecies(speciesB);
            bottomwall.set(Vec3D(0,0,-1), Vec3D(0,0,0));
            wallHandler.copyAndAddObject(bottomwall);

            /* Periodic boundaries */
            PeriodicBoundary bounds;
            bounds.set(Vec3D(1,0,0), -pars.length/2, +pars.length/2);
            boundaryHandler.copyAndAddObject(bounds);

            PeriodicBoundary sides;
            sides.set(Vec3D(0,1,0), -pars.width/2, +pars.width/2);
            boundaryHandler.copyAndAddObject(sides);

            /* Rough base */
            BaseParticle basalp;
            basalp.setSpecies(speciesB);
            int numberBase = 
                pars.baseConc * (pars.length * pars.width) / (M_PI * pow(pars.baseRadius,2));
            for (int i = 0; i < numberBase; i++)
            {
                basalp.setRadius( pars.particleRadius 
                                    * generator.getRandomNumber(
                                        1-pars.baseDispersity,
                                        1+pars.baseDispersity) 
                    );
                basalp.setPosition(Vec3D(
                            pars.length * generator.getRandomNumber(-0.5, 0.5),
                            pars.width  * generator.getRandomNumber(-0.5, 0.5),
                            0.0)
                    );
                basalp.fixParticle();
                particleHandler.copyAndAddObject(basalp);
            }

            setGravity(Vec3D(0,0,-1));
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
            insb->set( p0, 400, 
                    Vec3D(  -pars.length/2  + 1*pars.particleRadius, 
                            -pars.width/2   + 1*pars.particleRadius, 
                            0*pars.particleRadius
                        ),
                    Vec3D(  +pars.length/2  - 2*pars.particleRadius,
                            +pars.width/2   - 1*pars.particleRadius, 
                             pars.height    + 0*pars.particleRadius
                        ),
                    Vec3D(0,0,0), Vec3D(0,0,0),
                    pars.particleRadius * (1-pars.dispersity),
                    pars.particleRadius * (1+pars.dispersity)
                    );
            insb = boundaryHandler.copyAndAddObject(insb);

            /* Dam and lid to block any initial flow */
            lid = new InfiniteWall;
            lid->setSpecies(speciesB);
            lid->set(Vec3D(0,0,1), Vec3D(0, 0, pars.height));
            lid = wallHandler.copyAndAddObject(lid);

            not_yet_removed_insb = true;
        }

        void actionsAfterTimeStep()
        {
            /* We remove the CubeInsertionBoundary so that it doesn't keep
             * giving new particles. After a few time steps, it should have
             * saturated the system. */
            if (getNumberOfTimeSteps() >= pars.saveEvery/2 && not_yet_removed_insb
                    && getKineticEnergy() < getTotalMass()*1e-4) 
            {
                boundaryHandler.removeObject(insb->getIndex());
                wallHandler.removeObject(lid->getIndex());
                setGravity(Vec3D(pars.g*sin(pars.theta), 0,
                                 -pars.g*cos(pars.theta)));

                not_yet_removed_insb = false;
                fprintf(stderr, "time %f, removed insb, dam and lid\n",
                        getTime());
                
                setTime(0);

                char fn[MAX_STRLEN];
                snprintf(fn, MAX_STRLEN, "%s.muI", getName().c_str());
                muICal3D_doubleperiodic_f = fopen(fn, "w");
                setbuf(muICal3D_doubleperiodic_f, NULL);
                fprintf(muICal3D_doubleperiodic_f, "time theta depth mass xmom ke\n"); 
                fprintf(stderr, "Started writing to .muI file\n");

            }

            if (!not_yet_removed_insb 
                    && getNumberOfTimeSteps() % pars.saveEvery == 0)
                fprintf(muICal3D_doubleperiodic_f, "%g %g %g %g %g %g\n",
                        getTime(), pars.theta,
                        2*getCentreOfMass().Z,
                        getTotalMass(), 
                        getTotalMomentum().X, getKineticEnergy());
        }

    private:
        muICal3D_pars_t pars;
        RNG generator;
        LinearViscoelasticFrictionSpecies* speciesP;
        LinearViscoelasticFrictionSpecies* speciesB;
        CubeInsertionBoundary* insb;
        InfiniteWall*          dam;
        InfiniteWall*          lid;
        bool not_yet_removed_insb;
        FILE * muICal3D_doubleperiodic_f;
};

int main(int argc, char ** argv) 
{
    if ((argc == 2) || (argc == 3))
    {
        muICal3D_pars_t pars;
        char parsfn[MAX_STRLEN];
        char* probname;
        if (argc == 2)
        {
            probname = argv[1];
            snprintf(parsfn, MAX_STRLEN, "%s.pars", probname);
            muICal3D_pars_read(&pars, parsfn);
        }
        else if (argc == 3)
        {
            char* speciesname = argv[1];
            snprintf(parsfn, MAX_STRLEN, "%s.pars", speciesname);
            muICal3D_pars_read(&pars, parsfn);
            double angle = atof(argv[2]);
            pars.theta = angle * M_PI / 180.;
            probname = (char*)calloc(MAX_STRLEN, sizeof(char));
            snprintf(probname, MAX_STRLEN, "%s-%.1f", speciesname, angle);
        }
        else 
            assert(false);

        fprintf(stderr, "About to run with the following parameters\n");
        muICal3D_pars_write(stderr, pars);
        muICal3D_doubleperiodic problem(pars);
        problem.setName(probname);
        problem.solve();
    }
    else
    {
        muICal3D_pars_t pars;
        muICal3D_pars_preset(&pars);
        muICal3D_pars_write(stdout, pars);
    }
    return 0;
}
