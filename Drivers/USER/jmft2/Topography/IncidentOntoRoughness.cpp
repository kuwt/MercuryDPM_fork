/* IncidentOntoRoughness - Flow on a slope that starts smooth and suddenly becomes rough */
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/FluxBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"

#include "IncidentOntoRoughness.h"

class IncidentOntoRoughness : public Mercury2D {
    public:
        IncidentOntoRoughness(IncidentOntoRoughness_pars_t pars_in)
        {
            pars = pars_in;
            setTimeStep(pars.timeStep);
            setTimeMax(pars.timeMax);
            setSaveCount(pars.saveEvery);

            setXMin(pars.xmin - pars.reservoirLength);
            setXMax(pars.xmax + pars.reservoirLength);
            setYMin(0);
            setYMax(pars.reservoirHeight);
            setZMin(0);
            setZMax(1);

            // Working in the tilted frame
            setGravity(Vec3D(pars.g*sin(pars.theta), -pars.g*cos(pars.theta), 0.0));

            /* Species */
            species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
            species->setDensity(pars.rho);
            // Mass is area of circle times density (note that in two
            // dimensions, density is ML^-2, not ML^-3. 
            species->setCollisionTimeAndRestitutionCoefficient(
                    pars.collisionTime,
                    pars.restitutionCoefficient, 
                    M_PI * pars.rho
                        * pow( pars.particleRadius, 2)
            );
            species->setSlidingFrictionCoefficient(tan(pars.betaslide)); 
            species->setSlidingStiffness(2.0/7.0 * species->getStiffness());
            species->setSlidingDissipation(2.0/7.0 * species->getDissipation());
            species->setRollingFrictionCoefficient(tan(pars.betaroll));
            species->setRollingStiffness(2.0/5.0 * species->getStiffness());
            species->setRollingDissipation(2.0/5.0 * species->getDissipation());
            species->setTorsionFrictionCoefficient(0);
            species->setTorsionStiffness(2.0/5.0 * species->getStiffness());
            species->setTorsionDissipation(2.0/5.0 * species->getDissipation());

            /* Walls */
            InfiniteWall backWall;
            backWall.setSpecies(species);
            backWall.set(Vec3D(-1,0,0), Vec3D(pars.xmin - pars.reservoirLength, 0, 0));
            wallHandler.copyAndAddObject(backWall);

            IntersectionOfWalls base;
            base.setSpecies(species);
            base.addObject(Vec3D(0, -1, 0), Vec3D(0,0,0));
            base.addObject(Vec3D(-1, 0, 0), Vec3D(pars.xmax, 0,0));
            wallHandler.copyAndAddObject(base);

            /* Rough base */
            int number_bp = floor( 
                 pars.baseConc * (pars.xmax) / (2 * pars.baseRadius) );
            for (int i = 0; i < number_bp; i++) 
            {
                double xpos = generator.getRandomNumber(0, pars.xmax);
                double ypos = 0;

                BaseParticle rbParticle;
                rbParticle.setSpecies(species);
                rbParticle.setRadius(pars.baseRadius  * 
                        (1 + pars.baseDispersity * generator.getRandomNumber(-1,1)));
                rbParticle.setPosition(Vec3D(xpos, ypos, 0));
                rbParticle.setVelocity(Vec3D(0,0,0));
                rbParticle.fixParticle();
                particleHandler.copyAndAddObject(rbParticle);
            }

            /* The deletion boundary below the chute takes away any particles
             * that have left the chute. */
            delb = new DeletionBoundary;
            // delb.set(Vec3D(0,-1,0), pars.reservoirHeight);
            delb->set(Vec3D(1,0,0), pars.xmax);
            delb = boundaryHandler.copyAndAddObject(delb);

            /* CubeInsertionBoundary for introducing new particles */
            BaseParticle* generandum = new BaseParticle;
            generandum->setSpecies(species);
            generandum->setRadius(pars.particleRadius);
            double velvar = pars.reservoirTemperature * sqrt(pars.g * pars.particleRadius);
            insb = new CubeInsertionBoundary();
            insb->set(
                generandum, 100, 
                    Vec3D(pars.xmin - pars.reservoirLength + pars.particleRadius, 
                          0, 0),
                    Vec3D(pars.xmin - 2*pars.particleRadius, 
                          pars.reservoirHeight, 0),
                    Vec3D(pars.reservoirVel - velvar, -pars.reservoirVel - velvar, 0),
                    Vec3D(pars.reservoirVel + velvar, -pars.reservoirVel + velvar, 0),
                    pars.particleRadius * (1 - pars.dispersity),
                    pars.particleRadius * (1 + pars.dispersity)
                );
            insb = boundaryHandler.copyAndAddObject(insb);

            /* Flux boundaries for measuring flow rates. */
            fluxb_upstream = new FluxBoundary;
            fluxb_upstream->set(Vec3D(1,0,0), pars.xmin);
            fluxb_upstream = boundaryHandler.copyAndAddObject(fluxb_upstream);
            fluxb_middle = new FluxBoundary;
            fluxb_middle->set(Vec3D(1,0,0), 0);
            fluxb_middle = boundaryHandler.copyAndAddObject(fluxb_middle);
            fluxb_downstream = new FluxBoundary;
            fluxb_downstream->set(Vec3D(1,0,0), pars.xmax / 2);
            fluxb_downstream = boundaryHandler.copyAndAddObject(fluxb_downstream);

            dataFile.setFileType(FileType::MULTIPLE_FILES);
            fStatFile.setFileType(FileType::MULTIPLE_FILES);
            restartFile.setFileType(FileType::ONE_FILE);
        }

        ~IncidentOntoRoughness(void) {
            fclose(flux_f);
        }

        void setupInitialConditions() 
        {
            char flux_fn[1024];
            snprintf(flux_fn, 1024, "%s.flux", getName().c_str());
            flux_f = fopen(flux_fn, "w");
            setlinebuf(flux_f);
            fprintf(flux_f, 
                    "time mass insb fluxb_upstream fluxb_middle fluxb_downstream delb\n");
        }

        void actionsAfterTimeStep() 
        {
            if (getNtimeSteps() % pars.saveEvery == 0)
            {
                fprintf(flux_f, "%f %f %f %f %f %f %f\n",
                        getTime(), 
                        getTotalMass(),
                        insb->getMassOfParticlesInserted(),
                        fluxb_upstream->getMassOfParticlesCrossedForw() - fluxb_upstream->getMassOfParticlesCrossedBack(),
                        fluxb_middle->getMassOfParticlesCrossedForw() - fluxb_middle->getMassOfParticlesCrossedBack(),
                        fluxb_downstream->getMassOfParticlesCrossedForw() - fluxb_downstream->getMassOfParticlesCrossedBack(),
                        delb->getMassOfParticlesDeleted());
            }
        }


    private:
        IncidentOntoRoughness_pars_t pars;
        CubeInsertionBoundary* insb;
        FluxBoundary* fluxb_upstream;
        FluxBoundary* fluxb_middle;
        FluxBoundary* fluxb_downstream;
        DeletionBoundary* delb;
        LinearViscoelasticFrictionSpecies* species;
        RNG generator;

        FILE* flux_f;

};

int main(const int argc, char* argv[]) {
    IncidentOntoRoughness_pars_t pars = IncidentOntoRoughness_pars_default();

    if (argc == 1 || strcmp(argv[1], "--default-pars") == 0)
    {
        IncidentOntoRoughness_pars_write(stdout, IncidentOntoRoughness_pars_default());
        exit(0);
    }

    char* probname = argv[1];
    char parsfn[1024]; 
    snprintf(parsfn, 1024, "%s.pars", probname);
    IncidentOntoRoughness_pars_read(&pars, parsfn);
    IncidentOntoRoughness problem(pars);
    problem.setName(probname);

    argv[1] = argv[0];
    problem.solve(argc-1,argv+1);
    return 0;
}
