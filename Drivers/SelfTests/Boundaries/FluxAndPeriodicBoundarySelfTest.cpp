#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
#include "Boundaries/FluxBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"

FILE * flux_f; 

class FluxAndPeriodicBoundarySelfTest : public Mercury2D
{
    public:
        FluxAndPeriodicBoundarySelfTest()
        {
            LinearViscoelasticFrictionSpecies* speciesP;

            setTimeStep(1e-3);
            setTimeMax(1.5);
            setSaveCount(100);
            
            speciesP = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
            speciesP->setDensity(1);

            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    5e-1,
                    0.5,
                    (4.0/3.0)*constants::pi*pow(0.05, 3)*1
            );

            speciesP->setSlidingFrictionCoefficient(tan(34 * constants::pi / 180));
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(0);
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(0);
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());

            perb = new PeriodicBoundary;
            perb->set(Vec3D(1,0,0), -1, 1);
            perb = boundaryHandler.copyAndAddObject(perb);

            fluxb = new FluxBoundary;
            fluxb->set(Vec3D(1,0,0), 0);
            fluxb = boundaryHandler.copyAndAddObject(fluxb);

            BaseParticle p;
            p.setSpecies(speciesP);
            for (int i = -20; i <= 20; i++)
            {
                p.setRadius( 0.5 + 0.01*i );
                p.setPosition(Vec3D(0.02*i, i, 0));
                p.setVelocity(Vec3D(i, 0, 0));
                particleHandler.copyAndAddObject(p);
            }

            /* Output options */
            setXBallsAdditionalArguments("-cmode 7 -v0");
            setXMin(0);
            setXMax(1);
            setYMin(-1);
            setYMax(1);
            setZMin(-1);
            setZMax(1);

            dataFile.setFileType(FileType::ONE_FILE);
            restartFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setGravity(Vec3D(1, 0, 0));
        }

        ~FluxAndPeriodicBoundarySelfTest() override = default;

    void actionsAfterTimeStep() override {
            fprintf(flux_f, "%f %f %d %f %f %d %f %f\n", 
                    getTime(),
                    particleHandler.getMass(),
                    fluxb->getNumberOfParticlesCrossedForw(),
                    fluxb->getMassOfParticlesCrossedForw(),
                    fluxb->getVolumeOfParticlesCrossedForw(),
                    fluxb->getNumberOfParticlesCrossedBack(),
                    fluxb->getMassOfParticlesCrossedBack(),
                    fluxb->getVolumeOfParticlesCrossedBack() );
        }

    private:
        FluxBoundary* fluxb;
        PeriodicBoundary* perb;
};

int main(int argc, char *argv[])
{
    flux_f = fopen("FluxAndPeriodicBoundarySelfTest.flux", "w");
    fprintf(flux_f, "time mass flux_num_forw flux_mass_forw flux_vol_forw flux_num_back flux_mass_back flux_vol_back\n");
    setbuf(flux_f, nullptr);
    FluxAndPeriodicBoundarySelfTest problem;
    problem.setName("FluxAndPeriodicBoundarySelfTest");
    fprintf(stdout, "Initialising the problem\n");
    problem.solve();
    fclose(flux_f);
}