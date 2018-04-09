/* A test of insertion, flux and deletion boundaries. In particular, can test if
 * they allocate and free memory properly, when used with valgrind. */

#include "Mercury2D.h"
#include "Particles/BaseParticle.h"
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

FILE * flux_f; 

class BoundariesSelfTest : public Mercury2D
{
    public:
        BoundariesSelfTest()
        {
            LinearViscoelasticFrictionSpecies* speciesP;

            setTimeStep(1e-3);
            setTimeMax(3.00);
            setSaveCount(1);
            
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

            /* Introduce the InsertionBoundary */
            insb = new CubeInsertionBoundary;
            BaseParticle* insertionBoundaryParticle = new BaseParticle; // Possibly evil!
            insertionBoundaryParticle->setSpecies(speciesP);
            insb->set(
                    insertionBoundaryParticle,
                    10,
                    Vec3D(0, 0, 0),
                    Vec3D(0, 1, 0),
                    Vec3D(-0.2,0,0),
                    Vec3D(0.2,0,0),
                    0.05 * 0.9, 
                    0.05 * 1.1 
            );
            insb = boundaryHandler.copyAndAddObject(insb);
            not_yet_deleted_insb = true;

            /* Introduce the DeletionBoundary */
            delb = new DeletionBoundary;
            delb->set(Vec3D(1,0,0), 1);
            delb = boundaryHandler.copyAndAddObject(delb);

            /* Introduce the FluxBoundary */
            fluxb = new FluxBoundary;
            fluxb->set(Vec3D(1,0,0), 0.5);
            fluxb = boundaryHandler.copyAndAddObject(fluxb);

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
            fStatFile.setFileType(FileType::NO_FILE);

            /* Set gravity */
            double g = 1;
            setGravity(Vec3D(g, 0, 0));
        }

        ~BoundariesSelfTest()
        {
        }

        void actionsAfterTimeStep()
        {
            if (not_yet_deleted_insb)
            {
                num_inserted = insb->getNumberOfParticlesInserted();
                mass_inserted = insb->getMassOfParticlesInserted();
                vol_inserted = insb->getVolumeOfParticlesInserted();
                if (getTime() > 1)
                {
                    boundaryHandler.removeObject(insb->getIndex());
                    not_yet_deleted_insb = false;
                }
            }

            fprintf(flux_f, "%f %f %d %f %f %d %f %f %d %f %f\n", 
                    getTime(),
                    particleHandler.getMass(),
                    num_inserted,
                    mass_inserted,
                    vol_inserted,
                    fluxb->getNumberOfParticlesCrossedNet(),
                    fluxb->getMassOfParticlesCrossedNet(),
                    fluxb->getVolumeOfParticlesCrossedNet(),
                    delb->getNumberOfParticlesDeleted(),
                    delb->getMassOfParticlesDeleted(),
                    delb->getVolumeOfParticlesDeleted() );
        }

    private:
        CubeInsertionBoundary* insb;
        DeletionBoundary* delb;
        FluxBoundary* fluxb;

        bool not_yet_deleted_insb;
        unsigned int num_inserted; 
        double mass_inserted, vol_inserted;
};

int main(int argc, char *argv[])
{
    flux_f = fopen("BoundariesSelfTest.flux", "w");
    fprintf(flux_f, "time mass insb_num insb_mass insb_vol fluxb_num fluxb_mass fluxb_vol delb_num delb_mass delb_vol\n");
    setbuf(flux_f, NULL);
    BoundariesSelfTest problem;
    problem.setName("BoundariesSelfTest");
    fprintf(stdout, "Initialising the problem\n");
    problem.solve();
    fclose(flux_f);
}
