//
// Created by paolo on 13-12-19.
//

#include <iostream>
#include "Mercury3D.h"
#include <Boundaries/RandomClusterInsertionBoundary.h>
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include <Walls/InfiniteWall.h>

/*!
 * In this the RandomClusterInsertionBoundary is used: in particular here radius cluster and radius particle will be
 *      used to create the cluster. Knowing those two parameters is enough to know the number of particles as
 *      the number of particles (N) per cluster given the relative cluster radius (hatR) and penetration depth max (phi)
        can be computed as: N = ( hatR / (1 - eps0*phi) )^3 * eps0, with eps0 being the mass fraction in the limit
        of penetration depth max -> 0.
        Final mass fraction can be computed as eps0 + 3*pow(eps0,2)*phi.
        IMPORTANT: There will be accordance among this formula and the obtained values
        as long as the user will set friction parameters similar to the ones set here.

 */
class RandomClusterInsertionBoundarySelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        //setting standard variables
        setName("ClusterInsertionBoundaryUnitTest");
        fStatFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        setXBallsAdditionalArguments("-solidf -v0");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, -0.0981));
        setTimeStep(1e-4);
        dataFile.setSaveCount(60);
        setTimeMax(0.5);
        setHGridMaxLevels(2);

        setMin(Vec3D(-0.4, -0.05, -0.4));
        setMax(Vec3D(0.4, 0.05, 0.4));


        //creating species
        Mdouble loadingStiffness = 1e3;
        Mdouble radiusParticle = 0.03;
        Mdouble restitutionCoefficient = 0.5;
        Mdouble constantMass = 2700 * 4 * constants::pi * pow(radiusParticle, 3) / 3;
        Mdouble collisionTimeIntra = sqrt(constantMass * (pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2)) /
                                          (2 * loadingStiffness));

        LinearPlasticViscoelasticFrictionSpecies* species;
        species = new LinearPlasticViscoelasticFrictionSpecies;
        species->setConstantRestitution(true);
        species->setDensity(2700);
        species->setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, 0.5, 1);
        species->setUnloadingStiffnessMax(species->getLoadingStiffness() * 5);
        species->setCohesionStiffness(species->getLoadingStiffness());
        species->setPenetrationDepthMax(0.1);

        species->setSlidingFrictionCoefficient(0.5);
        species->setSlidingStiffness(species->getLoadingStiffness() * 2.0/7.0);
        species->setSlidingDissipation(species->getDissipation() * 2.0 / 7.0);
        species->setRollingFrictionCoefficient(0.3);
        species->setRollingStiffness(species->getLoadingStiffness() * 2.0/7.0);
        species->setRollingDissipation(species->getDissipation() * 2.0 / 7.0);
        speciesHandler.copyAndAddObject(species);

        //spherical particle needed for the insertion
        SphericalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        //setting insertion boundary properties
        RandomClusterInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle, 100, {-0.4, -0.01, 0.2}, {-0.2, 0.01, 0.4}, Vec3D(0.2, 0, 0),
                              Vec3D(0.2, 0, 0), 0.05, 0.05, radiusParticle);
        insertionBoundary.setRandomised(false);
        boundaryHandler.copyAndAddObject(insertionBoundary);

    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }
};

int main()
{
    logger(INFO,"Simple box for creating particles");

    RandomClusterInsertionBoundarySelfTest insertionBoundary_problem;
    insertionBoundary_problem.solve();


    helpers::check(24, insertionBoundary_problem.particleHandler.getSize(), 0.1, "Number of particles check");
    helpers::check(26, insertionBoundary_problem.interactionHandler.getSize(), 0.1, "Number of interactions check");
}

