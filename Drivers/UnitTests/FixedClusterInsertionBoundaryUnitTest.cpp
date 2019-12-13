//
// Created by paolo on 13-12-19.
//

#include <iostream>
#include "Mercury3D.h"
#include "Boundaries/FixedClusterInstertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


/*!
 * In this the FixedClusterInsertionBoundary is used: in order to create clusters, radii and positions are needed.
 *      This works only creating clusters imposing radius cluster and radius particle.
 *      Knowing those two parameters is enough to know the number of particles as
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
        setName("FixedClusterInsertionBoundaryUnitTest");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, 0));
        setXBallsAdditionalArguments("-solidf -v0");
        setTimeStep(1e-4);
        dataFile.setSaveCount(10);
        setTimeMax(0.5);
        setHGridMaxLevels(2);

        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(0.4,0.4,0.4));

        //creating species
        Mdouble densityParticle = 2000;
        Mdouble loadingStiffness = 1e3;
        Mdouble radiusParticle = 0.03;
        Mdouble restitutionCoefficient = 0.5;
        Mdouble constantMass = densityParticle * 4 * constants::pi * pow(radiusParticle, 3) / 3;
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


        SphericalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        std::vector<Mdouble> radii;
        radii.push_back(0.05);
        radii.push_back(0.05);
        radii.push_back(0.05);

        std::vector<Vec3D> pos;
        pos.reserve(3);
        pos.push_back(Vec3D(0.1, 0.2, 0.1));
        pos.push_back(Vec3D(0.2, 0.2, 0.2));
        pos.push_back(Vec3D(0.1, 0.2, 0.3));

        FixedClusterInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle, pos, radii, Vec3D(1, 0, 0), Vec3D(1, 0, 0), radiusParticle);
        insertionBoundary.setRandomised(false);
        //Here the insertion boundary is not added to the boundary handler: this way clusters will be inserted just one time,
        // at the beginning of the simulation
        insertionBoundary.checkBoundaryBeforeTimeStep(this);

    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Fixed cluster insertion boundary.");

    RandomClusterInsertionBoundarySelfTest insertionBoundary_problem;
    insertionBoundary_problem.solve();

    helpers::check(9, insertionBoundary_problem.particleHandler.getSize(), 0.1, "Number of particles check");
    helpers::check(9, insertionBoundary_problem.interactionHandler.getSize(), 0.1, "Number of interactions check");
}
