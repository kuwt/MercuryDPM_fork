//ciao!

#include <iostream>
#include "Mercury3D.h"
#include <Boundaries/RandomClusterInstertionBoundary.h>
#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include <Walls/InfiniteWall.h>


class InsertionBoundarySelfTest1 : public Mercury3D
{
public:

    void setupInitialConditions() override {
        setName("InsertionBoundarySelfTest");
        fStatFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        setXBallsAdditionalArguments("-solidf -v0");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, -0.0981));
        setTimeStep(1e-4);
        dataFile.setSaveCount(60);
        setTimeMax(50);
        setHGridMaxLevels(2);

        setMin(Vec3D(-0.4, -0.05, -0.4));
        setMax(Vec3D(0.4, 0.05, 0.4));


        Mdouble constantMass = 2700 * 4 * constants::pi * pow(0.013, 3) / 3;
        Mdouble collisionTimeIntra = sqrt(constantMass * (pow(constants::pi, 2) + pow(log(0.5), 2)) /
                                          (2 * 1e3));

        LinearPlasticViscoelasticFrictionSpecies* species;
        species = new LinearPlasticViscoelasticFrictionSpecies;
        species->setConstantRestitution(true);
        species->setDensity(2700);
        species->setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, 0.5, 1);
        species->setUnloadingStiffnessMax(species->getLoadingStiffness() * 5);
        species->setCohesionStiffness(species->getUnloadingStiffnessMax() / 3);
        species->setPenetrationDepthMax(0.1);

        species->setSlidingFrictionCoefficient(0.1);
        species->setSlidingStiffness(species->getLoadingStiffness() * 2.0/7.0);
        species->setSlidingDissipation(species->getDissipation() * 2.0 / 7.0);
        speciesHandler.copyAndAddObject(species);


        LinearPlasticViscoelasticFrictionSpecies* speciesWall;
        speciesWall = new LinearPlasticViscoelasticFrictionSpecies;
        speciesWall->setConstantRestitution(true);
        speciesWall->setDensity(2700);
        speciesWall->setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, 1, 1);
        speciesWall->setUnloadingStiffnessMax(speciesWall->getLoadingStiffness());
        speciesWall->setCohesionStiffness(0);
        speciesWall->setPenetrationDepthMax(0.1);

        speciesWall->setSlidingFrictionCoefficient(0);
        speciesWall->setSlidingStiffness(0);
        speciesWall->setSlidingDissipation(0);
        speciesHandler.copyAndAddObject(speciesWall);

        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setLoadingStiffness(species->getLoadingStiffness());
        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setUnloadingStiffnessMax(species->getLoadingStiffness());
        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setCohesionStiffness(0);
        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setSlidingFrictionCoefficient(0);

        InfiniteWall wall;
        wall.set({0,0,-1}, {0,0,-0.4});
        wall.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(wall);
        wall.set({1,0,0}, {0.4,0,0});
        wall.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(wall);
        wall.set({-1,0,0}, {-0.4,0,0});
        wall.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(wall);
        wall.set({0,1,0}, {0.0,0.1,0.0});
        wall.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(wall);
        wall.set({0,-1,0}, {0.0,-0.1,0.0});
        wall.setSpecies(speciesHandler.getObject(1));
        wallHandler.copyAndAddObject(wall);


        SphericalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));


        RandomClusterInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle, 100, {-0.4, -0.01, 0.2}, {-0.2, 0.01, 0.4}, Vec3D(0.2, 0, 0),
                              Vec3D(0.2, 0, 0), 0.05, 0.05, 0.0125);
        //insertionBoundary.set(&insertionBoundaryParticle,100,{0.2, -0.05, 0.2},{0.4, 0.05, 0.4},Vec3D(0.2,0,0),Vec3D(0.2,0,0),0.05,0.05, 0.0125);
        //insertionBoundary.checkBoundaryBeforeTimeStep(this);
        boundaryHandler.copyAndAddObject(insertionBoundary);

    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Simple box for creating particles");

    InsertionBoundarySelfTest1 insertionBoundary_problem;
    insertionBoundary_problem.setMin(Vec3D(-0.5, -0.5, -0.5));
    insertionBoundary_problem.setMax(Vec3D(0.5, 0.5, 0.5));



    Mdouble constantMass = 2700 * 4 * constants::pi * pow(0.013, 3) / 3;
    Mdouble collisionTimeIntra = sqrt(constantMass * (pow(constants::pi, 2) + pow(log(0.5), 2)) /
    (2 * 1e3));

    LinearPlasticViscoelasticFrictionSpecies* species;
    species = new LinearPlasticViscoelasticFrictionSpecies;
    species->setConstantRestitution(true);
    species->setDensity(2700);
    species->setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, 0.5, 1);
    species->setUnloadingStiffnessMax(
    species->getLoadingStiffness() * 5);
    species->setCohesionStiffness(species->getUnloadingStiffnessMax() / 3);
    species->setPenetrationDepthMax(0.1);

    species->setSlidingFrictionCoefficient(0.1);
    species->setSlidingStiffness(species->getLoadingStiffness() * 2.0/7.0);
    species->setSlidingDissipation(species->getDissipation() * 2.0 / 7.0);
    insertionBoundary_problem.speciesHandler.copyAndAddObject(species);





    insertionBoundary_problem.setNumberOfDomains({1,2,1});
    insertionBoundary_problem.solve();
}
