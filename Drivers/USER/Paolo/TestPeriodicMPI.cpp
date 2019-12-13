//
// Created by paolo on 4-9-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>
#include <Particles/SphericalParticle.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>

class Compression : public Mercury3D
{
public:



    //! Initializing variables, boxSize, species, radii, walls, particles, timeStep, output, name
    void setupInitialConditions() override {

        /*
         * ----------------------------------------------------------------------
         *                        DEFINING VARIABLES
         * ----------------------------------------------------------------------
         */

        random.randomise();

        //! Strain velocity
        isoStrainDot = -1e-8 / getTimeStep();

        //!CREATING BOUNDARIES
        createBoundaries();

        //!SETTING PARTICLES POSITIONS
        setParticlePositions();

        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(8e-6);
        p0.setVelocity({0,0,0});

        for (int i = 0; i < nParticles; ++i) {
            p0.setPosition(particlePositions[i]);
            particleHandler.copyAndAddObject(p0);
        }

    }

    //! Makes data analysis and writes to isotropic compression file
    void actionsBeforeTimeStep() override
    {
        moveBoundaries();
        computeStress();

        if ( uniStress > 5e5 ) {
            setTimeMax(getTime());
        }
    }

    //! Prints time and values of interest
    void printTime() const override
    {
        std::cout <<

                  "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
                  ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
                  ", Stress: " << std:: fixed << std::setprecision(2) << std::setw(9) << totalStress <<
                  ", uniStress: " << std:: fixed << std::setprecision(2) << std::setw(9) << uniStress <<
                  ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy()

                  << std::endl;
    }

    //! inserts the particles in the simulation domainss
    void setParticlePositions()
    {
        particlePositions = {
                {-4.12471e-06, 1.50160e-10, -4.04206e-06},
                {1.17580e-05, -1.09613e-06, -2.49174e-06},
                {-1.17524e-05, 1.17537e-05, -1.17564e-05},
                {-1.11593e-05, -1.17501e-05, 4.23077e-06},
                {-2.22116e-07, -1.32469e-06, 1.14160e-05},
                {4.21677e-06, 1.17523e-05, 3.33631e-06},
                {-1.17517e-05, 1.13704e-05, 4.23706e-06},
                {4.24337e-06, -1.17532e-05, -1.17548e-05},
                {-1.17544e-05, -1.17532e-05, -1.17562e-05},
                {1.17504e-05, -1.17523e-05, 9.43951e-06},
                {9.52649e-06, 1.17528e-05, -1.17554e-05}
        };
    }

    //! Creating stress strain control walls
    void createBoundaries()
    {
        boundaryPositionX = 3.95e-05/2;
        boundaryPositionY = 3.95e-05/2;
        boundaryPositionZ = 3.95e-05/2;

        boundary.set(Vec3D(0, 0, 1), -boundaryPositionZ, boundaryPositionZ);
        BNS = boundaryHandler.copyAndAddObject(boundary);

        boundary.set(Vec3D(1, 0, 0), -boundaryPositionX, boundaryPositionX);
        BEW = boundaryHandler.copyAndAddObject(boundary);

        boundary.set(Vec3D(0, 1, 0), -boundaryPositionY, boundaryPositionY);
        BFB = boundaryHandler.copyAndAddObject(boundary);

    }

    //! Moves periodic boundaries according to isoStrainDot
    void moveBoundaries()
    {
        boundaryPositionX *= 1+isoStrainDot*getTimeStep();
        boundaryPositionY *= 1+isoStrainDot*getTimeStep();
        boundaryPositionZ *= 1+isoStrainDot*getTimeStep();

        volumeBox = boundaryPositionX * boundaryPositionY * boundaryPositionZ * 8;

        BNS->set(Vec3D(0, 0, 1), -boundaryPositionZ, boundaryPositionZ);

        BEW->set(Vec3D(1, 0, 0), -boundaryPositionX, boundaryPositionX);

        BFB->set(Vec3D(0, 1, 0), -boundaryPositionY, boundaryPositionY);
    }

    //! computes stress
    void computeStress()
    {
        stressXX = 0.0;
        stressYY = 0.0;
        stressZZ = 0.0;
        stressXY = 0.0;
        stressXZ = 0.0;
        stressYZ = 0.0;

        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {
            stressXX += (*i)->getForce().X * (*i)->getNormal().X * (*i)->getDistance();
            stressYY += (*i)->getForce().Y * (*i)->getNormal().Y * (*i)->getDistance();
            stressZZ += (*i)->getForce().Z * (*i)->getNormal().Z * (*i)->getDistance();
            stressXY += (*i)->getForce().X * (*i)->getNormal().Y * (*i)->getDistance();
            stressXZ += (*i)->getForce().X * (*i)->getNormal().Z * (*i)->getDistance();
            stressYZ += (*i)->getForce().Y * (*i)->getNormal().Z * (*i)->getDistance();
        }

        stressXX /= volumeBox;
        stressYY /= volumeBox;
        stressZZ /= volumeBox;
        stressXY /= volumeBox;
        stressXZ /= volumeBox;
        stressYZ /= volumeBox;

        totalStress = (stressXX+stressYY+stressZZ)/3;

        uniStress = stressZZ;

    }

    /*
     * -------------------------- VARIABLES ------------------------------------
     */


    // Simulation
    Mdouble cdatOutputTimeInterval;

    // domain
    Mdouble volumeBox;
    Mdouble initialSolidFraction;


    // Particles
    Mdouble radiusParticle;
    Mdouble volumeParticle;
    Mdouble massParticle;
    int nParticles;

    // species
    LinearPlasticViscoelasticFrictionSpecies* species;
    Mdouble collisionTimeSmallestMass;

    Mdouble densityParticle;

    Mdouble loadingStiffness;
    Mdouble unLoadingStiffnessMax;
    Mdouble cohesionStiffness;
    Mdouble slidingFrictionCoefficient;
    Mdouble slidingRollingCoefficient;
    Mdouble slidingTorsionCoefficient;
    Mdouble restitutionCoefficient;
    Mdouble penetrationDepthMax;


private:

    Mdouble isoStrainDot;
    std::vector<Vec3D> particlePositions;

    // Walls
    Mdouble boundaryPositionX;
    Mdouble boundaryPositionY;
    Mdouble boundaryPositionZ;
    PeriodicBoundary boundary;
    PeriodicBoundary* BNS;
    PeriodicBoundary* BEW;
    PeriodicBoundary* BFB;

    // Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble totalStress;
    Mdouble uniStress;

};



int main(){

    Compression script;

    script.dataFile.setFileType(FileType::ONE_FILE);
    script.restartFile.setFileType(FileType::ONE_FILE);
    script.fStatFile.setFileType(FileType::NO_FILE);
    script.eneFile.setFileType(FileType::NO_FILE);
    logger(INFO, "run number: %", script.dataFile.getCounter());

    Mdouble boxSize = 3.95e-05;
    script.setXMin(-0.5*boxSize);
    script.setYMin(-0.5*boxSize);
    script.setZMin(-0.5*boxSize);
    script.setXMax( 0.5*boxSize);
    script.setYMax( 0.5*boxSize);
    script.setZMax( 0.5*boxSize);


    // Particles properties
    script.nParticles = 11;
    script.radiusParticle =8e-6;
    script.volumeParticle = 4 * constants::pi * pow(script.radiusParticle, 3) / 3;



    //!Species parameters
    script.densityParticle = 2825; // ottenuta considerando la densità della bentonite = 2150 Kg/m^3 ed una vF di 0.37 (corrispondente al momento in cui E_ratio crolla)
    script.loadingStiffness = 1e3;
    script.unLoadingStiffnessMax = 5 * script.loadingStiffness;
    script.cohesionStiffness = 0.5 * script.loadingStiffness;
    script.slidingFrictionCoefficient = 0.0;
    script.slidingRollingCoefficient = 0.5 * script.slidingFrictionCoefficient;
    script.slidingTorsionCoefficient = 0.0;
    script.restitutionCoefficient = 0.5;
    script.penetrationDepthMax = 0.1;
    script.massParticle = script.volumeParticle*script.densityParticle;
    script.collisionTimeSmallestMass = sqrt( script.massParticle * ( pow(constants::pi, 2) + pow(log(script.restitutionCoefficient), 2) ) / ( 2 * script.loadingStiffness ) );


    // SINGLE_PARTICLE-SINGLE_PARTICLE (set here, added in class)

    script.species = new LinearPlasticViscoelasticFrictionSpecies;
    script.species -> setDensity(script.densityParticle);
    script.species -> setConstantRestitution(true);
    script.species -> setCollisionTimeAndRestitutionCoefficient(script.collisionTimeSmallestMass, script.restitutionCoefficient, 1);
    script.species -> setUnloadingStiffnessMax(script.species->getLoadingStiffness() * script.unLoadingStiffnessMax / script.loadingStiffness);
    script.species -> setCohesionStiffness(script.species-> getLoadingStiffness() * script.cohesionStiffness / script.loadingStiffness);
    script.species -> setPenetrationDepthMax(script.penetrationDepthMax);

    script.species -> setSlidingFrictionCoefficient(script.slidingFrictionCoefficient);
    script.species -> setSlidingStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setSlidingDissipation(script.species -> getDissipation()*2.0/7.0);
    script.species -> setRollingFrictionCoefficient(script.slidingRollingCoefficient);
    script.species -> setRollingStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setRollingDissipation(script.species -> getDissipation()*2.0/7.0);
    script.species -> setTorsionFrictionCoefficient(script.slidingTorsionCoefficient);
    script.species -> setTorsionStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setTorsionDissipation(script.species -> getDissipation()*2.0/7.0);

    script.speciesHandler.copyAndAddObject(script.species);


    Mdouble contactTimeOverTimeStep = 51;

    script.setTimeStep(script.species -> getCollisionTime(script.massParticle)/contactTimeOverTimeStep);

    script.cdatOutputTimeInterval = script.getTimeStep() * 5000;
    script.setSaveCount( floor(script.cdatOutputTimeInterval/script.getTimeStep()) );
    script.setTimeMax(1000);
    script.setXBallsAdditionalArguments("-v0 -p 10");

    // NAME SETTING
    std::ostringstream name;
    name << "TestPeriodicMPI";
    script.setName(name.str());

    script.setNumberOfDomains({2,2,2});

    script.solve();

    return 0;
}




