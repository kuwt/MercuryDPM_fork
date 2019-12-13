//
// Created by paolo on 9-9-19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>
#include <Particles/SphericalParticle.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <fstream>
#include <Boundaries/StressStrainControlBoundary.h>

class Compression : public Mercury3D
{


private:

    //! Initializing variables, boxSize, species, radii, walls, particles, timeStep, output, name
    void setupInitialConditions() override {

        /*
         * ----------------------------------------------------------------------
         *                        DEFINING VARIABLES
         * ----------------------------------------------------------------------
         */

        //random.randomise();

        readRestartFile("BaseSinglePorosityR2");

        std::cout << std::endl << std::endl << "SETTING SPECIES" << std::endl << std::endl;
        setSpecies();

        // Strain
        isoStrainDot = -1e-8/getTimeStep(); // should be something like 1e-8/timeStep but for the first part it is set faster

        totalParticleVolume = particleHandler.getVolume();
        std::cout << std::endl << std::endl << "CREATING WALLS" << std::endl << std::endl;
        createWalls();

        t0 = getTime();

        cdatOutputTimeInterval = getTimeStep() * 3000;

        setSaveCount( floor(cdatOutputTimeInterval/getTimeStep()) );

        setTimeMax(1e5*cdatOutputTimeInterval);

        setXBallsAdditionalArguments("-v0 -p 10");

        // NAME SETTING
        std::ostringstream name;
        name << "SinglePorosityR2S3F1";
        setName(name.str());

        makeIsoCFile();

        stage = 1;

    }

    //! Makes data analysis and writes to isotropic compression file
    void actionsBeforeTimeStep() override
    {
        if (stage == 1) {
            makeDataAnalysis();

            //dampVelocities();
            if (fmod(getTime() - t0, cdatOutputTimeInterval) < getTimeStep()) {
                writeToIsoCFile();
            }

            if ( uniStress > 5e5)
            {
                boundaryHandler.clear();
                w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), true);
                boundaryHandler.copyAndAddObject(w0);
                stage = 2;
            }
        }

        if (stage == 2) {
            makeDataAnalysis();

            //dampVelocities();
            if (fmod(getTime() - t0, cdatOutputTimeInterval) < getTimeStep()) {
                writeToIsoCFile();
            }

            if ( uniStress > 2e6 )
                setTimeMax(getTime());
        }
    }

    //! Prints time and values of interest
    void printTime() const override
    {

        switch (stage)
        {
            case 1: std::cout << "Isotropic compression: ";
                break;

            case 2: std::cout << "Isotropic compression: ";
                break;

            case 3: std::cout << "Uniaxial compression: ";
                break;

            default: std::cout << "Final values: ";
                break;
        }
        std::cout <<

                  "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
                  ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
                  ", Stress: " << std:: fixed << std::setprecision(2) << std::setw(9) << totalStress <<
                  ", uniStress: " << std:: fixed << std::setprecision(2) << std::setw(9) << uniStress <<
                  ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
                  ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber <<
                  ", vR = " << std::fixed << std::setprecision(3) << voidRatio <<
                  ", Porosity = " << std::fixed << std::setprecision(3) << e <<
                  ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlap <<
                  ", dMean = " << std:: fixed << std::setw(7) << meanRelativeOverlap <<
                  ", dAvCu = " << std:: fixed << std::setw(7) << averageCubicRelativeOverlap <<
                  ", dMax = " << maxRelativeOverlap

                  << std::endl;
    }

    //! Creating stress strain control walls
    void createWalls()
    {
        boundaryHandler.clear();
        w0.setHandler(&boundaryHandler);
        w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(isoStrainDot, 0, 0, 0, isoStrainDot, 0, 0, 0, isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), true);
        boundaryHandler.copyAndAddObject(w0);
    }

    void setSpecies(){

        speciesHandler.clear();

        Mdouble densityParticle = 2825; // ottenuta considerando la densità della bentonite = 2150 Kg/m^3 ed una vF di 0.37 (corrispondente al momento in cui E_ratio crolla)

        Mdouble loadingStiffness = 1e3;
        Mdouble slidingFrictionCoefficient = 0.0;
        Mdouble slidingRollingCoefficient = slidingFrictionCoefficient/2;
        Mdouble slidingTorsionCoefficient = 0.0;
        Mdouble restitutionCoefficient = 0.5;
        Mdouble penetrationDepthMax = 0.5;


        Mdouble contactTimeOverTimeStep = 51;

        // Particles properties
        radiusParticle =5.5e-5;
        dispersity = 0.1;

        Mdouble smallestMass = densityParticle * 4 * constants::pi * pow(radiusParticle*(1-dispersity), 3) / 3;
        Mdouble collisionTimeSmallestMass = sqrt( smallestMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );

        // Species
        LinearPlasticViscoelasticFrictionSpecies* species;


        // SINGLE_PARTICLE-SINGLE_PARTICLE

        species = new LinearPlasticViscoelasticFrictionSpecies;
        species -> setDensity(densityParticle);
        species -> setConstantRestitution(true);
        species -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
        species -> setUnloadingStiffnessMax(species->getLoadingStiffness());
        species -> setCohesionStiffness(0);
        species -> setPenetrationDepthMax(penetrationDepthMax);

        species -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
        species -> setSlidingStiffness(species -> getLoadingStiffness()*2.0/7.0);
        species -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
        species -> setRollingFrictionCoefficient(slidingRollingCoefficient);
        species -> setRollingStiffness(species -> getLoadingStiffness()*2.0/7.0);
        species -> setRollingDissipation(species -> getDissipation()*2.0/7.0);
        species -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
        species -> setTorsionStiffness(species -> getLoadingStiffness()*2.0/7.0);
        species -> setTorsionDissipation(species -> getDissipation()*2.0/7.0);

        speciesHandler.copyAndAddObject(species);


        std::cout << "collision time: " << species -> getCollisionTime(smallestMass) << std::endl;
        setTimeStep(species -> getCollisionTime(smallestMass)/contactTimeOverTimeStep);
        // prints the timeStep
        std::cout << "timeStep: " << std::setprecision(4) << getTimeStep() << std::endl;
        // prints the ratio between collision time and time step for a pure particle-particle interaction (should be > 50 to be precise, it is 51)
        std::cout << "tC_PP/dt, at least: " << std::setprecision(4) << species -> getCollisionTime(smallestMass)/getTimeStep() << std::endl << std::endl;



        for (int i = 0; i < particleHandler.getNumberOfRealObjects(); ++i) {
            particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(0));
        }
    }

    //! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
    void makeDataAnalysis()
    {
        //! Computes coordination number and overlaps.
        computeCoordinationNumberAndOverlaps();

        //! Computes void ratio (and volumeBox)
        computeVoidRatioAndPorosity();

        //! Computes the stress.
        computeStress();
    }

    //! Computes coordination number and overlaps
    void computeCoordinationNumberAndOverlaps()
    {
        meanForceOnInteraction = 0.0;
        meanCoordinationNumber = 0.0;
        maxRelativeOverlap = 0.0;
        averageCubicRelativeOverlap = 0.0;
        meanRelativeOverlap = 0.0;
        minRelativeOverlap = 2 * ( 1 + dispersity );
        Mdouble relativeOverlap;

        //! loops over each particle to compute mean coordination number and center of mass.
        for (int i=0; i < particleHandler.getNumberOfRealObjects(); i++)
        {
            meanCoordinationNumber += (particleHandler.getObject(i)->getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfRealObjects();

        //! loops over every interaction to compute mean force acting on interaction, maximum, mean and minimum relative particle overlap.
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            meanForceOnInteraction += ((*i) -> getForce()).getLength();

            /*!
             * \details the relative overlap is computed as an average of the relative overlap on the two particles.
             *          rO = ( O/R1 + O/R2 ) / 2.
             */
            relativeOverlap = ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                              ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
            relativeOverlap /= 2;
            meanRelativeOverlap += relativeOverlap;
            averageCubicRelativeOverlap += pow( relativeOverlap, 3 );
            if (relativeOverlap > maxRelativeOverlap)
                maxRelativeOverlap = relativeOverlap;

            if (relativeOverlap < minRelativeOverlap)
                minRelativeOverlap = relativeOverlap;
        }
        meanForceOnInteraction/= interactionHandler.getSize();
        meanRelativeOverlap /= interactionHandler.getSize();
        averageCubicRelativeOverlap = cbrt( averageCubicRelativeOverlap);
        averageCubicRelativeOverlap /= interactionHandler.getSize();
    }

    //! computes voidRatio and porosity (e)
    void computeVoidRatioAndPorosity()
    {
        volumeBox = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
        voidRatio = 1 - totalParticleVolume / volumeBox;
        e = (volumeBox - totalParticleVolume) / totalParticleVolume;
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

    //! makes a .isoc file in wich all ifo of isotropic compression are saved
    void makeIsoCFile()
    {

        std::ostringstream cdatName;
        cdatName << getName() << ".isoc";

        isoCFile.open(cdatName.str(), std::ios::out);

        isoCFile << "CLUSTER DATA AND INFORMATION" << std::endl << std::endl;
        isoCFile << "radiusParticle: " << std:: scientific << std::setprecision(2) << radiusParticle << std::endl;
        isoCFile << "sizeDispersityParticle: " << std::defaultfloat << dispersity << std::endl;
        isoCFile << "boxSize: " << std::defaultfloat << (getXMax()-getXMin()) << std::endl;
        isoCFile << "initial solid fraction: " << std::defaultfloat << initialSolidFraction << std::endl;
        // If constantRestitution(true) loading, unloading, and cohesion stiffness are multiplied by the mass of a particle whose radius is radiusParticle_*(1-sizeDispersity),
        // which is the mass that should be used to compute collision time.
        isoCFile << "velocityDampingModulus: " << std::defaultfloat << velocityDampingModulus << std::endl << std::endl;

        isoCFile << "COMPRESSION INFORMATION" << std::endl << std::endl;
        isoCFile <<

                 "   Time" << std::setw(13) <<
                 "Time Max" << std::setw(16) <<
                 "Stress" << std::setw(16) <<
                 "uniStress" << std::setw(15) <<
                 "E_ratio" << std::setw(12) <<
                 "cN" << std::setw(12) <<
                 "vR" << std::setw(12) <<
                 "e" << std::setw(14) <<
                 "dMin" << std::setw(14) <<
                 "dMean" << std::setw(14) <<
                 "dMax"

                 << std::endl;
    }

    //! writes on the .isoc file
    void writeToIsoCFile()
    {
        isoCFile <<

                 std::setprecision(2) << std::setw(7) << getTime() <<
                 std::setprecision(2) << std::setw(13) << getTimeMax() <<
                 std::fixed << std::setprecision(2) << std::setw(16) << totalStress <<
                 std::fixed << std::setprecision(2) << std::setw(16) << uniStress <<
                 std::scientific << std::setprecision(2) << std::setw(15) << getKineticEnergy()/getElasticEnergy() <<
                 std::fixed << std::setw(12) << meanCoordinationNumber <<
                 std::fixed << std::setprecision(3) << std::setw(12) << voidRatio <<
                 std::fixed << std::setprecision(3) << std::setw(12) << e <<
                 std::fixed << std::setprecision(5) << std::setw(14) << minRelativeOverlap <<
                 std::setw(14) << averageCubicRelativeOverlap <<
                 std::setw(14) << maxRelativeOverlap

                 << std::endl;

    }

    /*
     * -------------------------- VARIABLES ------------------------------------
     */


    // Simulation
    int stage;
    Mdouble t0;
    Mdouble cdatOutputTimeInterval;


    // Domain
    Mdouble isoStrainDot;
    Mdouble volumeBox;
    Mdouble initialSolidFraction;

    // Particles
    Mdouble volumeParticle;
    Mdouble massParticle;
    Mdouble collisionTime;
    Mdouble nParticles;
    Mdouble velocityDampingModulus;
    Mdouble totalParticleVolume;
    std::vector<double> radii;
    Mdouble dispersity;

    // Walls
    StressStrainControlBoundary w0;

    // Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble totalStress;
    Mdouble uniStress;

    // Data analysis
    Mdouble  meanForceOnInteraction;
    Mdouble  meanCoordinationNumber;
    Mdouble  maxRelativeOverlap;
    Mdouble  averageCubicRelativeOverlap;
    Mdouble  meanRelativeOverlap;
    Mdouble  minRelativeOverlap;
    Mdouble  voidRatio;
    Mdouble  e;

    // Output
    //!\brief cluster data file.
    std::ofstream isoCFile;

public:

    Mdouble radiusParticle;
};



int main(){

    Compression script;

    script.dataFile.setFileType(FileType::ONE_FILE);
    script.restartFile.setFileType(FileType::ONE_FILE);
    script.fStatFile.setFileType(FileType::NO_FILE);
    script.eneFile.setFileType(FileType::NO_FILE);
    logger(INFO, "run number: %", script.dataFile.getCounter());



    Mdouble boxSize = 7e-4;
    script.setXMin(-0.5*boxSize);
    script.setYMin(-0.5*boxSize);
    script.setZMin(-0.5*boxSize);
    script.setXMax( 0.5*boxSize);
    script.setYMax( 0.5*boxSize);
    script.setZMax( 0.5*boxSize);
    script.solve();

    return 0;
}