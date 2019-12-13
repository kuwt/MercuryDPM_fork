//
// Created by paolo on 17/05/19.
//

#include <iostream>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h> // Actually this is now useless (ClusterGenerator.h would be enough) but as a temporary solution it is fine.
#include <BaseCluster.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Walls/InfiniteWall.h>

class Script : public Mercury3D{

public:

    Mdouble radiusParticle = 1e-6;

    LinearPlasticViscoelasticFrictionSpecies* species = new LinearPlasticViscoelasticFrictionSpecies;

    //Used in print time (not really needed for the computation)
    Mdouble meanForceOnInteraction_;
    Mdouble relativeOverlap_;
    Mdouble minRelativeOverlap_;
    Mdouble meanRelativeOverlap_;
    Mdouble maxRelativeOverlap_;

    int nParticles = 300;
    Mdouble sizeDispersityParticle = 1.2;
    Mdouble densityParticle = 2500;
    Mdouble smallestMass = densityParticle * 4 * constants::pi * pow(2 * radiusParticle * (1+sizeDispersityParticle), 3) / 3;
    Mdouble wallRelativePosition;
    Mdouble initialPositionWall;
    Mdouble recordPosition;
    bool firstTime = true;
    bool start = true;

    Mdouble forceL;
    Mdouble forceU;
    InfiniteWall* wall1;
    InfiniteWall* wall2;
    Mdouble previousForceValue;
    Mdouble kEq;
    int fileOutputCount = 10000;
    Mdouble fileOutputTime;

    std::ofstream isoCFile;
    std::ofstream kFile;

    void setupInitialConditions() override {

            BaseCluster cluster;
            //cluster.random.randomise();
            cluster.setCollisionTimeOverTimeStep(50);
            cluster.setVelocityDampingModulus(0.9);
            cluster.setInternalStructureGridLength(300);
            cluster.setEnergyRatioTolerance(1e-8);
            cluster.doCdatOutput(false);
            cluster.doOverlOutput(false);
            cluster.doAmatOutput(false);
            cluster.doIntStrucOutput(false);
            cluster.doVtkOutput(false);
            cluster.doRestartOutput(true);
            cluster.doFStatOutput(false);
            cluster.doEneOutput(false);
            cluster.setNumberOfParticles(nParticles);
            cluster.setPosition({0,0,0});
            cluster.setClusterId(particleHandler.getNextGroupId());
            cluster.setSizeDispersityParticle(sizeDispersityParticle);
            cluster.setRadiusParticle(radiusParticle);
            cluster.setVelocity( Vec3D( radiusParticle/(2e1*species->getCollisionTime(1)) , 0 , 0 ) );
            cluster.setParticleSpecies(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)) );
            cluster.solve();


            importParticlesAs(cluster.particleHandler,
                    cluster.interactionHandler,
                    dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)));



        std::cout << "Mean cluster radius: " << cluster.getMeanClusterRadius() << std::endl;
        std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
        std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;

        initialPositionWall = 1.1 * cluster.getMeanClusterRadius();

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        wall1 = wallHandler.copyAndAddObject(w0);
        wall1->set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, - initialPositionWall));
        wall2 = wallHandler.copyAndAddObject(w0);
        wall2->set(Vec3D(0.0, 0.0,  1.0), Vec3D(0, 0,   initialPositionWall));

        previousForceValue = 0;
        fileOutputTime = fileOutputCount*getTimeStep();

        makeIsoCFile();
    }

    void actionsBeforeTimeLoop() override {

        std::cout << "N inter before = " << interactionHandler.getSize() << std::endl;

    }

    void actionsAfterTimeStep() override {

        computeData();

        if ( fmod(getTime(),fileOutputTime) < getTimeStep() )
            writeIsoCFile();

        moveWalls();

    }


    void actionsAfterSolve() override {
        std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
        std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;
    }

    void printTime() const override {
        std::cout
                << std::fixed << std::setprecision(2)
                << "Simulation progress: " << std::setw(6)  << 100*wallRelativePosition/initialPositionWall << "%"
                << ", Compression progress: " << std::setw(6) << (firstTime ? 0 : 100 * (wallRelativePosition - recordPosition)/(initialPositionWall-recordPosition)) << "%"
                << std::scientific
                << ", eRT: " << getKineticEnergy()/getElasticEnergy()
                << std::fixed << std::setprecision(5)
                << ", dMin: " << minRelativeOverlap_
                << ", dMean: " << meanRelativeOverlap_
                << ", dMax: " << maxRelativeOverlap_
                << std::scientific <<std::setprecision(2)
                << ", mFOI: " << meanForceOnInteraction_
                << ", uWP: " << initialPositionWall - wallRelativePosition
                << ", uWF: " << std::setw(10)  << forceU
                << ", lWF: " << std::setw(10)  << forceL
                << ", kEq: " << kEq << std::endl;
    }

    void makeIsoCFile() {
        std::ostringstream isoCFileName;
        isoCFileName << getName() << ".isoc";

        isoCFile.open(isoCFileName.str(), std::ios::out);
        isoCFile << "ISOTROPIC COMPRESSION OF SINGLE CLUSTER INFORMATION" << std::endl << std::endl;

        isoCFile << "radiusParticle: " << std:: scientific << std::setprecision(2) << radiusParticle << std::endl;
        isoCFile << "sizeDispersityParticle: " << std::defaultfloat << sizeDispersityParticle << std::endl;
        isoCFile << "densityParticle: " << std:: scientific << species -> getDensity() << std::endl;
        isoCFile << "nParticles: " << std::defaultfloat << nParticles << std::endl;
        isoCFile << "slidingFrictionCoeff: " << species -> getSlidingFrictionCoefficient() << std::endl;
        isoCFile << "rollingFrictionCoeff: " << species -> getRollingFrictionCoefficient() << std::endl;
        isoCFile << "torsionFrictionCoeff: " << species -> getTorsionFrictionCoefficient() << std::endl;
        isoCFile << "loadingStiffness: " << std::scientific << species -> getLoadingStiffness() * smallestMass << std::endl;
        isoCFile << "unloadingStiffnessMax: " << species -> getUnloadingStiffnessMax() * smallestMass << std::endl;
        isoCFile << "cohesionStiffness: " << species -> getCohesionStiffness() * smallestMass << std::endl;
        isoCFile << "restitutionCoefficient: " << species -> getRestitutionCoefficient(smallestMass) << std::endl;
        isoCFile << "penetrationDepthMax: " << species -> getPenetrationDepthMax() << std::endl;
        isoCFile << "collisionTime: " << std::scientific << std::setprecision(3) << species -> getCollisionTime(smallestMass) << std::endl;

        std::ostringstream kName;
        kName << getName() << ".k";
        kFile.open(kName.str(), std::ios::out);
    }

    void moveWalls() {

        wallRelativePosition = getTime() * initialPositionWall / getTimeMax();

        //  Move the boundary in z direction to next time step
        wall1->set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, - initialPositionWall + wallRelativePosition));

        wall2->set(Vec3D(0.0, 0.0,  1.0), Vec3D(0, 0,   initialPositionWall - wallRelativePosition));

    }

    void computeData() {

        if (firstTime && fabs(wall1->getForce().Z) > 0 && fabs(wall2->getForce().Z) > 0){
            recordPosition = wallRelativePosition;
            firstTime = false;
        }
        Mdouble Dx;
        meanForceOnInteraction_ = 0;
        relativeOverlap_ = 0;
        minRelativeOverlap_ = 2;
        meanRelativeOverlap_ = 0;
        maxRelativeOverlap_ = 0;
        //! loops over every interaction to compute mean force acting on interaction, maximum, mean and minimum relative particle overlap.
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            meanForceOnInteraction_ += ((*i) -> getForce()).getLength();

            /*!
             * \details the relative overlap is computed as an average of the relative overlap on the two particles.
             *          rO = ( O/R1 + O/R2 ) / 2.
             */
            relativeOverlap_ = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                               ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
            relativeOverlap_ /= 2;
            meanRelativeOverlap_ += relativeOverlap_;
            if (relativeOverlap_ > maxRelativeOverlap_)
                maxRelativeOverlap_ = relativeOverlap_;

            if (relativeOverlap_ < minRelativeOverlap_)
                minRelativeOverlap_ = relativeOverlap_;
        }
        meanForceOnInteraction_/= interactionHandler.getSize();
        meanRelativeOverlap_ /= interactionHandler.getSize();

        //Delta x, Ricordando che i muri che si muovono sono 2.
        Dx = 2 * fileOutputTime * getZMax() / getTimeMax();
        //Delta F / Delta Z.
        forceL = wall1->getForce().Z;
        forceU = wall2->getForce().Z;
        //kEq = (forceU - forceL - previousForceValue) / Dx;
        //previousForceValue = forceU - forceL;
        kEq = firstTime ? 0.0 : 0.5 * (forceU-forceL)/(wallRelativePosition-recordPosition);
    }

    void writeIsoCFile() {

        isoCFile
                << std::fixed << std::setprecision(2)
                << "Simulation progress: " << std::setw(6) << 100*wallRelativePosition/initialPositionWall << " %"
                << ", Compression progress: " << std::setw(6) << (firstTime ? 0 : 100 * (wallRelativePosition - recordPosition)/(initialPositionWall-recordPosition)) << " %"
                << std::scientific
                << ", eRT: " << getKineticEnergy()/getElasticEnergy()
                << std::fixed << std::setprecision(5)
                << ", dMin: " << minRelativeOverlap_
                << ", dMean: " << meanRelativeOverlap_
                << ", dMax: " << maxRelativeOverlap_
                << std::scientific <<std::setprecision(2)
                << ", mFOI: " << meanForceOnInteraction_
                << ", uWP: " << getXMax() - wallRelativePosition
                << ", uWF: " << std::setw(10)  << forceU
                << ", lWF: " << std::setw(10)  << forceL
                << ", kEq: " << kEq << std::endl;

        kFile
                << kEq << std::endl;
    }

};


int main(){

    Script script;

    Mdouble restitutionCoefficient = 0.5;
    Mdouble penetrationDepthMax = 0.1;
    Mdouble slidingFrictionCoefficient = 0.5;
    Mdouble slidingRollingCoefficient = 0.3;
    Mdouble slidingTorsionCoefficient = 0.0;

    Mdouble loadingStiffness = 1.0e4;
    Mdouble unLoadingStiffnessMax = 5.0e4;
    Mdouble cohesionStiffness = 5e3;

    Mdouble collisionTimeSmallestMass = sqrt( script.smallestMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );

    script.species -> setConstantRestitution(true);
    script.species -> setDensity(script.densityParticle);
    script.species -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
    script.species -> setUnloadingStiffnessMax(script.species -> getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
    script.species -> setCohesionStiffness(script.species -> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
    script.species -> setPenetrationDepthMax(penetrationDepthMax);

    script.species -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    script.species -> setSlidingStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setSlidingDissipation(script.species -> getDissipation()*2.0/7.0);
    script.species -> setRollingFrictionCoefficient(slidingRollingCoefficient);
    script.species -> setRollingStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setRollingDissipation(script.species -> getDissipation()*2.0/7.0);
    script.species -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
    script.species -> setTorsionStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setTorsionDissipation(script.species -> getDissipation()*2.0/7.0);

    script.speciesHandler.copyAndAddObject(script.species);

//! [T1:problemSetup]
    std::ostringstream name;
    name << "FirstTryClassic" << script.nParticles;
    script.setName(name.str());
    script.setSystemDimensions(3);

    script.setGravity(Vec3D(0.0, 0.0, 0.0));
    script.setXMax(1e-5);
    script.setYMax(1e-5);
    script.setZMax(1e-5);
    script.setXMin(-1e-5);
    script.setYMin(-1e-5);
    script.setZMin(-1e-5);

//! [T1:solve]
    //script.setGravity({0,0,-9.81});
    script.setTimeStep(script.species -> getCollisionTime(1)/50); // (collision time)/50.0
    script.setTimeMax(script.getTimeStep()*5e6);

//! [T1:output]
    script.setSaveCount( script.fileOutputCount );
    script.dataFile.setFileType(FileType::ONE_FILE);
    script.restartFile.setFileType(FileType::ONE_FILE);
    script.fStatFile.setFileType(FileType::NO_FILE);
    script.eneFile.setFileType(FileType::NO_FILE);
    logger(INFO, "run number: %", script.dataFile.getCounter());
//! [T1:output]

//! [T1:visualOutput]
    script.setXBallsAdditionalArguments("-solidf -v0");
    script.setNumberOfDomains({2,2,2});
//! [T1:visualOutput]
    script.solve();

}
