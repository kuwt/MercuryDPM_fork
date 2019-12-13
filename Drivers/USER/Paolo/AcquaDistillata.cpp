//
// Created by paolo on 17/05/19.
//

#include <iostream>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h> // Actually this is now useless (ClusterGenerator.h would be enough) but as a temporary solution it is fine.
#include <ClusterGenerator.h>
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
    int nClusters = 20;
    Mdouble sizeDispersityParticle = 0.1;
    Mdouble densityParticle = 2500;
    Mdouble smallestMass = densityParticle * 4 * constants::pi * pow(radiusParticle*(1-sizeDispersityParticle), 3) / 3;
    bool start = true;

    InfiniteWall* wall1;
    AxisymmetricIntersectionOfWalls* cylinder;
    int fileOutputCount = 10000;
    Mdouble fileOutputTime;

    std::ofstream isoCFile;

    void setupInitialConditions() override {

        Mdouble Theta;
        Mdouble X;
        Mdouble Y;
        Mdouble Z;

        for (int i = 0; i < nClusters; ++i) {

            Theta = constants::pi * ( floor(i/4) / 4 + fmod(i,4) /2 );

            X = 10 * radiusParticle * cos(Theta + 2 * constants::pi * fmod(i,4) );

            Y = 10 * radiusParticle * sin(Theta);

            Z = floor(i/4) * 16 * radiusParticle + 10 * radiusParticle;

            clusterGenerator = new ClusterGenerator;

            clusterGenerator->clusterProperties.random.randomise();

            clusterGenerator->clusterProperties.setPosition({X, Y, Z});

            clusterGenerator->clusterProperties.setNumberOfParticles(nParticles);

            clusterGenerator->clusterProperties.setSizeDispersityParticle(sizeDispersityParticle);

            clusterGenerator->clusterProperties.setRadiusParticle(radiusParticle);

            clusterGenerator->clusterProperties.setParticleSpecies(species);

            clusterGenerator->create();

            std::cout << "Creato Cluster Numero " << i << std::endl;

            particleHandler.copyContentsFromOtherHandler(clusterGenerator->clusterProperties.particleHandler);
            //interactionHandler.copyContentsFromOtherHandler(clusterGenerator->clusterProperties.interactionHandler);
        }

        std::cout << "Mean cluster radius: " << clusterGenerator->clusterProperties.getMeanClusterRadius() << std::endl;
        std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
        std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        wall1 = wallHandler.copyAndAddObject(w0);
        wall1->set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, 0));

        AxisymmetricIntersectionOfWalls c0;
        c0.setSpecies(species);
        c0.setPosition(Vec3D(0.0,0.0,0.0));
        c0.setOrientation(Vec3D(0.0,0.0,1.0));
        c0.addObject(Vec3D(1.0,0.0,0.0),Vec3D( 20 * radiusParticle , 0.0 , 0.0 ) );
        cylinder = wallHandler.copyAndAddObject(c0);



        setGravity({0,0,-1e9});

        fileOutputTime = fileOutputCount*getTimeStep();

        makeIsoCFile();
    }

    void actionsBeforeTimeLoop() override {
    }

    void actionsBeforeTimeStep() override {
        if(start){
            LinearPlasticViscoelasticInteraction* interaction;
            Mdouble r;
            Mdouble deltaStar = 0.125;
            Mdouble delta;
            BaseInteractable* P;
            BaseInteractable* I;
            for (int i = 0; i < interactionHandler.getSize(); ++i) {
                interaction =  dynamic_cast<LinearPlasticViscoelasticInteraction*>(interactionHandler.getObject(i));
                P=interactionHandler.getObject(i)->getP();
                I=interactionHandler.getObject(i)->getI();
                r = 2 * particleHandler.getObject(P->getIndex())->getRadius() * particleHandler.getObject(I->getIndex())->getRadius() / ( particleHandler.getObject(P->getIndex())->getRadius() + particleHandler.getObject(I->getIndex())->getRadius() );
                delta = interactionHandler.getObject(i)->getOverlap() / r;
                interaction->setMaxOverlap( 0.5*interaction->getMaxOverlap()*(1 + sqrt( 1 + deltaStar / delta ) ) );
            }
            start = false;
        }
    }

    void actionsAfterTimeStep() override {

        computeData();

        if ( fmod(getTime(),fileOutputTime) < getTimeStep() )
            writeIsoCFile();

        dampVelocities();



    }


    void actionsAfterSolve() override {
        std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
        std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;
    }

    void printTime() const override {
        std::cout
                << std::fixed << std::setprecision(2)
                << "Simulation progress: " << std::setw(6)  << 100 * getTime() / getTimeMax() << "%"
                << std::scientific
                << ", eRT: " << getKineticEnergy()/getElasticEnergy()
                << std::fixed << std::setprecision(5)
                << ", dMin: " << minRelativeOverlap_
                << ", dMean: " << meanRelativeOverlap_
                << ", dMax: " << maxRelativeOverlap_
                << std::scientific <<std::setprecision(2)
                << ", mFOI: " << meanForceOnInteraction_ << std::endl;
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
        isoCFile << "meanClusterRadius: " << std::scientific << std::setprecision(3) << clusterGenerator -> clusterProperties.getMeanClusterRadius() << std::endl << std::endl << std::endl;

    }

    void computeData() {

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
    }

    void writeIsoCFile() {

        isoCFile
                << std::fixed << std::setprecision(2)
                << "Simulation progress: " << std::setw(6) << 100 * getTime() / getTimeMax() << "%"
                << std::scientific
                << ", eRT: " << getKineticEnergy()/getElasticEnergy()
                << std::fixed << std::setprecision(5)
                << ", dMin: " << minRelativeOverlap_
                << ", dMean: " << meanRelativeOverlap_
                << ", dMax: " << maxRelativeOverlap_
                << std::scientific <<std::setprecision(2)
                << ", mFOI: " << meanForceOnInteraction_ << std::endl;
    }

    void dampVelocities()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            particleHandler.getObject(i) -> setVelocity(0.9*(particleHandler.getObject(i) -> getVelocity()));
        }
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
    Mdouble cohesionStiffness = 0.5e4;

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
    name << "Acqua_DistillataRifatto" << script.nParticles;
    script.setName(name.str());
    script.setSystemDimensions(3);

    script.setGravity(Vec3D(0.0, 0.0, 0.0));
    script.setXMax(2.1e-5);
    script.setYMax(2.1e-5);
    script.setZMax(1e-4);
    script.setXMin(-2.1e-5);
    script.setYMin(-2.1e-5);
    script.setZMin(0.0);

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
//! [T1:visualOutput]
    script.solve();

}
