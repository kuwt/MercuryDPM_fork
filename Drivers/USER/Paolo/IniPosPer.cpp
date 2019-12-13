//
// Created by paolo on 31-8-19.
//
//Ciao!

#include <iostream>
#include <Species/Species.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h> // Actually this is now useless (ClusterGenerator.h would be enough) but as a temporary solution it is fine.
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>

class Script : public Mercury3D{

public:

    void setupInitialConditions() override {

        random.randomise();

        setClusterPositions();

        InfiniteWall w0;
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, -a/2));
        w0.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(w0);


        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(a/2, 0, 0));
        w0.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(w0);


        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(-a/2, 0, 0));
        w0.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(w0);


        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, a/2, 0));
        w0.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(w0);


        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, -a/2, 0));
        w0.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(w0);


        setGravity({0,0,-1e8});

        isoStrainDot = -1e-7 / getTimeStep();

        velocityDM = 0.9;

        stage = 1;
    }

    void actionsAfterTimeStep() override {

        if (stage == 1) {

            computeData();

            dampVelocities();

            if (getKineticEnergy() / getElasticEnergy() < 1e-20) {

                setGravity({0, 0, 0});
                createBoundaries();
                velocityDM = 0.1;
                t = getTime();

                InfiniteWall w0;
                w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0, 0, a/2));
                w0.setSpecies(speciesHandler.getObject(0));
                wallHandler.copyAndAddObject(w0);

                std::cout << "Agito" << std::endl;

                for (int i = 0; i < particleHandler.getSize(); ++i) {
                    particleHandler.getObject(i)->setVelocity({0,0,0.04*random.getRandomNumber()});
                }

                fileOutputCount /= 1000;
                setSaveCount(fileOutputCount);

                stage = 2;

            }
        }

        if (stage == 2) {

            dampVelocities();

            computeData();

            if (getTime()-t > 100*fileOutputCount*getTimeStep()) {

                for (int i = 0; i < nClusters; ++i) {
                    std::cout << std::scientific << std::setprecision(5) << "{"
                              << particleHandler.getObject(i)->getPosition().X << ","
                              << particleHandler.getObject(i)->getPosition().Y << ","
                              << particleHandler.getObject(i)->getPosition().Z << "},"
                              << std::endl;
                }
                zu = particleHandler.getHighestPositionComponentParticle(2)->getPosition().Z + 5.65e-6;
                zl = particleHandler.getLowestPositionComponentParticle(2)->getPosition().Z - 5.65e-6;
                std::cout << std::endl << std::scientific << std::setprecision(5) << zu << std::endl << std::endl;
                std::cout << std::endl << std::scientific << std::setprecision(5) << zl << std::endl << std::endl;
                std::cout << wallPositionX << std::endl;
                std::cout << wallPositionY << std::endl;
                std::cout << wallPositionZ << std::endl;
                setTimeMax(getTime());
            }
        }

    }

    void actionsAfterSolve() override {
        std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
        std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;
    }

    void printTime() const override {
        std::cout
                << std::scientific << std::setprecision(2)
                << "Time: " << std::setw(6)  << getTime()
                << ", eR: " << getKineticEnergy()/getElasticEnergy()
                << ", vR: " << vR
                << ", vRPost: " << vRPost
                << std::fixed << std::setprecision(5)
                << ", dMin: " << minRelativeOverlap_
                << ", dMean: " << meanRelativeOverlap_
                << ", dMax: " << maxRelativeOverlap_
                << std::scientific <<std::setprecision(2)
                << ", mFOI: " << meanForceOnInteraction_ << std::endl;
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
        zu = particleHandler.getHighestPositionComponentParticle(2)->getPosition().Z + 5.65e-6;
        zl = particleHandler.getLowestPositionComponentParticle(2)->getPosition().Z - 5.65e-6;
        vR = 1-particleHandler.getVolume()/pow(a,3);
        vRPost = 1-(11 * 300 * 4 * constants::pi * 1e-18 / 3) / pow(a,3);
    }

    //! inserts the particles in the simulation domainss
    void setClusterPositions()
    {
        int nClustersInserted = 0;
        // insert particles till the desired amount is reached
        while (nClustersInserted < nClusters)
        {    // nParticle inserted corresponds to the index of the particle which is being inserted
            if (clusterInsertionSuccessful(nClustersInserted))
            {
                nClustersInserted++;
                std::cout << "Set position of cluster n. " << nClustersInserted << "/" << nClusters << std::endl;
            }
            else
            {
                logger(ERROR, "Cannot insert all clusters, try to decrease the value of initialSolidFraction.");
            }
        }
        std::cout << "CLUSTERS POSITION SETTING TERMINATED" << std::endl << std::endl;
    }

    //! Check if the particle is inserted whitout interactions
    bool clusterInsertionSuccessful(int n)
    {
        // initialization of parameters, n is the particle index
        int insertionFailCounter = 0;
        Vec3D clusterPosition;
        SphericalParticle c0;

        // setup of particle properties and initial conditions (besides position)
        c0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        c0.setRadius(radiusCluster);
        c0.setSpecies(speciesHandler.getObject(0));


        // in this cycle a random position inside of a sphere contained in the bounding box is taken
        // if the particle is not in contact with the others then insert it, otherwise increment the fail counter
        // the maximum number of failed attempts is capped to 1000

        while (insertionFailCounter < 1000)
        {

            clusterPosition.X = (getXMax() - getXMin() - 2*radiusCluster)*random.getRandomNumber(0.0,1.0) + getXMin() + radiusCluster;
            clusterPosition.Y = (getYMax() - getYMin() - 2*radiusCluster)*random.getRandomNumber(0.0,1.0) + getYMin() + radiusCluster;
            clusterPosition.Z = (getZMax() - getZMin() - 2*radiusCluster)*random.getRandomNumber(0.0,1.0) + getZMin() + radiusCluster;

            c0.setPosition(clusterPosition);

            if (checkParticleForInteraction(c0))
            {
                particleHandler.copyAndAddObject(c0);
                clusterPositions.push_back(clusterPosition);
                return true;
            }
            insertionFailCounter++;
        }
        return false;
    }

    void dampVelocities()
    {
        for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            particleHandler.getObject(i) -> setVelocity(velocityDM*(particleHandler.getObject(i) -> getVelocity()));
        }
    }

    void createBoundaries() {

        wallHandler.clear();

        wallPositionX = 3.95e-05/2;
        wallPositionY = 3.95e-05/2;
        wallPositionZ = 3.95e-05/2;

        wall.set(Vec3D(0, 0, 1), -wallPositionZ, wallPositionZ);
        WNS = boundaryHandler.copyAndAddObject(wall);

        wall.set(Vec3D(1, 0, 0), -wallPositionX, wallPositionX);
        WEW = boundaryHandler.copyAndAddObject(wall);

        wall.set(Vec3D(0, 1, 0), -wallPositionY, wallPositionY);
        WFB = boundaryHandler.copyAndAddObject(wall);

    }

    void moveBoundaries() {

        wallPositionX *= 1+isoStrainDot*getTimeStep();
        wallPositionY *= 1+isoStrainDot*getTimeStep();
        wallPositionZ *= 1+isoStrainDot*getTimeStep();

        WNS->set(Vec3D(0, 0, 1), -wallPositionZ, wallPositionZ);

        WEW->set(Vec3D(1, 0, 0), -wallPositionX, wallPositionX);

        WFB->set(Vec3D(0, 1, 0), -wallPositionY, wallPositionY);
    }



    int fileOutputCount;

    Mdouble a;
    Mdouble b;
    std::vector<Vec3D> clusterPositions;

    Mdouble wallPositionX;
    Mdouble wallPositionY;
    Mdouble wallPositionZ;
    PeriodicBoundary wall;
    PeriodicBoundary* WNS;
    PeriodicBoundary* WEW;
    PeriodicBoundary* WFB;
    Mdouble isoStrainDot;


    Mdouble radiusCluster = 5.65e-6;
    Mdouble velocityDM;
    Mdouble t;


    //Used in print time (not really needed for the computation)
    Mdouble meanForceOnInteraction_;
    Mdouble relativeOverlap_;
    Mdouble minRelativeOverlap_;
    Mdouble meanRelativeOverlap_;
    Mdouble maxRelativeOverlap_;
    Mdouble zu;
    Mdouble zl;
    Mdouble vR;
    Mdouble vRPost;

    int nClusters = 33;
    int stage;



};


int main(){

    Script script;

    LinearViscoelasticSpecies species;
    species.setDensity(2500.0); //sets the species type_0 density
    species.setStiffness(258000.5);//sets the spring stiffness.
    species.setDissipation(0.0); //sets the dissipation.
    script.speciesHandler.copyAndAddObject(species);

//! [T1:problemSetup]
    std::ostringstream name;
    name << "IniPosPer";
    script.setName(name.str());
    script.setSystemDimensions(3);

    //script.a = 3.50e-5;
    //script.a = 3.95e-5;
    script.a = 3.85e-5;
    //script.b = 6.42e-5;
    script.b = 1e-4;

    script.setXMax(script.a/2);
    script.setYMax(script.a/2);
    script.setZMax(script.b-script.a/2);
    script.setXMin(-script.a/2);
    script.setYMin(-script.a/2);
    script.setZMin(-script.a/2);

//! [T1:solve]
    //script.setGravity({0,0,-9.81});
    std::cout << "Il collision time è:" << species.getCollisionTime(5.36e-12) << std::endl;
    script.setTimeStep(species.getCollisionTime(5.36e-12) / 50.0); // (collision time)/50.0
    script.setTimeMax(1000);

    script.fileOutputCount = 50000;

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

