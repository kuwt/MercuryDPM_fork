//
// Created by paolo on 4-10-19.
//
//Ciao!


#include <BaseCluster.h>
#include "Boundaries/StressStrainControlBoundary.h"
#include "Walls/InfiniteWall.h"
#include <ctime>

class script : public Mercury3D {

    void setupInitialConditions() override {

        //radiusParticle = 5e-4;
        radiusParticle = 0.05;
        restitutionCoefficient = 0.5;
        penetrationDepthMax = 0.1;
        densityParticle = 2000.0;
        slidingFrictionCoefficient = 0.0;
        loadingStiffnessIntra = 1000.0;
        unLoadingStiffnessMaxIntra = 5000.0;
        cohesionStiffness = 100000;
        sizeDispersityParticleIntra = 1;
        relativeTangentialStiffness = 0.0;



        //! Standard Output
        dataFile.setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        setXBallsAdditionalArguments("-solidf -v0");

        setSpecies();

        totalParticleVolume = 0;

        Mdouble tstart, tend;
        tstart = time(0);


        //! Creation of the cluster
        BaseCluster cluster;
        cluster.random.randomise();
        cluster.setCollisionTimeOverTimeStep(50);
        cluster.setVelocityDampingModulus(0.9);
        cluster.setInternalStructureGridLength(300);
        cluster.setEnergyRatioTolerance(1e-8);
        cluster.doCdatOutput(true);
        cluster.doOverlOutput(true);
        cluster.doAmatOutput(true);
        cluster.doIntStrucOutput(true);
        cluster.doVtkOutput(true);
        cluster.doRestartOutput(true);
        cluster.doFStatOutput(true);
        cluster.doEneOutput(true);
        //cluster.setNumberOfParticles(44);
        cluster.setRadiusCluster(4*radiusParticle);
        cluster.setPosition({0,0,0});
        cluster.setClusterId(0);
        cluster.setSizeDispersityParticle(sizeDispersityParticleIntra);
        cluster.setRadiusParticle(radiusParticle);
        cluster.setVelocity(Vec3D(1, 0, 0));//radiusParticle/(2e1*species->getCollisionTime(1)) , 0 , 0 ) );
        cluster.setParticleSpecies(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)) );
        cluster.solve();


        tend = time(0);
        std::cout << "It took "<< std::scientific << difftime(tend, tstart) <<" second(s)."<< std::endl;

        importParticlesAs(cluster.particleHandler,
                cluster.interactionHandler,
                dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)));




        std::cout << "Final mass fraction: " << cluster.getFinalMassFraction() << std::endl;


        //!total particle volume
        Mdouble rCluster = cluster.getMeanClusterRadius();
        Mdouble rForMassFraction = rCluster - 0.5 * radiusParticle;
        totalParticleVolume += constants::pi * 4 * pow(rForMassFraction, 3) *
                               (0.58 + 3 * pow(0.58, 2) * species->getPenetrationDepthMax()) / 3;

        //! Setting timeStep and timeMax
        setTimeStep(species->getCollisionTime(1) / 50); // (collision time)/50.0
        setTimeMax(species->getCollisionTime(1));
        fileOutputCount = 50;
        setSaveCount(10*fileOutputCount);
        fileOutputTime = fileOutputCount * getTimeStep();

        //! Setting Domain
        boxSize = 2*cluster.getMeanClusterRadius();
        setDomain({-boxSize, -boxSize, -boxSize}, {boxSize, boxSize, boxSize});

        LinearPlasticViscoelasticFrictionSpecies* speciesWall;
        speciesWall = species;
        speciesHandler.copyAndAddObject(speciesWall);

        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setLoadingStiffness(10*species->getLoadingStiffness());
        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setUnloadingStiffnessMax(10*species->getLoadingStiffness());
        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setCohesionStiffness(0);
        speciesHandler.getMixedObject(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)), dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)) )->setSlidingFrictionCoefficient(0);

        auto w0 = new InfiniteWall;
        w0->setSpecies(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)));
        w0->set({1,0,0}, {getXMax(), 0, 0});
        wallHandler.copyAndAddObject(w0);
        w0->set({-1,0,0}, {getXMin(), 0, 0});
        wallHandler.copyAndAddObject(w0);
        w0->set({0,1,0}, {0, getYMax(), 0});
        wallHandler.copyAndAddObject(w0);
        w0->set({0,-1,0}, {0, getYMin(), 0});
        wallHandler.copyAndAddObject(w0);
        w0->set({0,0,1}, {0, 0, getZMax()});
        wallHandler.copyAndAddObject(w0);
        w0->set({0,0,-1}, {0, 0, getZMin()});
        wallHandler.copyAndAddObject(w0);

        //!Simulation variables
        t0 = getTime();
        stage = 1;

        //! Printing some values of interest
        std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
        std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;

        setName("script");


        makeOverlFile();
        writeOverlFile();
        makeContactModelFile();
    }

    void actionsAfterTimeStep() override {

        if (fmod(getTime()-t0, fileOutputTime)<getTimeStep()) {
            makeDataAnalysis();
        }

        //makeOverlFile();
        //writeOverlFile();
        //setTimeMax(getTime());
    }

    void printTime() const override {

        std::cout <<

                  "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
                  ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
                  ", Stress: " << std::fixed << std::setprecision(2) << std::setw(11) << totalStress <<
                  ", uniStress: " << std::fixed << std::setprecision(2) << std::setw(11) << uniStress <<
                  ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8)
                  << getKineticEnergy() / getElasticEnergy() <<
                  ", Inertia = " << std::scientific << std::setprecision(2) << std::setw(8) << inertia <<
                  ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber_ <<
                  ", n = " << std::fixed << std::setprecision(3) << voidRatio <<
                  ", e = " << std::fixed << std::setprecision(3) << e <<
                  ", em = " << std::fixed << std::setprecision(3) << em <<
                  ", eM = " << std::fixed << std::setprecision(3) << eM <<
                  ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlapIntra_ <<
                  ", dMean = " << std::fixed << std::setw(7) << meanRelativeOverlapIntra_ <<
                  ", dMax = " << maxRelativeOverlapIntra_ <<
                  ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlapInter_ <<
                  ", dMean = " << std::fixed << std::setw(7) << meanRelativeOverlapInter_ <<
                  ", dMax = " << maxRelativeOverlapInter_

                  << std::endl;
    }

    void setSpecies() {


        speciesHandler.clear();


        Mdouble constantMass = densityParticle * 4 * constants::pi * pow(radiusParticle, 3) / 3;
        Mdouble collisionTimeIntra = sqrt(constantMass * (pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2)) /
                                          (2 * loadingStiffnessIntra));

        //! Setting Species
        species = new LinearPlasticViscoelasticFrictionSpecies;
        species->setConstantRestitution(true);
        species->setDensity(densityParticle);
        species->setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, restitutionCoefficient, 1);
        species->setUnloadingStiffnessMax(
                species->getLoadingStiffness() * unLoadingStiffnessMaxIntra / loadingStiffnessIntra);
        species->setCohesionStiffness(species->getUnloadingStiffnessMax() * cohesionStiffness / unLoadingStiffnessMaxIntra);
        species->setPenetrationDepthMax(penetrationDepthMax);

        species->setSlidingFrictionCoefficient(slidingFrictionCoefficient);
        species->setSlidingStiffness(species->getLoadingStiffness() * relativeTangentialStiffness);
        species->setSlidingDissipation(species->getDissipation() * 2.0 / 7.0);
        speciesHandler.copyAndAddObject(species);

    }


//! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
    void makeDataAnalysis() {
        //! Computes coordination number and overlaps.
        computeCoordinationNumberAndOverlaps();

        //! Computes void ratio (and volumeBox)
        computeVoidRatioAndPorosity();

        //! Computes stress and inertia.
        computeStress();
    }

//! Computes coordination number and overlaps
    void computeCoordinationNumberAndOverlaps() {
        meanCoordinationNumber_ = 0.0;
        maxRelativeOverlapIntra_ = 0.0;
        meanRelativeOverlapIntra_ = 0.0;
        minRelativeOverlapIntra_ = inf;
        maxRelativeOverlapInter_ = 0.0;
        meanRelativeOverlapInter_ = 0.0;
        minRelativeOverlapInter_ = inf;
        Mdouble relativeOverlap;
        int counterOverlapIntra = 0;
        int counterOverlapInter = 0;

        //! loops over each particle to compute mean coordination number.
        for (int i = 0; i < particleHandler.getNumberOfRealObjects(); i++) {
            meanCoordinationNumber_ += (particleHandler.getObject(i)->getInteractions()).size();
        }
        meanCoordinationNumber_ /= particleHandler.getSize();

        //! loops over every interaction to compute maximum, mean and minimum relative particle overlap.
        for (std::vector<BaseInteraction *>::const_iterator i = interactionHandler.begin();
             i != interactionHandler.end(); ++i) {
            /*!
             * \details the relative overlap is computed as an average of the relative overlap on the two particles.
             *          rO = ( O/R1 + O/R2 ) / 2.
             */
            relativeOverlap =
                    ((*i)->getOverlap()) / (particleHandler.getObject((*i)->getP()->getIndex())->getRadius()) +
                    ((*i)->getOverlap()) / (particleHandler.getObject((*i)->getI()->getIndex())->getRadius());
            relativeOverlap /= 2;
            if ((*i)->getP()->getGroupId() == (*i)->getI()->getGroupId()) {
                meanRelativeOverlapIntra_ += relativeOverlap;
                if (relativeOverlap > maxRelativeOverlapIntra_)
                    maxRelativeOverlapIntra_ = relativeOverlap;
                if (relativeOverlap < minRelativeOverlapIntra_)
                    minRelativeOverlapIntra_ = relativeOverlap;
                counterOverlapIntra++;
            } else {
                meanRelativeOverlapInter_ += relativeOverlap;
                if (relativeOverlap > maxRelativeOverlapInter_)
                    maxRelativeOverlapInter_ = relativeOverlap;
                if (relativeOverlap < minRelativeOverlapInter_)
                    minRelativeOverlapInter_ = relativeOverlap;
                counterOverlapInter++;
            }

        }
        meanRelativeOverlapIntra_ /= counterOverlapIntra;
        meanRelativeOverlapInter_ /= counterOverlapInter;
    }

//! computes voidRatio and porosity (e)
    void computeVoidRatioAndPorosity() {
        volumeBox = (getXMax() - getXMin()) * (getYMax() - getYMin()) * (getZMax() - getZMin());
        voidRatio = 1 - totalParticleVolume / volumeBox;
        e = (volumeBox - totalParticleVolume) / totalParticleVolume;
        Mdouble mFIntra = 0.58 + 3 * pow(0.58, 2) * meanRelativeOverlapIntra_;
        em = (1 - mFIntra) / mFIntra;
        eM = (e - em) / (1 + em);
    }

//! computes stress
    void computeStress() {
        stressXX = 0.0;
        stressYY = 0.0;
        stressZZ = 0.0;
        stressXY = 0.0;
        stressXZ = 0.0;
        stressYZ = 0.0;

        for (std::vector<BaseInteraction *>::const_iterator i = interactionHandler.begin();
             i != interactionHandler.end(); ++i) {
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

        totalStress = (stressXX + stressYY + stressZZ) / 3;
        uniStress = stressZZ;

        inertia = -2 * radiusParticle * isoStrainDot / sqrt(totalStress / densityParticle);

    }

    void makeOverlFile() {
        std::ostringstream overlName;
        overlName << getName() << ".overl";

        overlFile.open(overlName.str(), std::ios::out);

        overlFile << "OVERLAP VS DISTANCE" << std::endl << std::endl;

    }

    void writeOverlFile() {

        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {
            Mdouble force = (*i)->getForce().X*(*i)->getNormal().X + (*i)->getForce().Y*(*i)->getNormal().Y + (*i)->getForce().Z*(*i)->getNormal().Z;
            Mdouble distance = (*i)->getContactPoint().getLength()/radiusParticle;
            if ( (*i)->getContactPoint().X > -radiusParticle && (*i)->getContactPoint().X < radiusParticle )
                overlFile << (*i)->getOverlap()/radiusParticle
                          << std::setw(18) << std::scientific << std::setprecision(5) << distance
                          << std::setw(18) << force <<std::endl;
        }
        overlFile.close();
    }

    void makeContactModelFile(){

        Mdouble deltaStar = penetrationDepthMax * unLoadingStiffnessMaxIntra / (unLoadingStiffnessMaxIntra - loadingStiffnessIntra);

        std::ostringstream cmName;
        cmName << getName() << ".cm";

        cmFile.open(cmName.str(), std::ios::out);

        cmFile << 0 << "    " << 0 << std::endl;
        cmFile << penetrationDepthMax << "    " << 0 << std::endl;
        cmFile << penetrationDepthMax*1.5 << "    " << 0.5*penetrationDepthMax*radiusParticle*unLoadingStiffnessMaxIntra << std::endl;
        cmFile << deltaStar << "    " << deltaStar*radiusParticle*loadingStiffnessIntra << std::endl;
        cmFile << 0 << "    " << 0 << std::endl;

    }






    //! Geometry
    Mdouble radiusParticle;
    Mdouble boxSize;
    Mdouble isoStrainDot;

    //! Species and cluister
    LinearPlasticViscoelasticFrictionSpecies* species;
    int nClusters;
    Mdouble restitutionCoefficient;
    Mdouble penetrationDepthMax;
    Mdouble densityParticle;
    Mdouble slidingFrictionCoefficient;
    Mdouble loadingStiffnessIntra;
    Mdouble unLoadingStiffnessMaxIntra;
    Mdouble cohesionStiffness;
    Mdouble relativeTangentialStiffness;

    int nParticles;
    Mdouble sizeDispersityParticleIntra;

    //! Name and output
    std::ostringstream name;
    std::ofstream isoCFile;
    std::ofstream overlFile;
    std::ofstream cmFile;
    int fileOutputCount;
    Mdouble fileOutputTime;

    //! Data of interest
    Mdouble relativeOverlap_;
    Mdouble minRelativeOverlapIntra_;
    Mdouble meanRelativeOverlapIntra_;
    Mdouble maxRelativeOverlapIntra_;
    Mdouble minRelativeOverlapInter_;
    Mdouble meanRelativeOverlapInter_;
    Mdouble maxRelativeOverlapInter_;
    Mdouble meanCoordinationNumber_;
    Mdouble volumeBox;
    Mdouble voidRatio;
    Mdouble e;
    Mdouble em;
    Mdouble eM;
    Mdouble totalParticleVolume;
    Mdouble inertia;

    //! Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble totalStress;
    Mdouble uniStress;


    //! Computation
    Mdouble t0;
    int stage;
    Mdouble dampingCoefficient;
    std::vector<Mdouble> results;

    std::vector<Mdouble> eValues = {
            0.96,
            0.95919,
            0.95837,
            0.95755,
            0.95672,
            0.95588,
            0.95504,
            0.95419,
            0.95334,
            0.95248,
            0.95161,
            0.95074,
            0.94986,
            0.94898,
            0.94809,
            0.94719,
            0.94629,
            0.94538,
            0.94446,
            0.94354,
            0.94261,
            0.94168,
            0.94074,
            0.93979,
            0.93884,
            0.93788,
            0.93692,
            0.93595,
            0.93497,
            0.93399,
            0.933,
            0.93201,
            0.93101,
            0.93,
            0.92899,
            0.92797,
            0.92694,
            0.92591,
            0.92487,
            0.92383,
            0.92278,
            0.92172,
            0.92066,
            0.91959,
            0.91852,
            0.91744,
            0.91635,
            0.91526,
            0.91416,
            0.91305,
            0.91194,
            0.91083,
            0.9097,
            0.90857,
            0.90744,
            0.9063,
            0.90515,
            0.90399,
            0.90283,
            0.90167,
            0.9005,
            0.89932,
            0.89813,
            0.89694,
            0.89575,
            0.89454,
            0.89333,
            0.89212,
            0.8909,
            0.88967,
            0.88844,
            0.8872,
            0.88595,
            0.8847,
            0.88344,
            0.88218,
            0.88091,
            0.87963,
            0.87835,
            0.87706,
            0.87576,
            0.87446,
            0.87316,
            0.87184,
            0.87052,
            0.8692,
            0.86787,
            0.86653,
            0.86519,
            0.86384,
            0.86248,
            0.86112,
            0.85975,
            0.85837,
            0.85699,
            0.85561,
            0.85421,
            0.85282,
            0.85141,
            0.85,
    };

};




int main(){
    script problem;
    problem.solve();
}

