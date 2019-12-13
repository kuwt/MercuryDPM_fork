//
// Created by paolo on 21-9-19.
//

#include "SingleClusterCompression.h"

SingleClusterCompression::SingleClusterCompression(){

    radiusParticle = 6.64e-6;
    restitutionCoefficient = 0.5;
    penetrationDepthMax = 0.1231;
    densityParticle = 2700;
    
    loadingStiffness = 3881.90;
    unLoadingStiffnessMax = 11117.57;
    cohesionStiffness = 13341.08;
    relativeTangentialStiffness = 0.419403;
    slidingFrictionCoefficient = 0.242656;
    
    /*
    loadingStiffness = 256.45;
    unLoadingStiffnessMax = 1032.95;
    cohesionStiffness = 34952.96;
    relativeTangentialStiffness = 1.34742;
    slidingFrictionCoefficient = 0.141247;
    */
    nParticles = 300;
    sizeDispersityParticle = 1;

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

SingleClusterCompression::SingleClusterCompression(int nP, Mdouble rP, Mdouble sDP, Mdouble dP, Mdouble rC, Mdouble pDM,
                                                   Mdouble sFC, Mdouble lS, Mdouble uSM, Mdouble cS, Mdouble rTS){

    nParticles = nP;
    radiusParticle = rP;
    sizeDispersityParticle = sDP;
    densityParticle = dP;
    restitutionCoefficient = rC;
    penetrationDepthMax = pDM;
    slidingFrictionCoefficient = sFC;
    loadingStiffness = lS;
    unLoadingStiffnessMax = uSM;
    cohesionStiffness = cS;
    relativeTangentialStiffness = rTS;

}

SingleClusterCompression::~SingleClusterCompression(){

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

void SingleClusterCompression::setupInitialConditions(){

    //! Standard Output
    dataFile.setFileType(FileType::ONE_FILE);
    restartFile.setFileType(FileType::ONE_FILE);
    fStatFile.setFileType(FileType::NO_FILE);
    eneFile.setFileType(FileType::NO_FILE);
    setXBallsAdditionalArguments("-solidf -v0");

    Mdouble constantMass = densityParticle * 4 * constants::pi * pow(radiusParticle, 3) / 3;
    collisionTimeSmallestMass = sqrt( constantMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );


    //! Setting Species
    species = new LinearPlasticViscoelasticFrictionSpecies;
    species -> setConstantRestitution(true);
    species -> setDensity(densityParticle);
    species -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
    species -> setUnloadingStiffnessMax(species -> getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
    species -> setCohesionStiffness(species -> getUnloadingStiffnessMax() * cohesionStiffness / unLoadingStiffnessMax);
    species -> setPenetrationDepthMax(penetrationDepthMax);

    species -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    species -> setSlidingStiffness(species -> getLoadingStiffness()*relativeTangentialStiffness);
    species -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
    speciesHandler.copyAndAddObject(species);

    //! Creation of the cluster

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
    cluster.setParticleSpecies(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)));
    cluster.solve();

    importParticlesAs(cluster.particleHandler, cluster.interactionHandler, dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)));

    //! Setting timeStep and timeMax
    setTimeStep(species -> getCollisionTime(1)/50); // (collision time)/50.0
    setTimeMax(1000);
    fileOutputCount = 3000;
    setSaveCount(fileOutputCount);
    fileOutputTime = fileOutputCount*getTimeStep();

    //! Setting Domain
    boxSize = 1.2 * cluster.getMeanClusterRadius();
    setDomain({-boxSize,-boxSize,-boxSize},{boxSize,boxSize,boxSize});
    initialPositionWall = getZMax();

    //! Setting name
    setName("SingleClusterCompressionWorst");

    //! Printing some values of interest
    std::cout << "Mean cluster radius: " << cluster.getMeanClusterRadius() << std::endl;
    std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
    std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;

    //! Creating walls
    makeWalls();

    //! Creating IsoCFile
    makeIsoCFile();

    //! Setting starting booleans
    firstTime = true;

    t0 = 0;

}

void SingleClusterCompression::actionsAfterTimeStep(){

    computeData();

    if ( fmod(getTime(),fileOutputTime) < getTimeStep() ) {
        writeIsoCFile();
        record.push_back(kEq);
    }

    moveWalls();

    if ( (wallRelativePosition - recordPosition)/(initialPositionWall-recordPosition)>0.5 && !firstTime && getTime() - t0 > 20*fileOutputTime) setTimeMax(getTime());

}


void SingleClusterCompression::actionsAfterSolve(){

    Mdouble averageKeq = 0;
    for (int i = 0; i < record.size() ; ++i) {
        averageKeq += record[i];
    }
    averageKeq /= record.size();


    std::cout << "Average Keq =  " << averageKeq << std::endl;
    std::cout << "Keq/K2 =  " << averageKeq/unLoadingStiffnessMax << std::endl;
    std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
    std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;

    isoCFile << "Average Keq =  " << std::scientific << std::setprecision(4) << averageKeq << std::endl;
    isoCFile << "Keq/K2 =  " << std::scientific << std::setprecision(4) << averageKeq/unLoadingStiffnessMax << std::endl;

}

void SingleClusterCompression::printTime() const{
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

void SingleClusterCompression::makeWalls(){

    InfiniteWall w0;
    w0.setSpecies(speciesHandler.getObject(0));
    wall1 = wallHandler.copyAndAddObject(w0);
    wall1->set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, - initialPositionWall));
    wall2 = wallHandler.copyAndAddObject(w0);
    wall2->set(Vec3D(0.0, 0.0,  1.0), Vec3D(0, 0,   initialPositionWall));
}

void SingleClusterCompression::makeIsoCFile() {
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

}

void SingleClusterCompression::moveWalls() {

    setZMax(getZMax()*(1 - 1e-6));
    setZMin(getZMin()*(1 - 1e-6));
    wallRelativePosition = initialPositionWall-getZMax();

    //  Move the boundary in z direction to next time step
    wall1->set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));

    wall2->set(Vec3D(0.0, 0.0,  1.0), Vec3D(0, 0, getZMax()));

}

void SingleClusterCompression::computeData() {

    if (firstTime && fabs(wall1->getForce().Z) > 0 && fabs(wall2->getForce().Z) > 0){
        recordPosition = wallRelativePosition;
        firstTime = false;
        t0 = getTime();
    }
    //Mdouble Dx;
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
    //Dx = 2 * fileOutputTime * getZMax() / getTimeMax();
    //Delta F / Delta Z.
    forceL = wall1->getForce().Z;
    forceU = wall2->getForce().Z;
    //kEq = (forceU - forceL - previousForceValue) / Dx;
    //previousForceValue = forceU - forceL;
    kEq = firstTime ? 0.0 : 0.5 * (forceU-forceL)/( 2 * (wallRelativePosition-recordPosition) );

}

void SingleClusterCompression::writeIsoCFile() {

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
}
