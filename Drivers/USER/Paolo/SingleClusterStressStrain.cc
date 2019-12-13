//
// Created by paolo on 1-10-19.
//
//Ciao!

#include "SingleClusterStressStrain.h"

SingleClusterStressStrain::SingleClusterStressStrain(){

    radiusParticle = 6.64e-6;
    restitutionCoefficient = 0.5;
    penetrationDepthMax = 0.1231;
    densityParticle = 2700;
    slidingFrictionCoefficient = 0.25;
    loadingStiffness = 1e3;
    unLoadingStiffnessMax = 1.00000e+04;
    cohesionStiffness = unLoadingStiffnessMax/3;
    nParticles = 300;
    sizeDispersityParticle = 1;
    relativeTangentialStiffness = 0.25;

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

SingleClusterStressStrain::SingleClusterStressStrain(int nP, Mdouble rP, Mdouble sDP, Mdouble dP, Mdouble rC, Mdouble pDM,
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

SingleClusterStressStrain::~SingleClusterStressStrain(){

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

void SingleClusterStressStrain::setupInitialConditions(){

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
    species -> setCohesionStiffness(0);
    species -> setPenetrationDepthMax(penetrationDepthMax);
    species -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    species -> setSlidingStiffness(species -> getUnloadingStiffnessMax()*1e6*relativeTangentialStiffness);
    species -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
    speciesHandler.copyAndAddObject(species);

    speciesSimulation = new LinearPlasticViscoelasticFrictionSpecies;
    speciesSimulation = species;
    speciesSimulation -> setCohesionStiffness(speciesSimulation->getUnloadingStiffnessMax() * cohesionStiffness / unLoadingStiffnessMax);
    speciesHandler.copyAndAddObject(speciesSimulation);

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
    cluster.setSizeDispersityParticle(1);
    cluster.setRadiusParticle(radiusParticle);
    cluster.setParticleSpecies(dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(0)));
    cluster.solve();


    importParticlesAs(cluster.particleHandler, cluster.interactionHandler, dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(speciesHandler.getObject(1)));
    //! Setting timeStep and timeMax
    setTimeStep(species -> getCollisionTime(1)/50); // (collision time)/50.0
    setTimeMax(1000);
    fileOutputCount = 50;
    setSaveCount( 10*fileOutputCount );
    fileOutputTime = fileOutputCount*getTimeStep();

    //! Setting Domain
    rCluster = cluster.getMeanClusterRadius();
    boxSize = 1.1 * rCluster;
    setDomain({-boxSize,-boxSize,-boxSize},{boxSize,boxSize,boxSize});
    rForMassFraction = rCluster - 0.5 * radiusParticle;
    totalParticleVolume = constants::pi*4*pow(rForMassFraction,3)*(0.58+3*pow(0.58,2)*speciesSimulation->getPenetrationDepthMax())/3;

    //!Simulation variables
    relaxationTime = 10 * fileOutputTime;
    compressionTime = relaxationTime / 5;
    t0 = getTime();
    stage = 1;
    dampingCoefficient = 5e-5;
    relax = false;
    isStrainRateControlled = true;
    strainRatePerIteration = -1e-6;
    counterTotalForces = 1;
    counterPositiveForce = 0;


    //! Setting name
    name << "SingleClusterNewCM";
    setName(name.str());

    //! Printing some values of interest
    std::cout << "Mean cluster radius: " << rCluster << std::endl;
    std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
    std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;

    //! Creating walls
    makeWalls();

    //! Creating IsoCFile
    makeIsoCFile();
}

void SingleClusterStressStrain::actionsAfterTimeStep(){

    makeDataAnalysis();

    if (getNumberOfTimeSteps()==1) {

        printTime();
        auto i = interactionHandler.getObject(1);

        std::cout << std::scientific << std::setprecision(5)
        << (*i).getForce().X*(*i).getNormal().X + (*i).getForce().Y*(*i).getNormal().Y + (*i).getForce().Z*(*i).getNormal().Z << std::endl
        << speciesSimulation->getPenetrationDepthMax() - i->getOverlap()/radiusParticle
        << std::endl;

    }

    dampVelocities();

    if ( fmod(getTime(),fileOutputTime) < getTimeStep() ) writeIsoCFile();

    if (e < 0.85) setTimeMax(getTime());

    //increaseCohesiveForces();

}


void SingleClusterStressStrain::actionsAfterSolve(){
    std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
    std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;
}

void SingleClusterStressStrain::printTime() const{

    std::cout <<

              "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
              ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
              ", Stress: " << std:: fixed << std::setprecision(2) << std::setw(11) << totalStress <<
              ", uniStress: " << std:: fixed << std::setprecision(2) << std::setw(11) << uniStress <<
              ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
              ", Inertia = " << std::scientific << std::setprecision(2) << std::setw(8) << inertia <<
              ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber_ <<
              ", n = " << std::fixed << std::setprecision(3) << voidRatio <<
              ", e = " << std::fixed << std::setprecision(3) << e <<
              ", em = " << std::fixed << std::setprecision(3) << em <<
              ", eM = " << std::fixed << std::setprecision(3) << eM <<
              ", floaters = " << std::fixed << std::setw(3) << floaters <<
              ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlap_ <<
              ", dMean = " << std:: fixed << std::setw(7) << meanRelativeOverlap_ <<
              ", dMax = " << maxRelativeOverlap_ <<
              ", PercPos = " << std::setw(5) << 100*counterPositiveForce/counterTotalForces << "%"

              << std::endl;
}

//void SingleClusterStressStrain::computeExternalForces(BaseParticle *CI) {
//    DPMBase::computeExternalForces(CI);
//    for(auto i = CI->getInteractions().begin(); i != CI->getInteractions().end(); ++i){
//        if(Vec3D::dot( (*i)->getForce(), (*i)->getNormal() ) < 0 && (*i)->getP()->getGroupId() == (*i)->getI()->getGroupId() )
//            (*i)->getP()->addForce(0.25*(*i)->getForce());
//    }
//}

void SingleClusterStressStrain::makeWalls(){

    initialPositionWall = 1.1 * rCluster;



    isoStrainDot = strainRatePerIteration/getTimeStep();
    boundaryHandler.clear();
    w0.setHandler(&boundaryHandler);
    w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(isoStrainDot, 0, 0, 0, isoStrainDot, 0, 0, 0, isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), isStrainRateControlled);
    //w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), Matrix3D(0,0,0,0,0,0,0,0,0), false);
    boundaryHandler.copyAndAddObject(w0);
}

void SingleClusterStressStrain::makeIsoCFile() {
    std::ostringstream cdatName;
    cdatName << getName() << ".isoc";

    isoCFile.open(cdatName.str(), std::ios::out);

    isoCFile << "CLUSTER DATA AND INFORMATION" << std::endl << std::endl;
    isoCFile << "radiusParticle: " << std:: scientific << std::setprecision(2) << radiusParticle << std::endl;
    isoCFile << "sizeDispersityParticle: " << std::defaultfloat << sizeDispersityParticle << std::endl;
    isoCFile << "nParticles: " << nParticles << std::endl;
    isoCFile << "boxSize: " << std::defaultfloat << (getXMax()-getXMin()) << std::endl;
    isoCFile << "Damping coefficient: " << std::defaultfloat << dampingCoefficient << std::endl;
    isoCFile << "Strain Rate Controlled: " << isStrainRateControlled << std::endl;
    isoCFile << "relax: " << relax << std::endl;
    isoCFile << "strainRatePerIteration: " << strainRatePerIteration << std::endl;
    isoCFile << "Loading stiffness: " << loadingStiffness << std::endl;
    isoCFile << "Unloading stiffness max: " << unLoadingStiffnessMax << std::endl;
    isoCFile << "Cohesion stiffness: " << speciesSimulation->getCohesionStiffness()*smallestMass / 3 << std::endl;


    isoCFile << "COMPRESSION INFORMATION" << std::endl << std::endl;
    isoCFile <<

             "   Time" << std::setw(13) <<
             "Time Max" << std::setw(16) <<
             "Stress" << std::setw(16) <<
             "uniStress" << std::setw(15) <<
             "E_ratio" << std::setw(11) <<
             "Inertia" << std::setw(12) <<
             "cN" << std::setw(12) <<
             "n" << std::setw(12) <<
             "e" << std::setw(12) <<
             "em" << std::setw(12) <<
             "eM" << std::setw(10) <<
             "floaters" << std::setw(14) <<
             "dMin" << std::setw(14) <<
             "dMean" << std::setw(14) <<
             "dMax"

             << std::endl;
}


void SingleClusterStressStrain::dampVelocities()
{
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p)
    {
        (*p)->addForce(-dampingCoefficient * (*p)->getVelocity() * pow((*p)->getRadius()/radiusParticle,3) );
    }
}


//! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
void SingleClusterStressStrain::makeDataAnalysis()
{
    //! Computes coordination number and overlaps.
    computeCoordinationNumberAndOverlaps();

    //! Computes void ratio (and volumeBox)
    computeVoidRatioAndPorosity();

    //! Computes stress and inertia.
    computeStress();
}

//! Computes coordination number and overlaps
void SingleClusterStressStrain::computeCoordinationNumberAndOverlaps()
{
    meanForceOnInteraction_ = 0.0;
    meanCoordinationNumber_ = 0.0;
    maxRelativeOverlap_ = 0.0;
    meanRelativeOverlap_ = 0.0;
    minRelativeOverlap_ = inf;
    Mdouble relativeOverlap;
    floaters = 0;

    //! loops over each particle to compute mean coordination number and center of mass.
    for (int i=0; i < particleHandler.getNumberOfRealObjects(); i++)
    {
        meanCoordinationNumber_ += (particleHandler.getObject(i)->getInteractions()).size();
        if (particleHandler.getObject(i)->getInteractions().empty())
            floaters ++;
    }
    meanCoordinationNumber_ /= particleHandler.getNumberOfRealObjects();

    //! loops over every interaction to compute mean force acting on interaction, maximum, mean and minimum relative particle overlap.
    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
    {
        meanForceOnInteraction_ += ((*i) -> getForce()).getLength();

        /*!
         * \details the relative overlap is computed as an average of the relative overlap on the two particles.
         *          rO = ( O/R1 + O/R2 ) / 2.
         */
        relativeOverlap = ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                          ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
        relativeOverlap /= 2;
        meanRelativeOverlap_ += relativeOverlap;
        if (relativeOverlap > maxRelativeOverlap_)
            maxRelativeOverlap_ = relativeOverlap;

        if (relativeOverlap < minRelativeOverlap_)
            minRelativeOverlap_ = relativeOverlap;
    }
    meanForceOnInteraction_/= interactionHandler.getSize();
    meanRelativeOverlap_ /= interactionHandler.getSize();
}

//! computes voidRatio and porosity (e)
void SingleClusterStressStrain::computeVoidRatioAndPorosity()
{
    volumeBox = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
    voidRatio = 1 - totalParticleVolume / volumeBox;
    e = (volumeBox - totalParticleVolume) / totalParticleVolume;
    Mdouble mF = 0.58+3*pow(0.58,2)*meanRelativeOverlap_;
    em = (1-mF)/mF;
    eM = (e-em)/(1+em);
}

//! computes stress
void SingleClusterStressStrain::computeStress()
{
    stressXX = 0.0;
    stressYY = 0.0;
    stressZZ = 0.0;
    stressXY = 0.0;
    stressXZ = 0.0;
    stressYZ = 0.0;

    counterTotalForces = 0;
    counterPositiveForce = 0;

    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {
        stressXX += (*i)->getForce().X * (*i)->getNormal().X * (*i)->getDistance();
        stressYY += (*i)->getForce().Y * (*i)->getNormal().Y * (*i)->getDistance();
        stressZZ += (*i)->getForce().Z * (*i)->getNormal().Z * (*i)->getDistance();
        counterTotalForces ++;

        if ( Vec3D::dot( (*i)->getForce(), (*i)->getNormal()) < 0 ) counterPositiveForce++;
    }

    stressXX /= volumeBox;
    stressYY /= volumeBox;
    stressZZ /= volumeBox;
    stressXY /= volumeBox;
    stressXZ /= volumeBox;
    stressYZ /= volumeBox;

    totalStress = (stressXX+stressYY+stressZZ)/3;
    uniStress = stressZZ;

    inertia = -2*radiusParticle*isoStrainDot/sqrt(totalStress/densityParticle);

}

void SingleClusterStressStrain::increaseCohesiveForces(){

    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {
        Mdouble force;
        force = Vec3D::dot((*i)->getForce(),(*i)->getNormal());
        if (force < 0)
            (*i)->setForce(1e100*(*i)->getForce());
    }

}

void SingleClusterStressStrain::writeIsoCFile() {

    isoCFile <<

             std::setprecision(2) << std::setw(7) << getTime() <<
             std::setprecision(2) << std::setw(13) << getTimeMax() <<
             std::fixed << std::setprecision(2) << std::setw(16) << totalStress <<
             std::fixed << std::setprecision(2) << std::setw(16) << uniStress <<
             std::scientific << std::setprecision(2) << std::setw(15) << getKineticEnergy()/getElasticEnergy() <<
             std::scientific << std::setprecision(2) << std::setw(11) << inertia <<
             std::fixed << std::setw(12) << meanCoordinationNumber_ <<
             std::fixed << std::setprecision(3) << std::setw(12) << voidRatio <<
             std::fixed << std::setprecision(3) << std::setw(12) << e <<
             std::fixed << std::setprecision(3) << std::setw(12) << em <<
             std::fixed << std::setprecision(3) << std::setw(12) << eM <<
             std::fixed << std::setw(10) << floaters <<
             std::fixed << std::setprecision(5) << std::setw(14) << minRelativeOverlap_ <<
             std::setw(14) << meanRelativeOverlap_ <<
             std::setw(14) << maxRelativeOverlap_ <<
             std::setw(5) << 100*counterPositiveForce/counterTotalForces << "%"

             << std::endl;
}
