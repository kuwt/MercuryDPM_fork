//
// Created by paolo on 4-10-19.
//
//Ciao!

#include "DoublePorosityStressStrainC.h"

DoublePorosityStressStrainC::DoublePorosityStressStrainC(){

    radiusParticle = 6.875e-6;
    restitutionCoefficient = 0.5;
    penetrationDepthMax = 0.0832;
    densityParticle = 2700;
    slidingFrictionCoefficient = 0.25;
    loadingStiffnessIntra = 1e3;
    unLoadingStiffnessMaxIntra = 5e3;
    cohesionStiffness = 10*loadingStiffnessIntra;
    nParticles = 800;
    sizeDispersityParticleIntra = 1.2;
    relativeTangentialStiffness = 0.25;
    nClusters = 1000;
    //boxSize = ?;

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

DoublePorosityStressStrainC::DoublePorosityStressStrainC(std::vector<Mdouble> radC, std::vector<Vec3D> posC, Mdouble rappCompress, Mdouble rExpCSIntra, Mdouble sFC, Mdouble relKt, Mdouble bSX, Mdouble bSY, Mdouble bSZ, std::string n){
    radiiCluster = radC;
    positionsCluster = posC;
    //radiusParticle = 6.29e-6;
    sizeDispersityParticleIntra = 1; //No dispersity per avere risultati più precisi col modello di luca
    densityParticle = 2700;
    restitutionCoefficient = 0.5;
    penetrationDepthMax = 0.251;
    Mdouble Nl = pow(50/0.58, 1.0/3.0);
    Mdouble RHat = Nl * (1 - 0.58 * penetrationDepthMax);
    radiusParticle = 0.9 * 5.5e-5 / RHat; // 0.8 safety factor per non farli toccare all'inizio
    unLoadingStiffnessMaxIntra = pow(10,4); //! Fisso tanto mi serve solo il rapporto e poi riscalo per trovare i veri valori di stress.
    compressionRatio = rappCompress;
    loadingStiffnessIntra      = 0.2 * unLoadingStiffnessMaxIntra;
    cohesionStiffness          = pow(10,rExpCSIntra) * unLoadingStiffnessMaxIntra;
    loadingStiffnessInter      = 5 * unLoadingStiffnessMaxIntra;
    relativeTangentialStiffness = relKt;
    slidingFrictionCoefficient = sFC;
    nClusters = 300;
    setDomain({-bSX, -bSY, -bSZ}, {bSX, bSY, bSZ});
    setName(n);
}

DoublePorosityStressStrainC::~DoublePorosityStressStrainC(){

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

void DoublePorosityStressStrainC::setupInitialConditions(){

    //! Standard Output
    dataFile.setFileType(FileType::ONE_FILE);
    restartFile.setFileType(FileType::ONE_FILE);
    fStatFile.setFileType(FileType::NO_FILE);
    eneFile.setFileType(FileType::NO_FILE);
    setXBallsAdditionalArguments("-solidf -v0");

    //!File per gnuplot
    makeContactModelFile();

    setSpecies();

    totalParticleVolume = 0;

    for (int i = 0; i < nClusters; ++i) {

        //! Creation of the cluster
        BaseCluster cluster;
        //cluster.random.randomise();
        cluster.setCollisionTimeOverTimeStep(50);
        cluster.setVelocityDampingModulus(0.9);
        cluster.setInternalStructureGridLength(300);
        cluster.setEnergyRatioTolerance(2e-7);
        cluster.doCdatOutput(false);
        cluster.doOverlOutput(false);
        cluster.doAmatOutput(false);
        cluster.doIntStrucOutput(false);
        cluster.doVtkOutput(false);
        cluster.doRestartOutput(true);
        cluster.doFStatOutput(false);
        cluster.doEneOutput(false);
        cluster.setRadiusCluster(0.9 * radiiCluster[i]); // 0.8* è un safety factor
        cluster.setPosition(positionsCluster[i]);
        cluster.setClusterId(particleHandler.getNextGroupId());
        cluster.setSizeDispersityParticle(sizeDispersityParticleIntra);
        cluster.setRadiusParticle(radiusParticle);
        cluster.setParticleSpecies(species);
        cluster.solve();


        importParticlesAs(cluster.particleHandler, cluster.interactionHandler, speciesVector[i]);

        //!total particle volume
        Mdouble rCluster = cluster.getMeanClusterRadius();
        Mdouble rForMassFraction = rCluster - 0.5 * radiusParticle;
        totalParticleVolume += constants::pi*4*pow(rForMassFraction,3)*(0.58+3*pow(0.58,2)*speciesVector[i]->getPenetrationDepthMax())/3;

        std::cout << "Created cluster number " << i+1 << "/" << nClusters << std::endl;

    }

    nParticles = particleHandler.getSize();

    //! Setting timeStep and timeMax
    //setTimeStep(species -> getCollisionTime(1)/50); // (collision time)/50.0
    setTimeStep( speciesHandler.getMixedObject(speciesVector[0], speciesVector[1])->getCollisionTime(1) / 50 );
    setTimeMax(1000);
    fileOutputCount = 50;
    setSaveCount( 100*fileOutputCount );
    fileOutputTime = fileOutputCount*getTimeStep();

    //! Setting Domain (already set in the constructor)
    //setDomain({-boxSize,-boxSize,-boxSize},{boxSize,boxSize,boxSize});

    //!Simulation variables
    t0 = getTime();
    compressionTime = 0.5 * fileOutputTime;
    relaxationTime = 2.5 * fileOutputTime;
    stage = 1;
    dampingCoefficient = 8e-6;
    isStrainRateControlled = true;
    strainRatePerIteration = -2e-6;
    relaxPoint1 = 1.5;
    relaxPoint2 = 0.98;
    counterTotalForces = 1;
    counterTractionForce = 0;

    //! needed for results
    counterCompression = 0;
    primoValore = true;
    results.reserve(100);
    stress.reserve(100);

    //! Printing some values of interest
    std::cout << "number of particles: " << particleHandler.getSize() << std::endl;
    std::cout << "number of interactions: " << interactionHandler.getSize() << std::endl;



    std::cout << getXMax() << "   " << getXMin() << std::endl;
    std::cout << getYMax() << "   " << getYMin() << std::endl;
    std::cout << getZMax() << "   " << getZMin() << std::endl;


    //! Creating walls
    makeWalls();

    //! Creating IsoCFile
    makeIsoCFile();

    random.randomise();
}

void DoublePorosityStressStrainC::actionsAfterTimeStep(){

    //! Computes void ratio (and volumeBox)
    computeVoidRatioAndPorosity();

    //dampVelocities();

    if (getNumberOfTimeSteps()==1) {
        makeDataAnalysis();
        printTime();
        writeIsoCFile();
        makeOverlFile();
        writeOverlFile();
    }

    if ( fmod(getTime(),fileOutputTime) < getTimeStep() ) {
        makeDataAnalysis();
        writeIsoCFile();
        printTime();
    }

    if (stage == 1 && e<=eValues[0]){
        makeDataAnalysis();
        makeOverlFile();
        writeOverlFile();
        stressE0pre = uniStress;
        stage = 2;
    }

    //! Relax (si espande il dominio verticale, no strain rate control)
    if ( stage == 2 && uniStress > compressionRatio * stressE0pre)
    {
        makeOverlFile();
        writeOverlFile();
        boundaryHandler.clear();
        w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(-isoStrainDot, 0, 0, 0, -isoStrainDot, 0, 0, 0, -isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), isStrainRateControlled);
        boundaryHandler.copyAndAddObject(w0);
        stage = 3;
    }

    if (stage == 3 && e > 0.78){
        makeOverlFile();
        writeOverlFile();
        boundaryHandler.clear();
        w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), isStrainRateControlled);
        boundaryHandler.copyAndAddObject(w0);
        t0 = getTime();
        stage = 4;
    }

    if ( stage == 4 && e <= eValues[counterCompression] )
    {
        if (primoValore){
            makeDataAnalysis();
            stressE0 = uniStress;
            primoValore = false;
            makeOverlFile();
            writeOverlFile();
        }
        makeDataAnalysis();
        writeToIsoCFileTarget();
        printTime();
        results.push_back(uniStress/stressE0);
        stress.push_back(uniStress);
        counterCompression++;
    }

    if ( stage == 4 && e <= eValues[99]) {
        makeDataAnalysis();
        printTime();
        writeIsoCFile();
        makeOverlFile();
        writeOverlFile();

        std::ostringstream resultsName;
        resultsName << getName() << ".results";
        resultsFile.open(resultsName.str(), std::ios::out);

        std::ostringstream stressName;
        stressName << getName() << ".stress";
        stressFile.open(stressName.str(), std::ios::out);

        for (int i = 0; i < 100; ++i) {
            resultsFile <<
                        std::fixed << std::setprecision(3) << results[i] <<
                        std::endl;

            stressFile <<
                       std::fixed << std::setprecision(3) << stress[i] <<
                       std::endl;

        }

        setTimeMax(getTime());
    }





/*
    if ( stage == 1 && (getTime() - t0 > compressionTime || e<=eValues[0]) ){

        makeDataAnalysis();
        printTime();
        writeIsoCFile();
        //! Sets relaxation boundary
        boundaryHandler.clear();
        w0.set(Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0),
               Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), true);
        boundaryHandler.copyAndAddObject(w0);
        t0 = getTime();
        stage = 2;
    }

    if ( stage == 2 && getTime() - t0 > relaxationTime ){

        makeDataAnalysis();
        printTime();
        writeIsoCFile();
        //! Sets relaxation boundary
        boundaryHandler.clear();
        if (e<=eValues[0]) {
            w0.set(Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, isoStrainDot),
                   Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), true);
            boundaryHandler.copyAndAddObject(w0);
            stage = 3;
        }
        else {
            w0.set(Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0),
                   Matrix3D(isoStrainDot, 0, 0, 0, isoStrainDot, 0, 0, 0, isoStrainDot),
                   Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), true);
            boundaryHandler.copyAndAddObject(w0);
            stage = 1;
        }
        t0 = getTime();
    }

    if (stage == 3){

        if (e <= eValues[counterCompression]) {
            makeDataAnalysis();
            writeToIsoCFileTarget();
            printTime();
            results.push_back(uniStress);
            counterCompression++;
        }

        if (getTime() - t0 > compressionTime) {
            makeDataAnalysis();
            printTime();
            writeIsoCFile();
            //! Sets relaxation boundary
            boundaryHandler.clear();
            w0.set(Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0),
                   Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), true);
            boundaryHandler.copyAndAddObject(w0);
            t0 = getTime();
            stage = 2;
        }
    }

    if (e <= eValues[99]) {
        makeDataAnalysis();
        printTime();
        writeIsoCFile();

        std::ostringstream resultsName;
        resultsName << getName() << ".results";

        resultsFile.open(resultsName.str(), std::ios::out);

        for (int i = 0; i < 100; ++i) {
            resultsFile <<
                        std::fixed << std::setprecision(3) << results[i] <<
                        std::endl;
        }
        setTimeMax(getTime());
    }
*/
}

void DoublePorosityStressStrainC::printTime() const{

    std::cout <<

              "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
              ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
              ", Stress: " << std:: fixed << std::setprecision(2) << std::setw(11) << totalStress <<
              ", uniStress: " << std:: fixed << std::setprecision(2) << std::setw(11) << uniStress <<
              ", uniStressIntra: " << std:: fixed << std::setprecision(2) << std::setw(11) << uniStressIntra <<
              ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
              ", Inertia = " << std::scientific << std::setprecision(2) << std::setw(8) << inertia <<
              ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber_ <<
              ", n = " << std::fixed << std::setprecision(3) << voidRatio <<
              ", e = " << std::fixed << std::setprecision(3) << e <<
              ", em = " << std::fixed << std::setprecision(3) << em <<
              ", eM = " << std::fixed << std::setprecision(3) << eM <<
              ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlapIntra_ <<
              ", dMean = " << std:: fixed << std::setw(7) << meanRelativeOverlapIntra_ <<
              ", dMax = " << maxRelativeOverlapIntra_ <<
              ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlapInter_ <<
              ", dMean = " << std::fixed << std::setw(7) << meanRelativeOverlapInter_ <<
              ", dMax = " << maxRelativeOverlapInter_ <<
              ", rPS = " <<  std::setprecision(2) << relativePoreSize <<
              ", PercPos = " << std::setw(5) << 100*counterTractionForce/counterTotalForces << "%" <<
              ", NC = " << std::scientific << std::setprecision(2) << std::setw(5) << NC

              << std::endl;
}

void DoublePorosityStressStrainC::setSpecies() {


    speciesHandler.clear();


    Mdouble constantMass = densityParticle * 4 * constants::pi * pow(radiusParticle, 3) / 3;
    Mdouble collisionTimeIntra = sqrt( constantMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffnessIntra ) );
    Mdouble collisionTimeInter = sqrt( constantMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffnessInter ) );

    std::cout << "HERE!" << loadingStiffnessIntra << std::endl;
    std::cout << "HERE!" << unLoadingStiffnessMaxIntra << std::endl;

    //! Setting Species
    species = new LinearPlasticViscoelasticFrictionSpecies;
    species -> setConstantRestitution(true);
    species -> setDensity(densityParticle);
    species -> setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, restitutionCoefficient, constantMass);
    species -> setUnloadingStiffnessMax(species->getLoadingStiffness()*unLoadingStiffnessMaxIntra/loadingStiffnessIntra);
    species -> setCohesionStiffness(0);
    species -> setPenetrationDepthMax(penetrationDepthMax);

    species -> setSlidingFrictionCoefficient(0);
    species -> setSlidingStiffness(0);
    species -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
    speciesHandler.copyAndAddObject(species);

    std::cout << "HERE!" << species->getUnloadingStiffnessMax() << std::endl;
    std::cout << "HERE!" << unLoadingStiffnessMaxIntra << std::endl;
    std::cout << "HERE!" << species->getUnloadingStiffnessMax()*cohesionStiffness/unLoadingStiffnessMaxIntra << std::endl;
    std::cout << "HERE!" << cohesionStiffness
              << std::endl;


    speciesVector.reserve(nClusters);
    // Single Cluster Species
    for (int i = 0; i < nClusters; ++i) {
        speciesVector[i] = new LinearPlasticViscoelasticFrictionSpecies;
        speciesVector[i] -> setDensity(densityParticle);
        speciesVector[i] -> setConstantRestitution(true);
        speciesVector[i] -> setCollisionTimeAndRestitutionCoefficient(collisionTimeIntra, restitutionCoefficient, constantMass);

        speciesVector[i] -> setUnloadingStiffnessMax(speciesVector[i]->getLoadingStiffness()*unLoadingStiffnessMaxIntra/loadingStiffnessIntra);
        speciesVector[i] -> setCohesionStiffness(speciesVector[i]->getUnloadingStiffnessMax()*cohesionStiffness/unLoadingStiffnessMaxIntra);
        speciesVector[i] -> setPenetrationDepthMax(penetrationDepthMax);



        speciesVector[i] -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
        speciesVector[i] -> setSlidingStiffness(speciesVector[i] -> getLoadingStiffness() * relativeTangentialStiffness);
        speciesVector[i] -> setSlidingDissipation(speciesVector[i] -> getDissipation()*2.0/7.0);
        speciesHandler.addObject(speciesVector[i]);
        speciesVector.push_back(speciesVector[i]);

    }
    std::cout << "Set single cluster species" << std::endl;


    // Cluster mixed species
    for (int i = 0; i < nClusters; ++i) {
        for (int j = 0; j < nClusters; ++j) {
            if (i!=j){

                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setConstantRestitution(true);
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setCollisionTimeAndRestitutionCoefficient(collisionTimeInter, restitutionCoefficient, constantMass);
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setUnloadingStiffnessMax( speciesHandler.getMixedObject(speciesVector[i], speciesVector[j])->getLoadingStiffness() );
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setCohesionStiffness(0);
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setPenetrationDepthMax(1);

                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingStiffness( speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> getLoadingStiffness()*relativeTangentialStiffness );
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingDissipation( speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> getDissipation()*2.0/7.0 );
            }
        }
    }

}

void DoublePorosityStressStrainC::setSpecies2() {


    std::cout << "Changing species: adding sliding friction." << std::endl;

    // Single Cluster Species
    for (int i = 0; i < nClusters; ++i) {
        speciesVector[i] -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    }

    // Cluster mixed species
    for (int i = 0; i < nClusters; ++i) {
        for (int j = 0; j < nClusters; ++j) {
            if (i!=j){
                speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
            }
        }
    }

}

void DoublePorosityStressStrainC::makeWalls(){

    isoStrainDot = strainRatePerIteration/getTimeStep();
    boundaryHandler.clear();
    w0.setHandler(&boundaryHandler);
    w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(isoStrainDot, 0, 0, 0, isoStrainDot, 0, 0, 0, isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), isStrainRateControlled);
    boundaryHandler.copyAndAddObject(w0);
}

void DoublePorosityStressStrainC::makeIsoCFile() {

    std::ostringstream cdatName;
    cdatName << getName() << ".isoc";

    isoCFile.open(cdatName.str(), std::ios::out);

    isoCFile << "CLUSTER DATA AND INFORMATION" << std::endl << std::endl;
    isoCFile << "radiusParticle: " << std:: scientific << std::setprecision(2) << radiusParticle << std::endl;
    isoCFile << "densità: " << std:: scientific << std::setprecision(2) << densityParticle << std::endl;
    isoCFile << "sizeDispersityParticle: " << std::defaultfloat << sizeDispersityParticleIntra << std::endl;
    isoCFile << "nParticles: " << nParticles << std::endl;
    isoCFile << "boxSize: " << std::defaultfloat << (getXMax()-getXMin()) << std::endl;
    isoCFile << "Damping coefficient: " << std::defaultfloat << dampingCoefficient << std::endl;
    isoCFile << "Strain Rate Controlled: " << isStrainRateControlled << std::endl;
    isoCFile << "strainRatePerIteration: " << strainRatePerIteration << std::endl;
    isoCFile << "Loading stiffness intra: " << std:: scientific << std::setprecision(2) << loadingStiffnessIntra << std::endl;
    isoCFile << "Unloading stiffness max intra: " << unLoadingStiffnessMaxIntra << std::endl;
    isoCFile << "Cohesion stiffness: " << cohesionStiffness << std::endl;
    isoCFile << "Loading stiffness inter: " << loadingStiffnessInter << std::endl;


    isoCFile << "COMPRESSION INFORMATION" << std::endl << std::endl;
    isoCFile <<

             "   Time" << std::setw(13) <<
             "Time Max" << std::setw(16) <<
             "Stress" << std::setw(16) <<
             "uniStress" << std::setw(16) <<
             "uniStressIntra" << std::setw(15) <<
             "E_ratio" << std::setw(11) <<
             "Inertia" << std::setw(12) <<
             "cN" << std::setw(12) <<
             "n" << std::setw(12) <<
             "e" << std::setw(12) <<
             "em" << std::setw(12) <<
             "eM" << std::setw(14) <<
             "dMinIntra" << std::setw(14) <<
             "dMeanIntra" << std::setw(14) <<
             "dMaxIntra"<< std::setw(14) <<
             "dMinInter" << std::setw(14) <<
             "dMeanInter" << std::setw(14) <<
             "dMaxInter" << std::setw(6) <<
             "rPS" << std::setw(6) <<
             "stage"

             << std::endl;
}

void DoublePorosityStressStrainC::makeOverlFile() {

    std::ostringstream overlName;
    overlName << getName() << "_e_" << e  << "_stage_" << stage << ".overl";

    overlFile.open(overlName.str(), std::ios::out);

    overlFile << "OVERLAP VS FORCE" << std::endl << std::endl;

}


void DoublePorosityStressStrainC::dampVelocities()
{
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p)
    {
        (*p)->addForce(-dampingCoefficient * (*p)->getVelocity() * pow((*p)->getRadius()/radiusParticle,3) );
    }
}


//! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
void DoublePorosityStressStrainC::makeDataAnalysis()
{
    //! Computes coordination number and overlaps.
    computeCoordinationNumberAndOverlaps();

    //! Computes stress and inertia.
    computeStress();

    //! Computes void dimension
    //computeVoidDimensions();
}

//! Computes coordination number and overlaps
void DoublePorosityStressStrainC::computeCoordinationNumberAndOverlaps()
{
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
    for (int i=0; i < particleHandler.getNumberOfRealObjects(); i++)
    {
        meanCoordinationNumber_ += (particleHandler.getObject(i)->getInteractions()).size();
    }
    meanCoordinationNumber_ /= particleHandler.getSize();

    //! loops over every interaction to compute maximum, mean and minimum relative particle overlap.
    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
    {
        /*!
         * \details the relative overlap is computed as an average of the relative overlap on the two particles.
         *          rO = ( O/R1 + O/R2 ) / 2.
         */
        relativeOverlap = ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                          ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
        relativeOverlap /= 2;
        if ( (*i) -> getP() -> getGroupId() == (*i) -> getI() -> getGroupId() ){
            meanRelativeOverlapIntra_ += relativeOverlap;
            if (relativeOverlap > maxRelativeOverlapIntra_)
                maxRelativeOverlapIntra_ = relativeOverlap;
            if (relativeOverlap < minRelativeOverlapIntra_)
                minRelativeOverlapIntra_ = relativeOverlap;
            counterOverlapIntra++;
        } else{
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
void DoublePorosityStressStrainC::computeVoidRatioAndPorosity()
{
    volumeBox = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
    voidRatio = 1 - totalParticleVolume / volumeBox;
    e = (volumeBox - totalParticleVolume) / totalParticleVolume;
    Mdouble mFIntra = 0.58+3*pow(0.58,2)*meanRelativeOverlapIntra_;
    em = (1-mFIntra)/mFIntra;
    eM = (e-em)/(1+em);
}

//! computes stress
void DoublePorosityStressStrainC::computeStress()
{
    stressXX = 0.0;
    stressYY = 0.0;
    stressZZ = 0.0;
    stressXY = 0.0;
    stressXZ = 0.0;
    stressYZ = 0.0;
    stressZZIntra = 0.0;
    inertia = 0.0;
    NC = 0.0;

    counterTotalForces = 0;
    counterTractionForce = 0;

    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {

        stressXX += (*i)->getForce().X * (*i)->getNormal().X * (*i)->getDistance();
        stressYY += (*i)->getForce().Y * (*i)->getNormal().Y * (*i)->getDistance();
        stressZZ += (*i)->getForce().Z * (*i)->getNormal().Z * (*i)->getDistance();
        //Vec3D v = (*i)->getI()->getPosition()-(*i)->getP()->getPosition();

        if ( (*i)->getP()->getGroupId() == (*i)->getI()->getGroupId() ) {
            stressZZIntra += (*i)->getForce().Z * (*i)->getNormal().Z * (*i)->getDistance();
            counterTotalForces++;
            Vec3D v = -(*i)->getNormal();
            if ( (*i)->getForce().X*v.X + (*i)->getForce().Y*v.Y + (*i)->getForce().Z*v.Z > 0 ) counterTractionForce++;
        }

    }

    stressXX /= volumeBox;
    stressYY /= volumeBox;
    stressZZ /= volumeBox;
    stressZZIntra /= volumeBox;

    totalStress = (stressXX+stressYY+stressZZ)/3;
    uniStress = stressZZ;
    uniStressIntra = stressZZIntra;

    inertia = -2*radiusParticle*isoStrainDot/sqrt(totalStress/densityParticle);

    NC = stressZZIntra * radiusParticle / unLoadingStiffnessMaxIntra;

}

//! Computes void dimension
void DoublePorosityStressStrainC::computeVoidDimensions()
{
    bool inserted = false;
    Mdouble radius = 2 * radiusParticle;

    while (!inserted)
    {
        if (particleInsertionSuccessful(radius))
        {
            inserted = true;
        }
        else
        {
            radius *= 0.99;
        }
    }
    relativePoreSize = radius/radiusParticle;
}

bool DoublePorosityStressStrainC::particleInsertionSuccessful(Mdouble r) {

    int insertionFailCounter = 0;
    Vec3D particlePosition;
    SphericalParticle p0;

    //! setup of particle properties and initial conditions (besides position)
    p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
    p0.setRadius(r);
    p0.setSpecies(speciesHandler.getObject(0));
    p0.setGroupId(0);

    while (insertionFailCounter < 1e6)
    {

        particlePosition.X = random.getRandomNumber(getXMin(), getXMax());
        particlePosition.Y = random.getRandomNumber(getYMin(), getYMax());
        particlePosition.Z = random.getRandomNumber(getZMin(), getZMax());

        p0.setPosition(particlePosition);

        if (checkParticleForInteraction(p0))
        {
            return true;
        }

        insertionFailCounter++;
    }

    return false;
}

void DoublePorosityStressStrainC::writeOverlFile() {

    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {
        Mdouble force = (*i)->getForce().X*(*i)->getNormal().X + (*i)->getForce().Y*(*i)->getNormal().Y + (*i)->getForce().Z*(*i)->getNormal().Z;
        overlFile << (*i)->getOverlap()/radiusParticle << std::setw(18) << std::scientific << std::setprecision(5) << force <<std::endl;
    }
    overlFile.close();
}


void DoublePorosityStressStrainC::writeIsoCFile() {

    isoCFile <<

             std::scientific << std::setprecision(2) << std::setw(7) << getTime() <<
             std::fixed << std::setprecision(2) << std::setw(13) << getTimeMax() <<
             std::fixed << std::setprecision(2) << std::setw(16) << totalStress <<
             std::fixed << std::setprecision(2) << std::setw(16) << uniStress <<
             std::fixed << std::setprecision(2) << std::setw(16) << uniStressIntra <<
             std::scientific << std::setprecision(2) << std::setw(15) << getKineticEnergy()/getElasticEnergy() <<
             std::scientific << std::setprecision(2) << std::setw(11) << inertia <<
             std::fixed << std::setw(12) << meanCoordinationNumber_ <<
             std::fixed << std::setprecision(3) << std::setw(12) << voidRatio <<
             std::fixed << std::setprecision(3) << std::setw(12) << e <<
             std::fixed << std::setprecision(3) << std::setw(12) << em <<
             std::fixed << std::setprecision(3) << std::setw(12) << eM <<
             std::fixed << std::setprecision(5) << std::setw(14) << minRelativeOverlapIntra_ <<
             std::setw(14) << meanRelativeOverlapIntra_ <<
             std::setw(14) << maxRelativeOverlapIntra_ <<
             std::setw(14) << minRelativeOverlapInter_ <<
             std::setw(14) << meanRelativeOverlapInter_ <<
             std::setw(14) << maxRelativeOverlapInter_ <<
             std::setw(6) << std::setprecision(2) << relativePoreSize <<
             std::setw(5) << 100*counterTractionForce/counterTotalForces << "%" <<
             std::scientific << std::setprecision(2) << std::setw(5) << NC <<
             std::fixed << std::setw(5) << stage

             << std::endl;
}

//! writes on the .isoc file
void DoublePorosityStressStrainC::writeToIsoCFileTarget()
{
    isoCFile <<

             std::scientific << std::setprecision(2) << std::setw(7) << getTime() <<
             std::fixed << std::setprecision(2) << std::setw(13) << getTimeMax() <<
             std::fixed << std::setprecision(2) << std::setw(16) << totalStress <<
             std::fixed << std::setprecision(2) << std::setw(16) << uniStress <<
             std::fixed << std::setprecision(2) << std::setw(16) << uniStressIntra <<
             std::scientific << std::setprecision(2) << std::setw(15) << getKineticEnergy()/getElasticEnergy() <<
             std::scientific << std::setprecision(2) << std::setw(11) << inertia <<
             std::fixed << std::setw(12) << meanCoordinationNumber_ <<
             std::fixed << std::setprecision(3) << std::setw(12) << voidRatio <<
             std::fixed << std::setprecision(3) << std::setw(12) << e <<
             std::fixed << std::setprecision(3) << std::setw(12) << em <<
             std::fixed << std::setprecision(3) << std::setw(12) << eM <<
             std::fixed << std::setprecision(5) << std::setw(14) << minRelativeOverlapIntra_ <<
             std::setw(14) << meanRelativeOverlapIntra_ <<
             std::setw(14) << maxRelativeOverlapIntra_ <<
             std::setw(14) << minRelativeOverlapInter_ <<
             std::setw(14) << meanRelativeOverlapInter_ <<
             std::setw(14) << maxRelativeOverlapInter_ <<
             std::setw(6) << std::setprecision(2) << relativePoreSize <<
             std::setw(5) << 100*counterTractionForce/counterTotalForces << "%" <<
             std::scientific << std::setprecision(2) << std::setw(5) << NC <<
             std::fixed << std::setw(5) << stage <<
             std::setw(23) << "Target: " << counterCompression + 1

             << std::endl;

}


void DoublePorosityStressStrainC::makeContactModelFile(){

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
