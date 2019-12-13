//
// Created by paolo on 6-10-19.
//
//Ciao!


#include "BaseDoublePorosityStressStrain.h"




BaseDoublePorosityStressStrain::BaseDoublePorosityStressStrain(){

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

BaseDoublePorosityStressStrain::BaseDoublePorosityStressStrain(Mdouble pD, Mdouble lS, Mdouble rTS, Mdouble sFC, Mdouble bS, Mdouble sRPTS, std::string n){

    particleDispersity = pD;
    loadingStiffness = lS;
    relativeTangentialStiffness = rTS;
    slidingFrictionCoefficient = sFC;
    boxSize = bS;
    nParticles = 100;
    strainRatePerTimeStep = sRPTS;
    setName(n);

}

BaseDoublePorosityStressStrain::BaseDoublePorosityStressStrain(Mdouble pD, Mdouble tC, Mdouble r, Mdouble rTS, Mdouble sFC, Mdouble bS, int nP, std::string n){

    particleDispersity = pD;
    collisionTime = tC;
    restitutionCoefficient = r;
    relativeTangentialStiffness = rTS;
    slidingFrictionCoefficient = sFC;
    boxSize = bS;
    nParticles = nP;
    setName(n);

}

BaseDoublePorosityStressStrain::BaseDoublePorosityStressStrain(Mdouble pD, Mdouble lS, Mdouble rTS, Mdouble sFC, Mdouble bS, std::string n){

    particleDispersity = pD;
    loadingStiffness = lS;
    relativeTangentialStiffness = rTS;
    slidingFrictionCoefficient = sFC;
    boxSize = bS;
    nParticles = 50;
    setName(n);

}

BaseDoublePorosityStressStrain::~BaseDoublePorosityStressStrain(){

    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

//! Initializing variables, boxSize, species, radii, walls, particles, timeStep, output, name
void BaseDoublePorosityStressStrain::setupInitialConditions()
{
    // OUTPUT FILES
    dataFile.setFileType(FileType::ONE_FILE);
    restartFile.setFileType(FileType::ONE_FILE);
    fStatFile.setFileType(FileType::NO_FILE);
    eneFile.setFileType(FileType::NO_FILE);

    //random.randomise();

    radiusParticle =5.5e-5;
    rMin = 10*radiusParticle;
    rMax = 0;

    setDomain({ -boxSize/2, -boxSize/2, -boxSize/2}, { boxSize/2, boxSize/2, boxSize/2});

    std::cout << std::endl << std::endl << "SETTING SPECIES" << std::endl << std::endl;
    setSpecies();

    std::cout << std::endl << std::endl << "INSERTING PARTICLES" << std::endl << std::endl;
    insertParticles();
    std::cout << std::endl << std::endl << "rMin: " << rMin << ", rMax: " << rMax << ", ratio: " << rMax/rMin << std::endl << std::endl;

    totalParticleVolume = particleHandler.getVolume();

    // Strain
    strainRatePerTimeStep = -1e-6;
    isoStrainDot = strainRatePerTimeStep/getTimeStep();
    std::cout << isoStrainDot << std::endl;
    std::cout << std::endl << std::endl << "CREATING WALLS" << std::endl << std::endl;
    createWalls();

    t0 = getTime();

    isocOutputTimeInterval = getTimeStep() * 50;
    std::cout << isocOutputTimeInterval << std::endl;
    std::cout << getTimeStep() << std::endl;
    setSaveCount( 100 * floor(isocOutputTimeInterval/getTimeStep()) );
    setTimeMax(1000);

    setXBallsAdditionalArguments("-v0 -p 10");

    dampingCoefficient = 0.001;
    relaxationTime = 100 * isocOutputTimeInterval;
    compressionTime = relaxationTime / 50;
    counterCompression = 0;
    counterRelaxation = 0;
    stage = 1;

}

//! Makes data analysis and writes to isotropic compression file
void BaseDoublePorosityStressStrain::actionsBeforeTimeStep()
{

    makeDataAnalysis();

    if (stage == 1) {


        if ( voidRatio < 0.45 )
        {
            printTime();
            //! Sets relaxation boundary
            boundaryHandler.clear();
            w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), Matrix3D(0,0,0,0,0,0,0,0,0), true);
            boundaryHandler.copyAndAddObject(w0);
            t0 = getTime();
            stage = 2;
        }
    }

    if (stage == 2) {

        if ( maxRelativeOverlap < 2e-4)
        {
            printTime();
            clusterPositions.reserve(nParticles);
            clusterRadii.reserve(nParticles);
            for (int i = 0; i < nParticles; ++i) {
                clusterPositions.push_back(particleHandler.getObject(i)->getPosition());
                clusterRadii.push_back(particleHandler.getObject(i)->getRadius());
            }
            setTimeMax(getTime());
        }
    }

}

//! Damps velocities
void BaseDoublePorosityStressStrain::actionsAfterTimeStep()
{
    dampVelocities();
}

//! Prints time and values of interest
void BaseDoublePorosityStressStrain::printTime() const
{

    switch (stage)
    {
        case 1: std::cout << "Isotropic compression: ";
            break;

        case 2: std::cout << "Relaxation " << counterCompression << ":         ";
            break;

        case 3: std::cout << "Uniaxial compression:  ";
            break;

        default: std::cout << "Final values:         ";
            break;
    }
    std::cout <<

              "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
              ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
              ", Stress: " << std:: fixed << std::setprecision(2) << std::setw(11) << totalStress <<
              ", uniStress: " << std:: fixed << std::setprecision(2) << std::setw(11) << uniStress <<
              ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
              ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber <<
              ", n = " << std::fixed << std::setprecision(3) << voidRatio <<
              ", e = " << std::fixed << std::setprecision(3) << e <<
              ", floaters = " << std::fixed << std::setw(3) << floaters <<
              ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlap <<
              ", dMean = " << std:: fixed << std::setw(7) << meanRelativeOverlap <<
              ", dMax = " << maxRelativeOverlap

              << std::endl;
}

void BaseDoublePorosityStressStrain::dampVelocities()
{
    for (int i=particleHandler.getNumberOfObjects()-1; i>=0; i--)
    {
        BaseParticle* P;
        P=particleHandler.getObject(i);
        P->addForce(-dampingCoefficient * P->getVelocity() * pow(P->getRadius()/radiusParticle,3) );
    }
}

void BaseDoublePorosityStressStrain::setSpecies(){

    speciesHandler.clear();

    //! Particles properties
    Mdouble densityParticle = 2700; //! ottenuta considerando la densità della bentonite = 2150 Kg/m^3 ed una vF di 0.37 (corrispondente al momento in cui E_ratio crolla)
    Mdouble penetrationDepthMax = 1;
    //! This the mass of a particle having radius radiusParticle (constant restitution will be true)
    Mdouble constantMass = densityParticle * 4 * constants::pi * pow(radiusParticle, 3) / 3;
    collisionTime = sqrt( constantMass * ( pow(constants::pi, 2) + pow(log(restitutionCoefficient), 2) ) / ( 2 * loadingStiffness ) );

    //! Linear Species
    LinearPlasticViscoelasticFrictionSpecies* linearSpecies;

    linearSpecies = new LinearPlasticViscoelasticFrictionSpecies;
    linearSpecies -> setConstantRestitution(true);
    linearSpecies -> setDensity(densityParticle);
    linearSpecies -> setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, 1);
    linearSpecies -> setUnloadingStiffnessMax(linearSpecies->getLoadingStiffness());
    linearSpecies -> setCohesionStiffness(0);
    linearSpecies -> setPenetrationDepthMax(penetrationDepthMax);

    linearSpecies -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    linearSpecies -> setSlidingStiffness(linearSpecies -> getLoadingStiffness()*relativeTangentialStiffness);
    linearSpecies -> setSlidingDissipation(linearSpecies -> getDissipation()*2.0/7.0);

    speciesHandler.copyAndAddObject(linearSpecies);

    Mdouble contactTimeOverTimeStep = 50;
    setTimeStep(linearSpecies -> getCollisionTime(constantMass)/contactTimeOverTimeStep);
    std::cout << "collision time: " << linearSpecies -> getCollisionTime(constantMass) << std::endl;
    std::cout << "timeStep: " << std::setprecision(4) << getTimeStep() << std::endl;
    std::cout << "tC_PP/dt, at least: " << std::setprecision(4) << linearSpecies -> getCollisionTime(constantMass)/getTimeStep() << std::endl << std::endl;



/*
    //! herzian Species
    HertzianViscoelasticFrictionSpecies* herzianSpecies;
    herzianSpecies = new HertzianViscoelasticFrictionSpecies;
    herzianSpecies -> setDensity(densityParticle);
    herzianSpecies -> setConstantRestitution(true);
    herzianSpecies -> setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(collisionTime, restitutionCoefficient, restitutionCoefficient, 1);
    herzianSpecies -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
    herzianSpecies -> setSlidingStiffness(loadingStiffness*relativeTangentialStiffness);
    herzianSpecies -> setSlidingDissipation(herzianSpecies -> getDissipation()*2.0/7.0);
    speciesHandler.copyAndAddObject(herzianSpecies);


    setTimeStep(herzianSpecies -> getCollisionTime(radiusParticle*2,densityParticle,(1-isoStrainDot)*boxSize/(2*getTimeStep()))/contactTimeOverTimeStep);
    std::cout << "collision time: " << linearSpecies -> getCollisionTime(constantMass) << std::endl;
    std::cout << "timeStep: " << std::setprecision(4) << getTimeStep() << std::endl;
    std::cout << "tC_PP/dt, at least: " << std::setprecision(4) << linearSpecies -> getCollisionTime(constantMass)/getTimeStep() << std::endl << std::endl;



*/


}

//! inserts the particles in the simulation domainss
void BaseDoublePorosityStressStrain::insertParticles()
{
    int nParticlesInserted = 0;
    // insert particles till the desired amount is reached
    while (nParticlesInserted < nParticles)
    {    // nParticle inserted corresponds to the index of the particle which is being inserted
        if (particleInsertionSuccessful(nParticlesInserted))
        {
            nParticlesInserted++;
            std::cout << "Set position of particle n. " << nParticlesInserted << "/" << nParticles<< std::endl;
        }
        else
        {
            logger(ERROR, "Cannot insert all particles, try to decrease the value of initialSolidFraction.");
        }
    }
    std::cout << "PARTICLE INSERTION TERMINATED" << std::endl << std::endl;
}

//! Check if the particle is inserted whitout interactions
bool BaseDoublePorosityStressStrain::particleInsertionSuccessful(int n)
{
    // initialization of parameters, n is the particle index
    int insertionFailCounter = 0;
    Vec3D particlePosition;
    SphericalParticle p0;

    // setup of particle properties and initial conditions (besides position)
    p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
    // Raggio calcolato in modo che il raggio medio sia radiusParticle
    Mdouble r = radiusParticle * random.getRandomNumber(2/(1+particleDispersity),2*particleDispersity/(1+particleDispersity));
    //Mdouble r = radiusParticle * random.getRandomNumber(1,particleDispersity);
    p0.setRadius(r);
    p0.setSpecies(speciesHandler.getObject(0));
    if (r>rMax)
        rMax = r;
    if (r<rMin)
        rMin = r;



    // in this cycle a random position inside of a sphere contained in the bounding box is taken
    // if the particle is not in contact with the others then insert it, otherwise increment the fail counter
    // the maximum number of failed attempts is capped to 1000

    while (insertionFailCounter < 1000)
    {

        //particlePosition.X = ( getXMax() - getXMin() - 2*radiusParticle*2*particleDispersity/(1+particleDispersity) )*random.getRandomNumber(0.0,1.0) + getXMin() + radiusParticle*2*particleDispersity/(1+particleDispersity);
        //particlePosition.Y = ( getYMax() - getYMin() - 2*radiusParticle*2*particleDispersity/(1+particleDispersity) )*random.getRandomNumber(0.0,1.0) + getYMin() + radiusParticle*2*particleDispersity/(1+particleDispersity);
        //particlePosition.Z = ( getZMax() - getZMin() - 2*radiusParticle*2*particleDispersity/(1+particleDispersity) )*random.getRandomNumber(0.0,1.0) + getZMin() + radiusParticle*2*particleDispersity/(1+particleDispersity);

        particlePosition.X = ( getXMax() - getXMin() )*random.getRandomNumber(0.0,1.0) + getXMin();
        particlePosition.Y = ( getYMax() - getYMin() )*random.getRandomNumber(0.0,1.0) + getYMin();
        particlePosition.Z = ( getZMax() - getZMin() )*random.getRandomNumber(0.0,1.0) + getZMin();


        p0.setPosition(particlePosition);

        if (checkParticleForInteraction(p0))
        {
            particleHandler.copyAndAddObject(p0);
            return true;
        }
        insertionFailCounter++;
    }
    return false;
}

//! Creating stress strain control walls
void BaseDoublePorosityStressStrain::createWalls()
{
    boundaryHandler.clear();
    w0.setHandler(&boundaryHandler);
    w0.set(Matrix3D(0,0,0,0,0,0,0,0,0), Matrix3D(isoStrainDot, 0, 0, 0, isoStrainDot, 0, 0, 0, isoStrainDot), Matrix3D(0,0,0,0,0,0,0,0,0), true);
    boundaryHandler.copyAndAddObject(w0);
}


//! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
void BaseDoublePorosityStressStrain::makeDataAnalysis()
{
    //! Computes coordination number and overlaps.
    computeCoordinationNumberAndOverlaps();

    //! Computes void ratio (and volumeBox)
    computeVoidRatioAndPorosity();

    //! Computes the stress.
    computeStress();
}

//! Computes coordination number and overlaps
void BaseDoublePorosityStressStrain::computeCoordinationNumberAndOverlaps()
{
    meanForceOnInteraction = 0.0;
    meanCoordinationNumber = 0.0;
    maxRelativeOverlap = 0.0;
    meanRelativeOverlap = 0.0;
    minRelativeOverlap = inf;
    Mdouble relativeOverlap;
    floaters = 0;

    //! loops over each particle to compute mean coordination number and center of mass.
    for (int i=0; i < particleHandler.getNumberOfRealObjects(); i++)
    {
        meanCoordinationNumber += (particleHandler.getObject(i)->getInteractions()).size();
        if (particleHandler.getObject(i)->getInteractions().empty())
            floaters ++;
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
        if (relativeOverlap > maxRelativeOverlap)
            maxRelativeOverlap = relativeOverlap;

        if (relativeOverlap < minRelativeOverlap)
            minRelativeOverlap = relativeOverlap;
    }
    meanForceOnInteraction/= interactionHandler.getSize();
    meanRelativeOverlap /= interactionHandler.getSize();
}

//! computes voidRatio and porosity (e)
void BaseDoublePorosityStressStrain::computeVoidRatioAndPorosity()
{
    volumeBox = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
    voidRatio = 1 - totalParticleVolume / volumeBox;
    e = (volumeBox - totalParticleVolume) / totalParticleVolume;
}

//! computes stress
void BaseDoublePorosityStressStrain::computeStress()
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


