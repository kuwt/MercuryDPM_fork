//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "BaseCluster.h"

/*!
 * \details Default constructor.
 */
BaseCluster::BaseCluster(){
    logger(DEBUG, "BaseCluster::BaseCluster() finished");
}

/*!
 * \details Default destructor.
 */
BaseCluster::~BaseCluster(){
    logger(DEBUG, "BaseCluster::BaseClusterr() finished");
}

/*
 * ----------------------------------------------------------------------
 *               FUNCTIONS: setters and getters
 * ----------------------------------------------------------------------
 */

/*!
 * \details This returns the value of position_, which is the position in which the cluster will be inserted (it is
 *          also the centre of mass of the cluster).
 */
Vec3D BaseCluster::getPosition() const{
    return position_;
}

/*!
 * \details This sets the value of position_, which is the position in which the cluster will be inserted (it is
 *          also the centre of mass of the cluster).
 */
void BaseCluster::setPosition(Vec3D p){
    position_ = p;
}


/*!
 * \details This returns the value of the ratio between collision time and time step
 */
Mdouble BaseCluster::getCollisionTimeOverTimeStep() const {
    return collisionTimeOverTimeStep_;
}

/*!
 * \details This sets the collisionTimeOverTimeStep_ number.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setCollisionTimeOverTimeStep(Mdouble cTOTS) {
    if (cTOTS < 45 && cTOTS > 0)
        logger(WARN, "collisionTimeOverTimeStep = % is too small: consider setting it greater or equal than 50.", cTOTS);
    else if (cTOTS <= 0)
        logger(ERROR, "collisionTimeOverTimeStep = % must be grater than zero.", cTOTS);
    collisionTimeOverTimeStep_ = cTOTS;
}

/*!
 * \details This returns the value of particles' radius if there's no dispersity in size.
 *          If the dispersity is not 1, this is the value on which all radii will be computed.
 */
Mdouble BaseCluster::getRadiusParticle() const {
    return radiusParticle_;
}

/*!
 * \details This sets the value of particles' radius if there's no dispersity in size.
 *          If the dispersity is not 1, this is the value on which all radii will be computed.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setRadiusParticle(Mdouble rP) {
    if (rP <= 0)
        logger(ERROR, "radiusParticle must be greater than zero. radiusParticle = %", rP);
    else {
        radiusParticle_ = rP;
        setRadiusParticle_ = true;
    }
}

/*!
 * \details This returns the value of particles' dispersity in size.
 */
Mdouble BaseCluster::getSizeDispersityParticle() const {
    return sizeDispersityParticle_;
}

/*!
 * \details This sets the value of particles' dispersity in size.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setSizeDispersityParticle(Mdouble sDP) {
    if (sDP < 1)
        logger(ERROR, "sizeDispersityParticle must be greater or equal than 1. sizeDispersityParticle = %", sDP);
    else sizeDispersityParticle_ = sDP;
}

/*!
 * \details This returns the value of the number of particles in the cluster.
 */
int BaseCluster::getNumberOfParticles() const {
    return nParticles_;
}

/*!
 * \details This sets the value of the number of particles in the cluster.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setNumberOfParticles(int nP) {
    if (nP < 0)
        logger(ERROR, "nParticles must be grater than zero. nParticles = %", nP);
    else {
        nParticles_ = nP;
        setNumberOfParticles_ = true;
    }
}

/*!
 * \details This sets the value of the Radius of the cluster.
 *          A boolean (setRadiusCluster_) is also set to true: this is done because the user has to set exactly two
 *          among radiusCluster_, radiusCluster_ and radiusParticle_.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setRadiusCluster(Mdouble rCR) {
    if (rCR <= 0)
        logger(ERROR,
               "relativeClusterRadius is smaller than 0. relativeClusterRadius = %",
               rCR);
    else
        radiusCluster_ = rCR;

    setRadiusCluster_ = true;
}

/*!
 * \details This return the value of the obtained mass fraction (could be a little different than the desired).
 *          In particular, if the internal structure output is computed, the value for the mass fraction will be solidFractionIntStruct_
 *          (computed within the internal structure and more precise), otherwise it will be solidFraction_.
 */
Mdouble BaseCluster::getFinalMassFraction(){
    if (isIntStrucOutputOn_)
        return solidFractionIntStruct_;
    else
        return solidFraction_;
}

/*!
 * \details This returns the value of the cluster ID.
 */
unsigned int BaseCluster::getClusterId() const {
    return idCluster_;
}

/*!
 * \details This sets the value of the cluster ID.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setClusterId(unsigned int iC) {
    if (iC < 0)
        logger(WARN, "idCluster = % is less than zero.", iC);
    idCluster_ = iC;
}

/*!
 * \details This returns the value of the velocity damping modulus.
 */
Mdouble BaseCluster::getVelocityDampingModulus() const {
    return velocityDampingModulus_;
}

/*!
 * \details This sets the value of the velocity damping modulus.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setVelocityDampingModulus(Mdouble vDM) {
    if (vDM < 0 || vDM > 1)
        logger(ERROR, "velocityDampingModulus must be grater than zero and less than 1. velocityDampingModulus = %", vDM);
    velocityDampingModulus_ = vDM;
}

/*!
 * \details This returns the value of the length of the number of particles used to compute the internal structure.
 */
int BaseCluster::getNumberOfInternalStructurePoints() const {
    return nInternalStructurePoints_;
}

/*!
 * \details This sets the value of the length of the number of particles used to compute the internal structure.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setNumberOfInternalStructurePoints(int gL) {
    if (gL <= 0)
        logger(ERROR, "nInternalStructurePoints_ must be grater than zero. nInternalStructurePoints_ = %", gL);
    nInternalStructurePoints_ = gL;
}

/*!
 * \details This returns the value of the value of the energy ratio threshold under which the process can be considered static, and so over.
 */
Mdouble BaseCluster::getEnergyRatioTolerance() const {
    return energyRatioTolerance_;
}

/*!
 * \details This sets the value of the value of the energy ratio threshold under which the process can be considered static, and so over.
 *          In addition to that, it checks if the value is acceptable.
 */
void BaseCluster::setEnergyRatioTolerance(Mdouble eRT) {
    if (eRT <= 0)
        logger(ERROR, "energyRatioTolerance must be grater than zero. energyRatioTolerance = %", eRT);
    energyRatioTolerance_ = eRT;
}

/*!
 * \details This returns the species of the particle.
 */
LinearPlasticViscoelasticFrictionSpecies *BaseCluster::getParticleSpecies() const {
    return particleSpecies_;
}

/*!
 * \details This sets the species of the particle.
 */
void BaseCluster::setParticleSpecies(LinearPlasticViscoelasticFrictionSpecies *pS) {
    particleSpecies_ = pS;
}

/*!
 * \details This returns the velocity of the cluster after creation.
 */
Vec3D BaseCluster::getVelocity() {
    return clusterVelocity_;
}

/*!
 * \details This sets the velocity of the cluster after creation.
 */
void BaseCluster::setVelocity(Vec3D v) {
    clusterVelocity_ = v;
}

/*!
 * \details This returns the bool variable that defines whether the cluster data output is written or not.
 */
bool BaseCluster::isCdatOutputOn() const {
    return isCdatOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster data output will be written or not.
 */
void BaseCluster::doCdatOutput(bool iCOO) {
    isCdatOutputOn_ = iCOO;
}

/*!
 * \details This returns the bool variable that defines whether the cluster overlap output is written or not.
 */
bool BaseCluster::isOverlOutputOn() const {
    return isOverlOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster overlap output will be written or not.
 */
void BaseCluster::doOverlOutput(bool iOOO) {
    isOverlOutputOn_ = iOOO;
}

/*!
 * \details This returns the bool variable that defines whether the cluster adjacency matrix output is written or not.
 */
bool BaseCluster::isAmatOutputOn() const {
    return isAmatOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster adjacency matrix output will be written or not.
 */
void BaseCluster::doAmatOutput(bool iAOO) {
    isAmatOutputOn_ = iAOO;
}

/*!
 * \details This returns the bool variable that defines whether the cluster internal structure output is written or not.
 */
bool BaseCluster::isIntStrucOutputOn() const {
    return isIntStrucOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster internal structure output will be written or not.
 */
void BaseCluster::doIntStrucOutput(bool iISOO) {
    isIntStrucOutputOn_ = iISOO;
}

/*!
 * \details This returns the bool variable that defines whether the cluster vtk output is written or not.
 */
bool BaseCluster::isVtkOutputOn() const {
    return isVtkOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster vtk output will be written or not.
 */
void BaseCluster::doVtkOutput(bool iVOO) {
    isVtkOutputOn_ = iVOO;
}

/*!
 * \details This returns the bool variable that defines whether the cluster restart output is written or not.
 */
bool BaseCluster::isRestartOutputOn() const {
    return isRestartOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster restart output will be written or not.
 */
void BaseCluster::doRestartOutput(bool r) {
    isRestartOutputOn_ = r;
}

/*!
 * \details This returns the bool variable that defines whether the cluster fStat output is written or not.
 */
bool BaseCluster::isFStatOutputOn() const {
    return isFStatOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster fStat output will be written or not.
 */
void BaseCluster::doFStatOutput(bool fS) {
    isFStatOutputOn_ = fS;
}

/*!
 * \details This returns the bool variable that defines whether the cluster ene output is written or not.
 */
bool BaseCluster::isEneOutputOn() const {
    return isEneOutputOn_;
}

/*!
 * \details This sets the bool variable that defines whether the cluster ene output will be written or not.
 */
void BaseCluster::doEneOutput(bool e) {
    isEneOutputOn_ = e;
}

/*!
 * \details this returns meanClusterRadius (radius of an ideal perfectly spherical cluster, there's no setter).
 */
Mdouble BaseCluster::getMeanClusterRadius() {
    return meanClusterRadius_;
}

/*!
 * \details this returns the average overlap.
 */
Mdouble BaseCluster::getAverageOverlap(){
    return meanRelativeOverlap_;
}

/*
 * ----------------------------------------------------------------------
 *               FUNCTIONS: overridden mercury3D functions
 * ----------------------------------------------------------------------
 */

/*!
 * \details Overrides DPMBase function setupInitialConditions(): in this initial conditions for the problem are set, and they are (in order):
 *          -checking that particleHandler and speciesHandler are empty,
 *          -defining cluster parameters starting from the user input,
 *          -all particles radii,
 *          -particle species,
 *          -domain limits,
 *          -time step,
 *          -particles initial positions (insertParticles),
 *          -initial time (t0),
 *          -maximum force applied on particles (maximumForceModulus),
 *          -time duration of applying (and releasing) force (forceTuningDuration) and dissipation time (dissipationDuration_),
 *          -time max,
 *          -output time for files and print time function (fileOutputTimeInterval),
 *          -time interval on which force is increased or decreased (forceTuningInterval),
 *          -value of forceDampingModulus_,
 *          -time interval on which velocity is damped (velocityTuningInterval),
 *          -xballs settings,
 *          -vtk settings,
 *          -standard mercury output data settings (dataFile, restartFile, fStatFile and eneFile),
 *          -name settings,
 *          -creation of cluster data file (if needed),
 *          -creation of overlap file (if needed),
 *          -stage of computation.
 */
void BaseCluster::setupInitialConditions()
{
    logger(VERBOSE, "CREATING CLUSTER");

    if (particleHandler.getSize()>0){
        logger(WARN, "ParticleHandler was not empty");
        particleHandler.clear();
    }

    if (speciesHandler.getSize()>0){
        logger(WARN, "speciesHandler was not empty");
        speciesHandler.clear();
    }
    // Defining cluster parameters starting from the user input.
    // In order to use this class the user has to set exactly 2 values among radiusParticle,
    // radiusCluster and numberOfParticles.
    if (setNumberOfParticles_*setRadiusParticle_ + setNumberOfParticles_*setRadiusCluster_ + setRadiusCluster_*setRadiusParticle_ != 1)
        logger(ERROR, "Please set exactly two values among radiusParticle, radiusCluster and numberOfParticles."
                      " radiusParticle = %, radiusCluster = %, numberOfParticles = %.",
                       radiusParticle_, radiusCluster_, nParticles_);

    // If the user sets the cluster radius and the radius of a single particle,
    // the number of particles has to be computed:
    if (setRadiusCluster_ && setRadiusParticle_) {
        if (radiusCluster_ < 2 * radiusParticle_)
            logger(VERBOSE, "clusterRadius is small compared to the radius of a single particle:"
                         " consider setting clusterRadius >= 2 * radiusParticle."
                          "clusterRadius = %, radiusParticle = %.", radiusCluster_, radiusParticle_);
        // relative cluster radius
        Mdouble hatR = radiusCluster_ / radiusParticle_;
        // mass fraction of the cluster in the limit of phi = 0.
        Mdouble eps0 = 0.58;
        // The number of particles (N) per cluster given the relative cluster radius (hatR) and penetration depth max (phi)
        // can be computed as: N = ( hatR / (1 - eps0*phi) )^3 * eps0.
        // It is very important to notice that this formula is accurate only if sliding friction is set to 0.5 and relative
        // tangential stiffness is set to 0.3 while creating the cluster. Different values do not guarantee accuracy.

        nParticles_ =
                round (std::pow( hatR / (1.0 - eps0 * particleSpecies_->getPenetrationDepthMax() ), 3) * eps0);
        logger(VERBOSE, "Number of particles: %.\n", nParticles_);
    }

    // If the user sets the cluster radius and the number of particles,
    // the radius of a single particle has to be computed:
    if (setRadiusCluster_ && setNumberOfParticles_) {
        // mass fraction of the cluster in the limit of phi = 0.
        Mdouble eps0 = 0.58;
        // The radius of a single particle (r) composig the cluster given the cluster radius (R), penetration depth max (phi)
        // and the number of particles (N) can be obtained as: r = R / ( cbrt(N/eps0) * (1-eps0*phi) ).
        // It is very important to notice that this formula is accurate only if sliding friction is set to 0.5 and relative
        // tangential stiffness is set to 0.3 while creating the cluster. Different values do not guarantee accuracy.
        radiusParticle_ = radiusCluster_ / (cbrt(nParticles_/eps0) * (1-eps0*particleSpecies_->getPenetrationDepthMax()));

        logger(VERBOSE, "radius particle: %.\n", radiusParticle_);
    }

    logger(VERBOSE, "SETTING RADII");
    setRadii();

    logger(VERBOSE, "SETTING SPECIES");
    setSpecies();

    logger(VERBOSE, "SETTING DOMAIN LIMITS");
    setDomainLimits();

    logger(VERBOSE, "COMPUTING TIME STEP");
    calculateTimeStep();

    logger(VERBOSE, "PARTICLE INSERTION");
    insertParticles();

    logger(VERBOSE, "Number of particles: %.", particleHandler.getSize());

    t0_ = getTime();

    // \details deltaStar
    Mdouble deltaStar = particleSpecies_->getPenetrationDepthMax() * particleSpecies_->getUnloadingStiffnessMax()
                        / (particleSpecies_->getUnloadingStiffnessMax() - particleSpecies_->getLoadingStiffness());

    /*
     * \brief maximum force modulus applied on particles (this value is then multiplied by the relative distance from
     *          force center d, which is d = D/r).
     *          It is the force necessary to get a overlap of deltaStar (computed above this description).
     *          If constant restitution is true then it is also multiplied by the mass of the biggest particle,
     *          or if  MERCURY_USE_MPI it is the biggest possible mass computed taking into account particle dispersity,
     *          i.e. the mass of a particle having radius
     *          r = radiusParticle_ * 2 * sizeDispersityParticle_ / (1 + sizeDispersityParticle_).
     *          (In order to get the right value of loading stiffness it should be multiplied by the mass of the smallest
     *          particle; multiplying it for the biggest mass here, instead, is for safety).
     */

#ifdef MERCURY_USE_MPI
    maximumForceModulus_ = radiusParticle_ * 2 * sizeDispersityParticle_ / (1 + sizeDispersityParticle_) * deltaStar *
                           particleSpecies_->getLoadingStiffness()
                           * (particleSpecies_->getConstantRestitution() ?
                              particleSpecies_->getMassFromRadius(radiusParticle_ * 2 * sizeDispersityParticle_ / (1 + sizeDispersityParticle_))
                              :
                              1);
#else
    maximumForceModulus_ = radiusParticle_ * 2 * sizeDispersityParticle_ / (1 + sizeDispersityParticle_) * deltaStar *
    particleSpecies_->getLoadingStiffness()
                           * (particleSpecies_->getConstantRestitution() ?
                              particleHandler.getLargestParticle()->getMass()
                              :
                              1);
#endif

    /*
     * \details
     * The time needed for a particle to cover a distance of x is  sqrt(2 * x * m / F), when a constant force F is applied.
     * In order to hit another particle, a single particle has to travel a distance of about boxSize/4 (boxSize is the domain length).
     * As a value for x it is used boxSize_ (instead of boxSize_/4 as a safety factor).
     * As a value for F it is used maximumForceModulus/50 which is half of maximumForceModulus/25 (in this way it is calculated
     * a measure of the time needed with a force which linearly increases from 0 to maximumForceModulus/25, again dividing for 25
     * is done for safety). As a value for m it is used the  mass of the biggest particle, or if MERCURY_USE_MPI it is
     * the mass computed taking into account size dispersity, so the mass
     * of a particle having radius r = radiusParticle_ * 2 * sizeDispersityParticle_ / (1 + sizeDispersityParticle_), for safety.
     * The factor 5* outside the sqrt is another safety factor empirically determined: with this values the computation is fast
     * and the obtained results are very similar to the ones obtained if longer times would be set.
     * All this safety factors are needed because this is not the exact value of time needed but a measure of it: the problem
     * indeed is quite complicated given that particles can also rearrange during compression and so they will eventually
     * move in a non radial direction and for this reason they will need more time to settle.
     */
#ifdef MERCURY_USE_MPI
    forceTuningDuration_ =
            5 * sqrt(100 * boxSize_ * particleSpecies_->getMassFromRadius(radiusParticle_ * 2 * sizeDispersityParticle_ /
                                                        (1 + sizeDispersityParticle_)) / maximumForceModulus_);
#else
    forceTuningDuration_ =
            5 * sqrt(100 * boxSize_ * particleHandler.getLargestParticle()->getMass() / maximumForceModulus_);
#endif

    //\details Maximum possible time duration of dissipation (i.e. duration of dissipation if energy ration tollerance not reached).
    dissipationDuration_ = forceTuningDuration_/2;

    // Compression + Decompression + Dissipation = 2 * Compression + Dissipation
    clusterTimeMax_ = 2 * forceTuningDuration_ + dissipationDuration_;
    setTimeMax(clusterTimeMax_);

    fileOutputTimeInterval_ = forceTuningDuration_ / 100;

    setSaveCount(floor(fileOutputTimeInterval_ / getTimeStep()));

    forceTuningInterval_ = getTimeStep();

    forceDampingModulus_ = 0.95;

    velocityDampingInterval_ = getTimeStep();

    setXBallsAdditionalArguments("-v0 -p 10");

    setParticlesWriteVTK(isVtkOutputOn());

    dataFile.setFileType(FileType::ONE_FILE);

    restartFile.setFileType(isRestartOutputOn() ? FileType::ONE_FILE : FileType::NO_FILE);

    fStatFile.setFileType(isFStatOutputOn() ? FileType::ONE_FILE : FileType::NO_FILE);

    eneFile.setFileType(isEneOutputOn() ? FileType::ONE_FILE : FileType::NO_FILE);

    // Name setting
    std::ostringstream name;
    name << "Cluster_ID_" << idCluster_;
    setName(name.str());

    if (isCdatOutputOn()) {
        logger(VERBOSE, "CREATING .cdat FILE\n");
        makeCdatFile();
    }

    if (isOverlOutputOn()) {
        logger(VERBOSE, "CREATING .overl FILE\n");
        makeOverlFile();
    }

    logger(VERBOSE, "ACTIVATING CENTRAL FORCES\n");

    /*
     * \details
     * Stage defines in which phase of the calculation the program is:
     * stage = 1: compressing particles and increasing force
     * stage = 2: releasing force
     * stage = 3: waiting for the system to be static.
     */
    stage_ = 1;
}

/*!
 * \details In this the process takes place: in particular after each time step the following tasks are computed (in order):
 *          -data analysis: calculating values of interest,
 *          -writing to cluster data file,
 *          -writing to overlap file,
 *          -increasing central force and damping velocity (stage == 1), decreasing central force and damp velocity (stage == 2),
 *              damp central force and damp velocity until the system is static (stage == 3).
 */
void BaseCluster::actionsAfterTimeStep()
{
    makeDataAnalysis();

    if ( isCdatOutputOn() && fmod(getTime()-t0_, fileOutputTimeInterval_) < getTimeStep() )
        writeToCdatFile();

    if ( isOverlOutputOn() && fmod(getTime()-t0_, fileOutputTimeInterval_) < getTimeStep() )
        writeToOverlFile();

    /*
     * \brief If stage == 1 force is linearly increased and velocities are damped for a time T = forceTuningDuration.
     *        If t > T, stage is set to 2.
     */
    if (stage_ == 1)
    {
        if (getTime() - t0_ < forceTuningDuration_)
        {
            if (fmod(getTime() - t0_, forceTuningInterval_) < getTimeStep())
                increaseForce();

            if (fmod(getTime() - t0_, velocityDampingInterval_) < getTimeStep())
                dampVelocities();

            applyCentralForce();
        }
        else
        {
            logger(VERBOSE, "DECREASING CENTRAL FORCE");

            t0_ = getTime();
            stage_++;
        }
    }

    /*
     * \brief If stage == 2 force is linearly decreased and velocities are damped for a time duration of T = forceTuningDuration.
     *        If t > T, stage is set to 3.
     */
    if (stage_ == 2)
    {
        if (getTime() - t0_ < forceTuningDuration_)
        {
            if (fmod(getTime() - t0_, forceTuningInterval_) < getTimeStep())
                decreaseForce();

            if (fmod(getTime() - t0_, velocityDampingInterval_) < getTimeStep())
                dampVelocities();

            applyCentralForce();
        }
        else
        {
            logger(VERBOSE, "DISSIPATING ENERGY");

            t0_ = getTime();
            fileOutputTimeInterval_ *= 2;
            setSaveCount( floor(fileOutputTimeInterval_/getTimeStep()) );
            stage_++;
        }
    }

    /*
     * \details If stage == 3 force is exponentially decreased and velocities are damped for a time duration of T = dissipationDuration_.
     *          If t>T or if the energy ratio is below the minimum threshold calculation is concluded and a few last operation are computed.
     *          They are:
     *          -timeMax is set to getTime(), in order to stop the calculation,
     *          -stage is set to 4 (if the energy threshold is not reached stage will remain 3 (because the simulation is stopped by the previous
     *           definition of timeMax): if this happens the user gets a warning, see actionsAfterSolve()).
     *          If MERCURY_USE_MPI this process lasts for a time T = dissipationDuration_ - getTimeStep().
     */
    if (stage_ == 3)
    {
        // \brief Now force is damped and not decreased.
        if (fmod(getTime() - t0_, forceTuningInterval_) < getTimeStep())
            dampForce();

        if (fmod(getTime() - t0_, velocityDampingInterval_) < getTimeStep())
            dampVelocities();

        applyCentralForce();
#ifdef MERCURY_USE_MPI
        if (getTime() - t0_ > dissipationDuration_ - getTimeStep())
#else
        if (getKineticEnergy() / getElasticEnergy() < energyRatioTolerance_ ||
            getTime() - t0_ > dissipationDuration_)
#endif
        {
            printTime();
            logger(VERBOSE, "ENERGY DISSIPATED\n");

            // stage++ now is a flag used to understand if the dissipation procedure has been completed.
            stage_++;

            setTimeMax(getTime());
        }
    }

}

/*!
 * \details In this:
 *          -data analysis is computed and cluster data and cluster overlap file are written for the last time.
 *          -if stage == 3 that means that the simulation is over but the energy ratio threshold is not reached:
 *              for this reason the user gets a warning.
 *          -adjacency matrix file is created and written (if needed),
 *          -internal structure file is created and written (if needed),
 *          -gnuplot file is created and written (if needed),
 *          -all particles are centred around the center of mass and their velocity is set,
 *           for this reason the user gets a warning,
 *          -finally cluster data file and cluster overlap file are closed.
 */
void BaseCluster::actionsAfterSolve()
{

    makeDataAnalysis();

    if ( isCdatOutputOn() )
        writeToCdatFile();

    if ( isOverlOutputOn() )
        writeToOverlFile();

    if (stage_ == 3)
        logger(WARN, "Dissipation process not completed: final energyRatioTollerance_ = %."
                     "Try to increase energyRatioTollerance_ or decreasing velocityDampingModulus.",
                     getKineticEnergy()/getElasticEnergy());

    if (isAmatOutputOn())
    {
        logger(VERBOSE, "CREATING ADJACENCY MATRIX FILE\n");
        createAdjacencyMatrix();
        makeAmatFile();
        writeAmatFile();
    }

    if (isIntStrucOutputOn())
    {
        logger(VERBOSE, "COMPUTING INTERNAL STRUCTURE FILE\n");
        computeInternalStructure();
    }

    if(isOverlOutputOn())
    {
        logger(VERBOSE, "COMPUTING GNUPLOT FILE\n");
        makeGnuplotFile();
    }


    /*
     * \brief with this loop all particles are moved so that center of mass == position_ and their velocity is set.
     */
    for (auto i = particleHandler.begin(); i != particleHandler.end(); ++i){
        (*i)->setPosition( (*i)->getPosition() + position_ - centerOfMass_ );
        (*i)->setVelocity( clusterVelocity_ );
    }

    logger(VERBOSE, "CLUSTER CREATED.\n");

    if (isCdatOutputOn())
        cdatFile_.close();

    if (isOverlOutputOn())
        overlFile_.close();
}

/*!
 * \details In this all variables needed by the program for restarting are written on a .restart file
 *          Most of the needed variables are already saved by MercuryBase::write, which for this reason is included.
 */
void BaseCluster::write(std::ostream& os, bool writeAllParticles) const
{
    writeAllParticles = true;
    MercuryBase::write(os, writeAllParticles);

    os <<
            "position " << position_ << " " <<
            "stage " << stage_ << " " <<
            "t0 " << t0_
            << std::endl <<
            "idCluster " << idCluster_ << " " <<
            "nParticles " << nParticles_ << " " <<
            "radiusParticle " << radiusParticle_ << " " <<
            "massParticle " << massParticle_ << " " <<
            "sizeDispersityParticle " << sizeDispersityParticle_ << " " <<
            "totalParticleVolume " << totalParticleVolume_
            << std::endl <<
            "maximumForceModulus " << maximumForceModulus_ << " " <<
            "forceTuningInterval " << forceTuningInterval_ << " " <<
            "forceTuningDuration " << forceTuningDuration_ << " " <<
            "velocityDampingInterval " << velocityDampingInterval_ << " " <<
            "velocityDampingModulus " << velocityDampingModulus_ << " " <<
            "energyRatioTolerance " << energyRatioTolerance_ << " " <<
            "forceDampingModulus " << forceDampingModulus_ << " " <<
            "forceModulus " << forceModulus_
            << std::endl <<
            "isCdatOutputON " << isCdatOutputOn_ << " " <<
            "isOverlOutputOn " << isOverlOutputOn_ << " " <<
            "isAmatOutputOn " << isAmatOutputOn_ << " " <<
            "isIntStrucOutputOn " << isIntStrucOutputOn_
            << std::endl;
    if( isIntStrucOutputOn() )
    os <<   "nInternalStructurePoints " << nInternalStructurePoints_ << std::endl;

    os <<   "isRestartOutputOn " << isRestartOutputOn_ << " " << //This is obviously on because otherwise restart
                                                                 // process would not take place.
                                                                 //For now it is saved but could eventually be removed.
            "isFStatOutputOn " << isFStatOutputOn_ << " " <<
            "isEneOutputOn " << isEneOutputOn_ << std::endl;
}

/*!
 * \details In this all variables needed by the program for restarting are read from a .restart file.
 *          As MercuryBase::write was included in ClusterDPM::write, here MercuryBase::read is included.
 */
void BaseCluster::read(std::istream& is, ReadOptions opt)
{
    MercuryBase::read(is);

    std::stringstream line;
    std::string dummy;

    helpers::getLineFromStringStream(is, line);
    line    >> dummy >> position_
            >> dummy >> stage_
            >> dummy >> t0_;
    helpers::getLineFromStringStream(is, line);
    line    >> dummy >> idCluster_
            >> dummy >> nParticles_
            >> dummy >> radiusParticle_
            >> dummy >> massParticle_
            >> dummy >> sizeDispersityParticle_
            >> dummy >> totalParticleVolume_;
    helpers::getLineFromStringStream(is, line);
    line    >> dummy >> maximumForceModulus_
            >> dummy >> forceTuningInterval_
            >> dummy >> forceTuningDuration_
            >> dummy >> velocityDampingInterval_
            >> dummy >> velocityDampingModulus_
            >> dummy >> energyRatioTolerance_
            >> dummy >> forceDampingModulus_
            >> dummy >> forceModulus_;
    helpers::getLineFromStringStream(is, line);
    line    >> dummy >> isCdatOutputOn_
            >> dummy >> isOverlOutputOn_
            >> dummy >> isAmatOutputOn_
            >> dummy >> isIntStrucOutputOn_;
    if(isIntStrucOutputOn() )
    {
    helpers::getLineFromStringStream(is, line);
    line >> dummy >> nInternalStructurePoints_;
    }
    line    >> dummy >> isRestartOutputOn_
            >> dummy >> isFStatOutputOn_
            >> dummy >> isEneOutputOn_;
}

/*!
 * \details In this all variables needed by the program for restarting are initialized.
 *          DPMBase::readRestartFile() is included.
 */
void BaseCluster::actionsOnRestart()
{
    readRestartFile();

    fileOutputTimeInterval_ = forceTuningDuration_ / 100;

    dissipationDuration_ = forceTuningDuration_/2;

    setSaveCount( floor(fileOutputTimeInterval_/getTimeStep()) );

    forceTuningInterval_ = getTimeStep();

    velocityDampingInterval_ = getTimeStep();

    setXBallsAdditionalArguments("-v0 -p 10");

    setParticlesWriteVTK(isVtkOutputOn());

    dataFile.setFileType(FileType::ONE_FILE);

    restartFile.setFileType(isRestartOutputOn()?FileType::ONE_FILE:FileType::NO_FILE);

    fStatFile.setFileType(isFStatOutputOn()?FileType::ONE_FILE:FileType::NO_FILE);

    eneFile.setFileType(isEneOutputOn()?FileType::ONE_FILE:FileType::NO_FILE);

    if (isCdatOutputOn())
    {
        std::ostringstream cdatName;
        cdatName << getName() << ".cdat";
        cdatFile_.open(cdatName.str(), std::ios::app);
    }

    if (isOverlOutputOn())
    {
        std::ostringstream overlName;
        overlName << getName() << ".overl";
        overlFile_.open(overlName.str(), std::ios::app);
    }

    logger(VERBOSE, "CALCULATION RESTARTED\n");
}

/*!
 * \details Overrides DPMBase printTime(): this way stage and variables of interest are shown. They are (in order):
 *          -computation progress,
 *          -energy ratio,
 *          -mean coordination number,
 *          -mean cluster radius,
 *          -mass fraction,
 *          -force modulus,
 *          -minimum relative overlap,
 *          -mean relative overlap,
 *          -maximum realative overlap,
 *          -center of mass,
 *          -number of particles.
 *
 *          printTime output is set to VERBOSE in order not to have too much output. If the user needs it, it is enough
 *          to set it to INFO.
 */
void BaseCluster::printTime() const
{
    std::ostringstream printTime;
    switch (stage_)
    {
        case 1: printTime << "Compression progress: " << std::setw(3) << int( ceil( 100 * (getTime() - t0_) / forceTuningDuration_ ) ) << "%, ";
            break;

        case 2: printTime << "Decompression progress: " << std::setw(3) << int( ceil( 100 * (getTime() - t0_) / forceTuningDuration_ ) )<< "%, ";
            break;

        case 3: printTime << "Dissipating energy: ";
            break;

        default: printTime << "Final values: ";
            break;
    }
    printTime <<
              "E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
              ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber_ << ", rMean = " << std::scientific << meanClusterRadius_ <<
              ", mF = " << std::fixed << std::setprecision(3) << solidFraction_ << ", Force Modulus = " << std::scientific << forceModulus_ <<
              ", dMin = " << std::fixed << std::setw(7) << std::setprecision(5) << minRelativeOverlap_ << ", dMean = " << std::setw(7) << meanRelativeOverlap_ <<
              ", dMax = " << maxRelativeOverlap_ << ", centerOfMass = " << std::scientific << std::setprecision(5) << std::setw(13) << centerOfMass_.X
                          << std::setw(14) << centerOfMass_.Y << std::setw(14) << centerOfMass_.Z <<
              " nParticles: " << particleHandler.getSize();
    logger(VERBOSE, printTime.str());
}



/*
 * ----------------------------------------------------------------------
 *           FUNCTIONS: functions inside setupInitialConditions
 * ----------------------------------------------------------------------
 */

/*!
 * \details Sets all radii according to particleRadius and sizeDispersityParticle_.
 */
void BaseCluster::setRadii()
{
    totalParticleVolume_ = 0;
    smallestRadius_ = constants::inf;
    for (int i = 0; i < nParticles_; ++i)
    {
        // This is the actual radius of the i-th particle
        radii_.push_back( random.getRandomNumber( radiusParticle_*2/(1+sizeDispersityParticle_) , radiusParticle_*2*sizeDispersityParticle_/(1+sizeDispersityParticle_) ) );
        // Computing totalParticleVolume (this is done because this value is needed before particle insertion).
        totalParticleVolume_ += 4.0*constants::pi*pow(radii_[i],3.0)/3.0;
        //computing the smallest radius (needed in computeTimeStep(), which is needed before particle insertion)
        if (radii_[i] < smallestRadius_) smallestRadius_ = radii_[i];
    }
}

/*!
 * \details adds particlePureSpaces (set by the user or by default) to the species handler.
 *          After this, values of stiffnesses, restitution coefficient and collision time are printed.
 */
void BaseCluster::setSpecies()
{
    if (particleSpecies_ == nullptr)
        logger(ERROR, "Species not set.");

    speciesHandler.copyAndAddObject(particleSpecies_);
    /*
     * \brief mass of the particle which has radius radiusParticle if constantRestitution(false) or
     * radiusParticle_*2/(1+sizeDispersityParticle_) otherwise.
     * It is set now and used various time in the code (for example for stiffnesses, collision time
     * and restitution coefficient below).
     */
     massParticle_ = particleSpecies_ -> getConstantRestitution()?
                    speciesHandler.getObject(0) -> getMassFromRadius(radiusParticle_*2/(1+sizeDispersityParticle_))
            :
                    speciesHandler.getObject(0) -> getMassFromRadius(radiusParticle_);

    std::ostringstream printSpecies;
    /*
     * If constantRestitution(true) loading, unloading, and cohesion stiffness are multiplied by the mass of a particle
     * whose radius is radiusParticle_*2/(1+sizeDispersityParticle_),
     * which is the mass that should be used to compute collision time.
     */
     printSpecies << "loadingStiffness: " << std::scientific << particleSpecies_ -> getLoadingStiffness()
                                                            * (particleSpecies_->getConstantRestitution()?massParticle_:1) << std::endl
                  << "unloadingStiffnessMax: " << particleSpecies_ -> getUnloadingStiffnessMax()
                                              * (particleSpecies_->getConstantRestitution()?massParticle_:1) << std::endl
                  << "cohesionStiffness: " << particleSpecies_ -> getCohesionStiffness()
                                          * (particleSpecies_->getConstantRestitution()?massParticle_:1) << std::endl
                  << "restitutionCoefficient: " << std::fixed << particleSpecies_ -> getRestitutionCoefficient(massParticle_) << std::endl
                  << "collisionTime: " << std::scientific << std::setprecision(3) << particleSpecies_ -> getCollisionTime(massParticle_) << std::endl;
    logger(VERBOSE, printSpecies.str());
}

/*!
 * \details Sets domain limits of the simulation.
 * Using spherical coordinates and the insertion process defined in particleInsertionSuccessful the final solid
 *  fraction obtained varies from 0.11 to 0.14 (empirically calculated, even for high size dispersity, which anyway
 *  shouldn't be contemplated in this field). For this reason the domain size (here called boxSize_) is computed as
 *  the cubic root of the total particle volume over 0.1.
 */
void BaseCluster::setDomainLimits()
{
    Mdouble initialSolidFraction = 0.1;
    boxSize_ = cbrt(totalParticleVolume_/initialSolidFraction);
    std::ostringstream printDomainLimits;
    printDomainLimits << "Cubic size " << boxSize_ << std::endl;
    logger(VERBOSE, printDomainLimits.str());

    setDomain(-0.5*boxSize_*Vec3D(1,1,1) + position_, 0.5*boxSize_*Vec3D(1,1,1) + position_);
}

/*!
 * \details Calculates the time step over the smallest particle.
 *          After this the time step and the ratio between collision time and time step are printed.
 */
void BaseCluster::calculateTimeStep()
{
    // if constantRestitution(true) mass for collision time will be automatically set to 1, otherwise the smallest particle's mass will be used
    setTimeStep(particleSpecies_ -> getCollisionTime(speciesHandler.getObject(0)->getMassFromRadius(smallestRadius_))/collisionTimeOverTimeStep_);

    // printing values of timeStep and ratio between collisionTime and timeStep
    std::ostringstream printTimeStep;
    printTimeStep << "timeStep: " << std::setprecision(4) << getTimeStep() << std::endl
                  << "cT/tS, at least: " << std::fixed << std::setprecision(1)
                  << particleSpecies_ -> getCollisionTime(speciesHandler.getObject(0)->getMassFromRadius(smallestRadius_))/getTimeStep()
                  << std::endl;
    logger(VERBOSE, printTimeStep.str());
}


/*!
 * \details Inserts particles inside the domain with a while cycle.
 */
void BaseCluster::insertParticles()
{
    int nParticlesInserted = 0;

    while (nParticlesInserted < nParticles_)
    {
        /* nParticleInserted corresponds to the index of the particle which is being inserted:
          * for example if no particle has been inserted yet nParticlesInserted=0 which is the
          * index of the first particle, and so on. For this reason this variable is the input
          * for particleInsertionSuccessful.
          */
        if (particleInsertionSuccessful(nParticlesInserted))
        {
            nParticlesInserted++;
        }
        else
        {
            logger(ERROR, "Cannot insert all particles, try to decrase the value of initialSolidFraction in "
                          "BaseCluster::setDomainLimits();\n"
                          "Inserted %/% particles.", nParticlesInserted, nParticles_);
        }
    }
    logger(VERBOSE, "PARTICLE INSERTION TERMINATED SUCCESSFULLY\n");
}

/*!
 * \details Creates the cluster data output file.
 *          In the first lines of the file some important infos about the cluster are written. They are (in order):
 *          -Ratio between collision time and time step,
 *          -radiusParticle,
 *          -size dispersity of particles,
 *          -density of the particles
 *          -number of particles,
 *          -cluster ID (group ID),
 *          -safety factor of the cluster radius,
 *          -sliding friction coefficient,
 *          -rolling friction coefficient,
 *          -torsion friction coefficient,
 *          -loading stiffness,
 *          -unloading stiffness,
 *          -cohesion stiffness,
 *          -number of particles for computing internal structure (if needed),
 *          -energy ratio threshold,
 *          -velocity damping modulus,
 *          -restitution cohefficient.
 *
 *          After this the indentation of the file is written. The values that will be written are (in order):
 *          -elastic energy,
 *          -energy ratio threshold,
 *          -mean coordination number,
 *          -mean cluster radius,
 *          -solid fraction,
 *          -force modulus,
 *          -mean force on interaction,
 *          -minimum relative overlap,
 *          -mean relative overlap,
 *          -maximum realative overlap,
 *          -center of mass.
 */
void BaseCluster::makeCdatFile()
{
    std::ostringstream cdatName;
    cdatName << getName() << ".cdat";

    cdatFile_.open(cdatName.str(), std::ios::out);

    cdatFile_ << "CLUSTER DATA AND INFORMATION" << std::endl << std::endl;
    cdatFile_ << "position: " << position_ << std::endl;
    cdatFile_ << "collisionTimeOverTimeStep: " << getCollisionTimeOverTimeStep() << std::endl;
    cdatFile_ << "radiusParticle: " << std:: scientific << std::setprecision(2) << getRadiusParticle() << std::endl;
    cdatFile_ << "sizeDispersityParticle: " << std::defaultfloat << getSizeDispersityParticle() << std::endl;
    cdatFile_ << "densityParticle: " << std:: scientific << particleSpecies_ -> getDensity() << std::endl;
    cdatFile_ << "nParticles: " << std::defaultfloat << getNumberOfParticles() << std::endl;
    cdatFile_ << "idCluster: " << getClusterId() << std::endl;
    cdatFile_ << "slidingFrictionCoeff: " << particleSpecies_ -> getSlidingFrictionCoefficient() << std::endl;
    cdatFile_ << "rollingFrictionCoeff: " << particleSpecies_ -> getRollingFrictionCoefficient() << std::endl;
    cdatFile_ << "torsionFrictionCoeff: " << particleSpecies_ -> getTorsionFrictionCoefficient() << std::endl;
    // If constantRestitution(true) loading, unloading, and cohesion stiffness are multiplied by the mass of a particle (massParticle_) whose radius is radiusParticle_*2/(1+sizeDispersityParticle_),
    // which is the mass that should be used to compute collision time.
    cdatFile_ << "loadingStiffness: " << std::scientific << particleSpecies_ -> getLoadingStiffness()
                                            * (particleSpecies_->getConstantRestitution()?massParticle_:1) << std::endl;
    cdatFile_ << "unloadingStiffnessMax: " << particleSpecies_ -> getUnloadingStiffnessMax()
                                            * (particleSpecies_->getConstantRestitution()?massParticle_:1) << std::endl;
    cdatFile_ << "cohesionStiffness: " << particleSpecies_ -> getCohesionStiffness()
                                            * (particleSpecies_->getConstantRestitution()?massParticle_:1) << std::endl;
    cdatFile_ << "restitutionCoefficient: " << particleSpecies_ -> getRestitutionCoefficient(massParticle_) << std::endl;
    cdatFile_ << "collisionTime: " << std::scientific << std::setprecision(3) << particleSpecies_ -> getCollisionTime(massParticle_) << std::endl;
    if (getNumberOfInternalStructurePoints())
        cdatFile_ << "nInternalStructurePoints: " << getNumberOfInternalStructurePoints() << std::endl;
    cdatFile_ << "energyRatioTolerance: " << getEnergyRatioTolerance() << std::endl;
    cdatFile_ << "velocityDampingModulus: " << std::defaultfloat << getVelocityDampingModulus() << std::endl << std::endl;

    cdatFile_ << "progress" << std::setw(16) << "ElastEne" << std::setw(23) << "Ekin/ElastEne" << std::setw(18) << "coord_number" << std::setw(14) << "meanRadius"
             << std::setw(19) << "solidFraction" << std::setw(16) << "forceModulus" << std::setw(22) << "AveFOnOverl" << std::setw(15) << "dMin" << std::setw(17) << "dMean"
             << std::setw(14) << "dMax" << std::setw(24) << "Mass Centre" << std::endl;
}

/*!
 * \detail Creates the cluster overlap output file. In this file forces vs overlaps will be written. Together with
 *         this a ".gnuplot" file is created, which reads from this, ready to be loaded (see
 *         BaseCluster::makeGnuplotFile() ).
 */
void BaseCluster::makeOverlFile()
{
    std::ostringstream overlName;
    overlName << getName() << ".overl";

    overlFile_.open(overlName.str(), std::ios::out);
    overlFile_ << "Overlap Vs Normal Force" << std::endl;
}

/*! \details Tries to insert a particle in a spherical area centered around the initial site of the cluster;
 *          With an insertion fail counter and a while cycle this function tries to insert a particle in the domain:
 *          in order to do this, after computing random spherical coordinates, checks for interaction. If no interaction is detected
 *          the particle is inserted and returns true, if it fails 1000 times returns false and the computation is over.
 *          The spherical coordinates are rescaled over the radius by a factor cbrt(rand(0,1) ) and also over vertical angle thanks to
 *          acos(rand (-1, 1) ): this has been done in order to ensure the most possible spherical shape to the cluster, which with
 *          classic spherical coordinates was not achieved.
 * @param n: index of the particle being inserted.
 * @return bool: true if particle inserted, false if not.
 */
bool BaseCluster::particleInsertionSuccessful(int n) {

    int insertionFailCounter = 0;
    Mdouble rad, theta, phi;
    Vec3D particlePosition;
    SphericalParticle p0;

    // setup of particle properties and initial conditions (besides position)
    p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
    p0.setRadius(radii_[n]);
    p0.setSpecies(speciesHandler.getObject(0));
    p0.setGroupId(idCluster_);

    while (insertionFailCounter < 1000)
    {
        theta = constants::pi * random.getRandomNumber(0, 2.0);
        phi = acos(random.getRandomNumber(-1.0, 1.0));
        rad = p0.getRadius() + cbrt( random.getRandomNumber( 0, 1 ) ) * ( 0.5 * boxSize_ - 2.01 * p0.getRadius() );

        particlePosition.X = position_.X + rad * sin(phi) * cos(theta);
        particlePosition.Y = position_.Y + rad * sin(phi) * sin(theta);
        particlePosition.Z = position_.Z + rad * cos(phi);

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



/*
 * ----------------------------------------------------------------------
 *        FUNCTIONS: functions inside actionsAfterTimeStep
 * ----------------------------------------------------------------------
 */

/*!
 * \details This functions computes some important cluster information needed by the program. They are (in order):
 *          -mean coordination number,
 *          -center of mass,
 *          -mean cluster radius,
 *          -cluster radius used to compute solid fraction (different from mean cluster radius, see \details below),
 *          -solid fraction,
 *          -mean force acting on interaction,
 *          -maximum overlap,
 *          -mean overlap,
 *          -minimum overlap.
 */
void BaseCluster::makeDataAnalysis()
{
    // resetting counters and variables
    Mdouble solidVolumeInsideRadius = 0;
    Vec3D localMin;
    Vec3D localMax;
    localMin.setZero();
    localMax.setZero();
    centerOfMass_.setZero();
    //\brief vector in which it is saved the relative position of a particle from the center of mass
    Vec3D distanceFromCenterOfMass;
    //\brief distance from the center of mass of the furthest particle
    Mdouble furthestParticleDistance = 0;
    // number of particles whose distance d from the center of mass is d > furthestParticleDistance - radiusParticle_
    int counter = 0;
    Mdouble relativeOverlap = 0;

    meanClusterRadius_ = 0.0;
    meanCoordinationNumber_ = 0.0;
    maxRelativeOverlap_ = 0.0;
    meanRelativeOverlap_ = 0.0;
    minRelativeOverlap_ = 2 * 2 * sizeDispersityParticle_/(1+sizeDispersityParticle_);

    // loops over each particle to compute mean coordination number and center of mass.
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p) {

        meanCoordinationNumber_ += ((*p)->getInteractions()).size();

        centerOfMass_ += ((*p)->getVolume()) * ((*p)->getPosition());
    }

    meanCoordinationNumber_ /= particleHandler.getSize();
    centerOfMass_ /= totalParticleVolume_;


    // loops over each particle to compute the furthest particle from the center of mass.
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p) {

        distanceFromCenterOfMass = (*p)->getPosition() - centerOfMass_;

        if (distanceFromCenterOfMass.getLength() > furthestParticleDistance)
            furthestParticleDistance = distanceFromCenterOfMass.getLength();

    }

    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p) {

        distanceFromCenterOfMass = (*p)->getPosition() - centerOfMass_;

        if (distanceFromCenterOfMass.getLength() > furthestParticleDistance - radiusParticle_) {
            meanClusterRadius_ += distanceFromCenterOfMass.getLength();
            counter++;
        }

    }

    meanClusterRadius_ /= counter;
    meanClusterRadius_ += radiusParticle_;


    /*
     * \details This is the radius used to compute solid fraction: it is smaller than the meanClusterRadius.
     */
    radiusForSolidFraction_ = meanClusterRadius_ - 1 * radiusParticle_ ;

    /*
     * \details With a for cycle the volume of the particles inside radiusForSolidFraction is computed and after this the value of solid
     *          fraction is calculated. This value is less precise as the maximum penetration depth increases and more precise as the
     *          number of particle increases.
     */
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p)
    {
        distanceFromCenterOfMass = (*p) -> getPosition() - centerOfMass_;

        if( distanceFromCenterOfMass.getLength() < radiusForSolidFraction_ )
            solidVolumeInsideRadius += (*p) -> getVolume();
    }

    solidFraction_ = 3 * solidVolumeInsideRadius / ( 4 * constants::pi * pow(radiusForSolidFraction_, 3) );

    // loops over every interaction to compute mean force acting on interaction, maximum, mean and minimum relative particle overlap.
    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
    {

        /*
         * \details the relative overlap is computed as an average of the relative overlap on the two particles.
         *          rO = ( O/R1 + O/R2 ) / 2.
         */
        relativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                                    ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
        relativeOverlap /= 2;
        meanRelativeOverlap_ += relativeOverlap;
        if (relativeOverlap > maxRelativeOverlap_)
            maxRelativeOverlap_ = relativeOverlap;

        if (relativeOverlap < minRelativeOverlap_)
            minRelativeOverlap_ = relativeOverlap;
    }
    meanRelativeOverlap_ /= interactionHandler.getSize();
}

/*!
 * \details This writes on the cluster data output file. Written values are (in order):
 *          -elastic energy,
 *          -energy ratio threshold,
 *          -mean coordination number,
 *          -mean cluster radius,
 *          -solid fraction,
 *          -force modulus,
 *          -mean force on interaction,
 *          -minimum relative overlap,
 *          -mean relative overlap,
 *          -maximum realative overlap,
 *          -center of mass.
 */
void BaseCluster::writeToCdatFile()
{
    switch (stage_)
    {
        case 1:
            cdatFile_ << "C, " << std::fixed << std::setprecision(0) << std::setw(4) << 100 * (getTime() - t0_) / forceTuningDuration_ << "%: ";
            break;

        case 2:
            cdatFile_ << "D, " << std::fixed << std::setprecision(0) << std::setw(4) << 100 * (getTime() - t0_) / forceTuningDuration_ << "%: ";
            break;

        case 3:
            cdatFile_ << "D-energy: ";
            break;

        default:
            cdatFile_ << "Final  v: ";
    }

    cdatFile_ <<
            std::scientific <<
            std::setprecision(2)<<
            std::setw(14) <<
            getElasticEnergy() <<
            std::setw(18) <<
            getKineticEnergy()/getElasticEnergy() <<
            std::fixed <<
            std::setprecision(2) <<
            std::setw(15) <<
            meanCoordinationNumber_ <<
            std::scientific <<
            std::setw(20) <<
            meanClusterRadius_ <<
            std::fixed <<
            std::setprecision(2) <<
            std::setw(13) <<
            solidFraction_ <<
            std::scientific <<
            std::setprecision(3) <<
            std::setw(24) <<
            forceModulus_ <<
            std::setw(22) <<
            std::fixed <<
            std::setprecision(5) <<
            std::setw(18) <<
            minRelativeOverlap_ <<
            std::setw(16) <<
            meanRelativeOverlap_ <<
            std::setw(15) <<
            maxRelativeOverlap_ <<
            std::scientific <<
            std::setprecision(2) <<
            std::setw(19) <<
            centerOfMass_.X <<
            std::setw(12) <<
            centerOfMass_.Y <<
            std::setw(12) <<
            centerOfMass_.Z <<
            "          " <<
            std::endl;
}

/*!
 * \details This writes on the cluster overlap output file. At each output time and  for each
 * interaction force vs overlap is written down.
 */
void BaseCluster::writeToOverlFile()
{
    //\brief force acting on overlap.
    Mdouble forceOnOverlap = 0;
    Mdouble relativeOverlap = 0;
    switch (stage_)
    {
        case 1:
            overlFile_ << "C, " << std::fixed << std::setprecision(0) << std::setw(4) << 100 * (getTime() - t0_) / forceTuningDuration_ << "%: ";
            break;

        case 2:
            overlFile_ << "D, " << std::fixed << std::setprecision(0) << std::setw(4) << 100 * (getTime() - t0_) / forceTuningDuration_ << "%: ";
            break;

        case 3:
            overlFile_ << "D energy: ";
            break;

        default:
            overlFile_ << "Final  v: ";
    }

    for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
    {
        forceOnOverlap = ((*i) -> getForce()).getLength();
        relativeOverlap = ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                            ((*i) -> getOverlap())/(particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
        relativeOverlap = relativeOverlap / 2;
        overlFile_ << std::setprecision(2) << std::scientific << std::setw(18) << forceOnOverlap
                    << std::defaultfloat << std::fixed << std::setprecision(4) << std::setw(9) << relativeOverlap;
    }

    overlFile_ << "     " << std::endl;
}

/*!
 * \details This applies force on each particle. The force applied is proportional to the distance from force center
 *          (which is position_) and normalized by a ( 1/(2r) + 1/norm ) factor: this way even if a particle is
 *          very close to the force center (so -forceModulus * distanceFromForceCenter / (2 * r) ~ 0 )
 *          a force F = -forceModulus * distanceFromForceCenter / norm is applied, which is the minimum guaranteed.
 *          (r is radiusParticle and norm is the 2-norm of distanceFromForceCenter vector).
 */
void BaseCluster::applyCentralForce()
{
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p)
    {
        //\brief distance from the center of forces (which is position_).
        Vec3D distanceFromForceCenter = (*p) -> getPosition() - position_;

        //\brief norm of distanceFromForceCenter vector
        Mdouble norm = distanceFromForceCenter.getLength();

        (*p) -> addForce( -forceModulus_ * distanceFromForceCenter * (2*radiusParticle_ + norm) / (2*radiusParticle_*norm) );

        }
}

/*!
 * \details This increases the value of forceModulus (stage = 1).
 *          ForceModulus varies from 0 to maximumForceModulus linearly with time.
 */
void BaseCluster::increaseForce()
{
    forceModulus_ = maximumForceModulus_*(getTime() - t0_)/forceTuningDuration_;
}

/*!
 * \details This damps values of each particle velocity (stage = 1, stage = 2, stage = 3).
 *          Damping is done by a factor velocityDampingModulus.
 */
void BaseCluster::dampVelocities()
{
    for (auto p = particleHandler.begin(); p != particleHandler.end(); ++p)
    {
        (*p) -> setVelocity(velocityDampingModulus_*( (*p) -> getVelocity() ));
    }
}

/*!
 * \details This linearly decreases values of forceModulus (stage = 2). forceModulus varies from maximumForceModulus to 0 linearly with time.
 *          Actually last value reached of forceModulus with this is maximumForceModulus * timeStep / forceTuningDuration because this
 *          process is discrete: starting from this value forceModulus will then be damped in stage = 3.
 */
void BaseCluster::decreaseForce()
{
    forceModulus_ = maximumForceModulus_ * (1 - (getTime() - t0_)/forceTuningDuration_);
}

/*!
 * \details This damps values of forceModulus (stage = 3).  Damping is done by a factor forceDampingModulus_.
 */
void BaseCluster::dampForce()
{
    forceModulus_ *= forceDampingModulus_;
}

/*!
 * \details This calculates the adjacency matrix of the cluster. Firstly a  NxN matrix is created and filled with zeros.
 *          with a for cycle on the interactions then, if there's a contact between particle i and j adjacencyMatrix(i,j) = adjacencyMatrix(j,i) = 1.
 *          (N is the number of particles).
 */
void BaseCluster::createAdjacencyMatrix()
{
    for (int i = 0; i < particleHandler.getSize(); i++)
    {
        std::vector<int> temporaryRowVector;
        temporaryRowVector.reserve(particleHandler.getSize());

        for (int j = 0; j < particleHandler.getSize(); j++)
            temporaryRowVector.push_back(0);

        adjacencyMatrix_.push_back(temporaryRowVector);
    }

    for (auto i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
    {
        adjacencyMatrix_[(*i) -> getP() -> getIndex()][(*i) -> getI() -> getIndex()] = 1;
        adjacencyMatrix_[(*i) -> getI() -> getIndex()][(*i) -> getP() -> getIndex()] = 1;
    }
}

/*!
 * \details This creates the adjacency matrix file.
 */
void BaseCluster::makeAmatFile()
{
    std::ostringstream amatName;
    amatName << getName() << ".amat";

    amatFile_.open(amatName.str(), std::ios::out);
    amatFile_ << "ADJACENCY MATRIX" << std::endl << std::endl;
}

/*!
 * \details This writes on the adjacency matrix file. With two nested for cycles adjacencyMatrix is written down.
 *          Finally number of intra-cluster bonds and mean coordination number are also written.
 *          After this the file is closed.
 */
void BaseCluster::writeAmatFile()
{
    for(int i=0; i < particleHandler.getSize(); i++)
    {
        for(int j=0; j < particleHandler.getSize(); j++)
        {
            amatFile_ << adjacencyMatrix_[i][j] << "  ";
        }
        amatFile_ << std::endl;
    }
    nIntraClusterBonds_ = interactionHandler.getSize();
    amatFile_ << std::endl;
    amatFile_ << "THE TOTAL NUMBER OF INTRACLUSTER BONDS IS: " << nIntraClusterBonds_ << std::endl;
    amatFile_ << "THE MEAN COORDINATION NUMBER IS: " << meanCoordinationNumber_ << std::endl;

    amatFile_.close();
}

/*!
 * \details This computes the internal structure  and solid fraction of the cluster.
 *          A total number of particles (equal to nInternalStructurePoints_) is tried to be inserted (with spherical
 * 			coordinates	identical to the ones used for inserting the actual particles) inside a sphere
 *          having radius radiusForSolidFraction_. If no interaction is found, nothing happens, otherwise the counuter
 *          nPointsInsideComponentsForMCTest is incremented and the corresponding position is written down in the
 *          internal structure file.
 */
void BaseCluster::computeInternalStructure()
{



    Vec3D mcPoint;
    SphericalParticle p0;
    Mdouble fictitiousGridPointRadiusRatio = 1.0e-5;
    p0.setSpecies(speciesHandler.getObject(0));
    p0.setRadius(radiusParticle_*fictitiousGridPointRadiusRatio);
    p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
    int nMonteCarloSamplingPoints = nInternalStructurePoints_;
    Mdouble nPointsInsideComponentsForMCTest = 0;

    makeIntenalStructureFile();

    for (int i = 0; i < nMonteCarloSamplingPoints; ++i)
    {

        Mdouble theta = constants::pi * random.getRandomNumber(0, 2.0);
        Mdouble phi = acos(random.getRandomNumber(-1.0, 1.0));
        Mdouble rad = radiusForSolidFraction_*cbrt( random.getRandomNumber( 0, 1 ) );

        mcPoint.X = rad * sin(phi) * cos(theta);
        mcPoint.Y = rad * sin(phi) * sin(theta);
        mcPoint.Z = rad * cos(phi);
        mcPoint += centerOfMass_;

        p0.setPosition(mcPoint);

        if (!checkParticleForInteraction(p0)) // collision -> the counter goes to the mass fraction
        {
            nPointsInsideComponentsForMCTest++;
            intStructFile_ << std::scientific << std::setprecision(5) << std::setw(12) << mcPoint.X
                           << std::setw(13) << mcPoint.Y << std::setw(13) << mcPoint.Z << std::setw(6) << 0 << std::endl;
        }

    }

    // Solid fraction (accordance between theoretical values and penetration depth max).
    // It is very important to notice that this value is accurate only if sliding friction is set to 0.5 and relative
    // tangential stiffness is set to 0.3 while creating the cluster. Different values do not guarantee accuracy.
    solidFractionIntStruct_ = nPointsInsideComponentsForMCTest/nMonteCarloSamplingPoints;
    Mdouble theoVal = 0.58 + 3*pow(0.58,2)*particleSpecies_->getPenetrationDepthMax();
    Mdouble diff = fabs(theoVal-solidFractionIntStruct_);
    Mdouble accordance = (theoVal - diff)/theoVal;
    // Solid fraction (accordance between theoretical values and average overlap).
    // It is very important to notice that this value is accurate only if sliding friction is set to 0.5 and relative
    // tangential stiffness is set to 0.3 while creating the cluster. Different values do not guarantee accuracy.
    solidFractionIntStruct_ = nPointsInsideComponentsForMCTest/nMonteCarloSamplingPoints;
    Mdouble theoValAvOverl = 0.58 + 3*pow(0.58,2)*meanRelativeOverlap_;
    Mdouble diffAvOverl = fabs(theoValAvOverl-solidFractionIntStruct_);
    Mdouble accordanceAvOverl = (theoValAvOverl - diffAvOverl)/theoValAvOverl;

    intStructFile_ << "n_points_inside_boundary: " << std::scientific << nMonteCarloSamplingPoints << std::endl;
    intStructFile_ << "n_points_inside_components: " << nPointsInsideComponentsForMCTest << std::endl;
    intStructFile_ << "solidFractionIntStruct_: " << std::fixed << std::setprecision(6) << solidFractionIntStruct_
                   << ", accordance with theoretical values: " << 100*accordance << "%." << std::endl
                   << "Accordance with average overlap: " << 100*accordanceAvOverl << "%." << std::endl
                   << "It is very important to notice that this formula is accurate only if sliding friction" << std::endl
                   << "is set to 0.5 and relative tangential stiffness is set to 0.3 while creating the cluster." << std::endl
                   << "Different values do not guarantee accuracy." << std::endl << std::endl;

    /*
     * computeInternalStructure output is set to VERBOSE in order not to have too much output. If the user needs it,
     * it is enough to set it to INFO.
     */

    std::ostringstream printResults;
    printResults << "n_points_inside_boundary: " << std::scientific << nMonteCarloSamplingPoints << std::endl;
    printResults << "n_points_inside_components: " << nPointsInsideComponentsForMCTest << std::endl;
    printResults << "solidFractionIntStruct_: " << std::fixed <<  std::setprecision(6) << solidFractionIntStruct_
                 << ", accordance with theoretical values: " << 100*accordance << "%." << std::endl
                 << "Accordance with average overlap: " << 100*accordanceAvOverl << "%." << std::endl
                 << "It is very important to notice that this formula is accurate only if sliding friction" << std::endl
                 << "is set to 0.5 and relative tangential stiffness is set to 0.3 while creating the cluster." << std::endl
                 << "Different values do not guarantee accuracy." << std::endl << std::endl;
    logger(VERBOSE, printResults.str());
}

/*!
 * \details This creates the gnuplot file needed for printing force vs overlaps values.
 *          After setting output tipe, title, labels and grid, all columns needed for overlap printing are written after the plot command.
 *          Only overlaps present at the end of the process will be printed. Loading this file with gnuplot a jpeg image is obtained
 *          showing all forces vs overlaps.
 */
void BaseCluster::makeGnuplotFile()
{
    std::ostringstream gnuplotname;
    gnuplotname << getName() << ".gnuplot";

    gnuplotFile_.open(gnuplotname.str(), std::ios::out);
    gnuplotFile_ << "set terminal jpeg" << std::endl;
    gnuplotFile_ << "set output \"" << getName() << "_overlaps" << ".jpeg\"" << std::endl;
    std::string titleLine=R"(set title "Overlap Vs Force")";// font ",14"
    std::string xLabel=R"(set xlabel "Overlap")";
    std::string yLabel=R"(set ylabel "Force")";
    gnuplotFile_ << titleLine << std::endl;
    gnuplotFile_ << xLabel << std::endl;
    gnuplotFile_ << yLabel << std::endl;
    gnuplotFile_ << "set grid" << std::endl;
    gnuplotFile_ << "plot ";
    for (int i = 0; i < interactionHandler.getSize(); ++i)
    {
        gnuplotFile_ << "\"" << getName() << ".overl" << "\"" << " using " << 2*i+4 << ":" << 2*i+3 << " title \"\" with lines lt 1  dashtype 2, ";
    }

    gnuplotFile_.close();

}

/*!
 * \details This creates the file needed for writing down datas from computeInternalStructure().
 *          In the first row nInternalStructurePoints is printed: this number is the total number of points
 *          tried to be inserted (and so is grater than the number of lines of the file, which corresponds to
 *          the number of points that gave interaction).
 */
void BaseCluster::makeIntenalStructureFile()
{
    std::ostringstream intStructName;
    intStructName << getName() << ".struct";

    intStructFile_.open(intStructName.str(), std::ios::out);
    intStructFile_ << "Number of Montecarlo points: " << nInternalStructurePoints_ << std::endl;
}
