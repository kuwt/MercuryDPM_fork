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

#include "InsertionBoundary.h"
#include "DPMBase.h"
#include "Particles/BaseParticle.h"
#include<iostream>

#ifdef MERCURY_USE_MPI
#include "MpiDataClass.h"
#endif

/*!
 * \details Default constructor, sets all data members to 0, nullptr or default.
 */
InsertionBoundary::InsertionBoundary()
{
    numberOfParticlesInserted_ = 0;
    massInserted_ = 0;
    volumeInserted_ = 0;
    particleToCopy_ = nullptr;
    maxFailed_ = 0;
    isActivated_ = true;
    volumeFlowRate_ = constants::inf;
    initialVolume_ = 0;
    samplingInterval_ = 0;
    checkParticleForInteraction_ = true;
    radMin_ = 0;
    radMax_ = 0;
    isManuallyInserting_ = false;
}

/*!
 * \details Copy constructor
 */
InsertionBoundary::InsertionBoundary(const InsertionBoundary& other)
        : BaseBoundary(other)
{
    numberOfParticlesInserted_ = other.numberOfParticlesInserted_;
    massInserted_ = other.massInserted_;
    volumeInserted_ = other.volumeInserted_;
    maxFailed_ = other.maxFailed_;
    isActivated_ = other.isActivated_;
    volumeFlowRate_ = other.volumeFlowRate_;
    initialVolume_ = other.initialVolume_;
    samplingInterval_ = other.samplingInterval_;
    variableCumulativeVolumeFlowRate_ = other.variableCumulativeVolumeFlowRate_;
    checkParticleForInteraction_ = other.checkParticleForInteraction_;
    distribution_ = other.distribution_;
    particleSizeDistribution_ = other.particleSizeDistribution_;
    radMin_ = other.radMin_;
    radMax_ = other.radMax_;
    isManuallyInserting_ = other.isManuallyInserting_;
    
    if (other.particleToCopy_ != nullptr)
    {
        particleToCopy_ = other.particleToCopy_->copy();
    }
    else
    {
        particleToCopy_ = nullptr;
    }
}

/*!
 * \details Destructor that deletes the BaseParticle that is copied and inserted
 *          at every insertion.
 */
InsertionBoundary::~InsertionBoundary()
{
    delete particleToCopy_;
}

/*!
 * \details Sets the particle that will be inserted and the maximum number of 
 * times for which insertion may fail.
 * \param[in] particleToCopy  Particle that will be copied and inserted in the domain
 * \param[in] maxFailed       Number of times that the wall may fail to insert a particle
 */
void InsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed)
{
    if (particleToCopy != nullptr)
        particleToCopy_ = particleToCopy->copy();
    maxFailed_ = maxFailed;
}


/*!
 * \details The default behaviour will be to return particleToCopy_, but this
 * can be overridden by the children (to get size or species dispersity).
 */
BaseParticle* InsertionBoundary::generateParticle(RNG& random)
{
    BaseParticle* P = getParticleToCopy()->copy();
    if (particleSizeDistribution_.getParticleSizeDistribution().empty())
    {
        switch (distribution_)
        {
            default:
                // default is a uniform distribution between radMin_ and radMax_
                P->setRadius(random.getRandomNumber(radMin_, radMax_));
                break;
            case Distribution::Normal_1_5:
                // normal distribution
                /// \todo TP: can we add a variable seed here? Let the user choose a seed and for reproducible results you always take the same seed.
                static std::mt19937 gen(0);
                Mdouble particleRadius = 0.5 * (radMax_ + radMin_);
                Mdouble polydispersity = 0.5 * (radMax_ - radMin_);
                static std::normal_distribution<> d(particleRadius, polydispersity);
                static const Mdouble radiusMin = particleRadius - 1.5 * polydispersity;
                static const Mdouble radiusMax = particleRadius + 1.5 * polydispersity;
                Mdouble radius = d(gen);
                while (radius > radiusMax || radius < radiusMin) radius = d(gen);
                P->setRadius(radius);
                break;
        }
    }
        // manual insertion routine to insert PSDs as accurate as possible into a given volume
    else if (isManuallyInserting_)
    {
        logger.assert(initialVolume_ > 0.0, "Use setInitialVolume to define the particle insertion volume");
        Mdouble radius;
        // getVolumeFlowRate() * time + initialVolume_ - volumeInserted_ lead to more inaccurate results, therfore
        // -volumeInserted was removed.
        radius = particleSizeDistribution_.insertManuallyByVolume(getVolumeFlowRate() *
                                                                  getHandler()->getDPMBase()->getTime() +
                                                                  initialVolume_);
        P->setRadius(radius);
    }
    else
    {
        Mdouble radius;
        radius = particleSizeDistribution_.drawSample();
        P->setRadius(radius);
    }
    return P;
}

bool InsertionBoundary::insertParticle(Mdouble time)
{
    // check if the flow rate limit has been reached
    if (variableCumulativeVolumeFlowRate_.empty())
    {
        return volumeInserted_ < getVolumeFlowRate() * time + initialVolume_;
    }
    else
    {
        const Mdouble iMax = (Mdouble) variableCumulativeVolumeFlowRate_.size() - 2;
        const Mdouble i = std::min(time / samplingInterval_, iMax);
        if (i == iMax)
        {
            static unsigned count = 0;
            if (count == 0) logger(WARN, "Reached end of volume flowrate function");
            ++count;
        }
        const size_t id = i;
        const Mdouble allowedVolume = variableCumulativeVolumeFlowRate_[id] +
                                      (variableCumulativeVolumeFlowRate_[id + 1] -
                                       variableCumulativeVolumeFlowRate_[id]) * (i - id);
        return volumeInserted_ < allowedVolume;
    }
}

/*!
 * \details Is used to fill the insides of the boundary with particles until 
 * it is filled up. 
 * \param[in,out] md    the problem's DPMBase object
 * \todo rename to something like "insertUntilMaxFailed"?
 */
void InsertionBoundary::checkBoundaryBeforeTimeStep(DPMBase* md)
{
    logger(VERBOSE, "In InsertionBoundary::checkBoundaryBeforeTimeStep\n");
    
    if (!isActivated_)
        return;
    
    /* Each timestep, the InsertionBoundary attempts to fill up a region with
     * particles.
     *
     * It first calls generateParticle() to get a randomised particle, subject
     * to a specified distribution over sizes and species. (The basic class
     * supports size dispersity only but PolydisperseInsertionBoundary will
     * support species dispersity.) 
     * Then it repeatedly calls placeParticle(), which gives the particle a
     * random location (and possibly velocity) in a specified,
     * geometry-dependent bound. Each time, it checks whether the new particle
     * would have an interaction with another particle or a wall. 
     *
     * If it manages to do that within maxFailed_ tries, then:
     *   * the new particle is inserted, 
     *   * the failure counter is reset, and
     * the processes is repeated with a new generateParticle().
     *
     * Otherwise, the processes terminates for this timestep.
     * */
    
    // Keep count of how many successive times we have failed to place a new
    // particle. 
    unsigned int failed = 0;
    while (failed <= maxFailed_ && insertParticle(md->getNextTime())) // 'generating' loop
    {
        /* Generate random *intrinsic* properties for the new particle. */
        logger(VERBOSE, "about to call generateParticle\n");
    
    
    
    
        //! HERE
    
    
    
    
        auto p0 = generateParticle(md->random);
        // Important for particle generation with a particle size distribution as it generates a particle with zero
        // radius. If a particle is not allowed to be inserted by the PSD criteria it will generate a particle with
        // zero diameter. This if statement prevents inserting particles with zero radius, which would else be a problem.
        /// \todo create a definition for zero-size particles. Right now the zero-size particle is only used as a stop criterion for the manual PSD insertion
        if (p0->getRadius() == 0)
        {
            logger(VERBOSE, "The PSD for the specified volume is fully set");
            // free up memory space
            delete p0;
            // out of the 'placing' loop
            failed = maxFailed_ + 1;
            continue;
        }
        logger(VERBOSE, "generated a particle with intrinsics %", p0);
    
        while (true) // 'placing' loop
        {
            /* Generate extrinsic properties (position and velocity) for this
             * new particle. */
        
        
        
            //! HERE
    
    
    
    
    
            placeParticle(p0, md->random);
            logger(VERBOSE, "attempting to place particle at %, vel %", p0->getPosition(), p0->getVelocity());

#ifdef MERCURY_USE_MPI
            /* Communicate the new particle's properties by setHandler (note
             * that this doesn't actually add the particle to the handler). */
            if (NUMBER_OF_PROCESSORS > 1)
            {
                MPIParticle particle;
                // 	//Every domain generates a particle (to get the species right etc)

                //Send particle data from root to other processors to sync the particle properties
                if (PROCESSOR_ID == 0)
                {
                    particle.copyDataFromParticleToMPIParticle(p0);
                }

                MPIContainer::Instance().broadcast(&particle,MercuryMPIType::PARTICLE);

                //Process the received data
                if (PROCESSOR_ID != 0)
                {
                    copyDataFromMPIParticleToParticle(&particle, p0, &(md->particleHandler));
                }
            }
#endif
            p0->setHandler(&md->particleHandler);
            /* Check whether the particle has any interactions. */
            if (!checkParticleForInteraction_ || md->checkParticleForInteraction(*p0))
            {
                //Note: in parallel only one of the domains will actually add the particle
                auto p = md->particleHandler.copyAndAddObject(p0);
                failed = 0;
    
                ++numberOfParticlesInserted_;
                const double volume = p0->getVolume();
                volumeInserted_ += volume;
                massInserted_ += p0->getSpecies()->getDensity() * volume;
                logger(VERBOSE, "successfully placed a particle %, with position: % after % fails.", p0,
                       p0->getPosition(), failed);
                
                /* JMFT: The generateParticle() routine allocates memory, so we should
                 * free it here. (Don't worry, the particle will have been copied to the
                 * particleHandler by this point iff we want it.) */
                delete p0;
    
                break; // out of the 'placing' loop
            }
            else
            {
                failed++;
                logger(VERBOSE, "failed to place a particle; have failed % times", failed);
            }
        
            if (failed > maxFailed_)
            {
                logger(VERBOSE, "failed too many times; giving up");
                break; // out of the 'placing' loop (and will leave the 'generating' loop too
            }
        }
        logger(VERBOSE, "failed % times, so breaking out of InsertionBoundary loop for this timestep.", failed);
    }
    // logger(INFO, "volumeInserted_ = %", volumeInserted_);
}

/*!
 * \details Returns the number of particles inserted in the boundary
 * \return  the number of particles inserted
 */
unsigned int InsertionBoundary::getNumberOfParticlesInserted() const
{
    return numberOfParticlesInserted_;
}

Mdouble InsertionBoundary::getMassOfParticlesInserted() const
{
    return massInserted_;
}

Mdouble InsertionBoundary::getVolumeOfParticlesInserted() const
{
    return volumeInserted_;
}

/*! 
 * \details reset() does not activate or deactivate the InsertionBoundary.
 */
void InsertionBoundary::reset()
{
    numberOfParticlesInserted_ = 0;
    massInserted_ = 0;
    volumeInserted_ = 0;
}

void InsertionBoundary::activate()
{
    isActivated_ = true;
}

void InsertionBoundary::deactivate()
{
    isActivated_ = false;
}

bool InsertionBoundary::isActivated()
{
    return isActivated_;
}

/*!
 * \details Sets the maximum number of times InsertionBoundary::checkBoundaryBeforeTimeStep()
 * may try to insert a particle and fail, before the insertion of particles stops.
 * \param[in] maxFailed     the maximum number of particle insertion trials
 */
void InsertionBoundary::setMaxFailed(unsigned int maxFailed)
{
    maxFailed_ = maxFailed;
}

/*!
 * \details Return maxFailed_ (see InsertionBoundary::setMaxFailed).
 * \return the maximum number of particle insertion trials
 */
unsigned int InsertionBoundary::getMaxFailed() const
{
    return maxFailed_;
}

/*!
 * \details Sets the pointer to the particle, copies of which are inserted 
 * \param[in] particleToCopy    pointer to the particle to be inserted
 */
void InsertionBoundary::setParticleToCopy(BaseParticle* particleToCopy)
{
    if (particleToCopy != nullptr)
        particleToCopy_ = particleToCopy->copy();
    else
        logger(ERROR, "Setting particleToCopy to be a null pointer?");
}

/*!
 * \details returns pointer to the particle copies of which are to be inserted
 */
BaseParticle* InsertionBoundary::getParticleToCopy() const
{
    ///\todo make this a debug check
    //if (particleToCopy_==nullptr)
    //    std::cerr << "Error: particleToCopy not set" << std::endl;
    return particleToCopy_;
}

/*!
 * \details reads the boundary's id_ and maxFailed_ from the given istream
 * \param[in,out] is    stream the data members are read from
 */
void InsertionBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy, type;
    is >> dummy >> maxFailed_ >> dummy;
    if (dummy == "volumeFlowRate")
        is >> volumeFlowRate_ >> dummy;
    is >> massInserted_;
    is >> dummy >> volumeInserted_;
    is >> dummy >> numberOfParticlesInserted_;
    is >> dummy >> isActivated_;
    if (!particleSizeDistribution_.getParticleSizeDistribution().empty())
    {
        for (auto p : particleSizeDistribution_.getParticleSizeDistribution())
        {
            is >> dummy >> p.radius;
            is >> dummy >> p.probability;
        }
    }
    else
    {
        is >> dummy >> distribution_;
    }
    ///\todo make theses reads non-optional
    helpers::readOptionalVariable(is, "checkParticleForInteraction", checkParticleForInteraction_);
    helpers::readOptionalVariable(is, "initialVolume", initialVolume_);
    if (helpers::readOptionalVariable(is, "samplingInterval", samplingInterval_))
    {
        size_t n;
        Mdouble flowRate;
        is >> dummy >> n;
        //variableCumulativeVolumeFlowRate_.clear();
        variableCumulativeVolumeFlowRate_.reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
            is >> flowRate;
            variableCumulativeVolumeFlowRate_.push_back(flowRate);
        }
    }
    is >> dummy;
    if (dummy != "noParticleToCopy")
    {
        delete particleToCopy_;
        particleToCopy_ = getHandler()->getDPMBase()->particleHandler.readAndCreateObject(is);
        // The .restart file records the index of the particle's species, but
        // doesn't record the pointer, i.e. the memory address of the species within
        // the speciesHandler. The latter needs to be reset now.
        particleToCopy_->setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(
                particleToCopy_->getIndSpecies()
        ));
    }
}

/*!
 * \details adds the boundary's id_ and maxFailed_ to the given ostream
 * \param[in,out] is    stream the data members are to be added to
 */
void InsertionBoundary::write(std::ostream& os) const
{
    logger(VERBOSE, "In InsertionBoundary::write\n");
    BaseBoundary::write(os);
    os << " maxFailed " << maxFailed_;
    if (std::isfinite(volumeFlowRate_))
        os << " volumeFlowRate " << volumeFlowRate_;
    os << " massInserted " << massInserted_;
    os << " volumeInserted " << volumeInserted_;
    os << " numberOfParticlesInserted " << numberOfParticlesInserted_;
    os << " isActivated " << isActivated_;
    if (!particleSizeDistribution_.getParticleSizeDistribution().empty())
    {
        os << " psd " << particleSizeDistribution_.getParticleSizeDistribution().size();
        for (auto p : particleSizeDistribution_.getParticleSizeDistribution())
        {
            os << " " << p.radius
               << " " << p.probability;
        }
    }
    else
    {
        os << " distribution " << distribution_;
    }
    if (checkParticleForInteraction_ == false)
    {
        os << " checkParticleForInteraction " << checkParticleForInteraction_;
    }
    os << " initialVolume " << initialVolume_;
    os << " samplingInterval " << samplingInterval_;
    os << " variableCumulativeVolumeFlowRate " << variableCumulativeVolumeFlowRate_.size();
    for (const auto flowRate : variableCumulativeVolumeFlowRate_)
    {
        os << ' ' << flowRate;
    }
    if (particleToCopy_ == nullptr)
        os << " noParticleToCopy";
    else
    {
        os << " particleToCopy " << *particleToCopy_;
    }
}

Mdouble InsertionBoundary::getVolumeFlowRate() const
{
    return volumeFlowRate_;
}

void InsertionBoundary::setVolumeFlowRate(Mdouble volumeFlowRate)
{
    volumeFlowRate_ = volumeFlowRate;
}

Mdouble InsertionBoundary::getInitialVolume() const
{
    return initialVolume_;
}

void InsertionBoundary::setInitialVolume(Mdouble initialVolume)
{
    initialVolume_ = initialVolume;
    if (!std::isfinite(volumeFlowRate_)) volumeFlowRate_ = 0;
}

void InsertionBoundary::setVariableVolumeFlowRate(const std::vector<Mdouble>& variableCumulativeVolumeFlowRate,
                                                  Mdouble samplingInterval)
{
    logger.assert(samplingInterval > 0, "sampling interval needs to be positive");
    const Mdouble endTime = variableCumulativeVolumeFlowRate.size() * samplingInterval;
    logger(INFO, "variable flowrate is defined up to %", endTime);
    logger.assert_always(getHandler()->getDPMBase()->getTimeMax() < endTime,
                         "variable flowrate is defined up to %, but tMax is set to %", endTime,
                         getHandler()->getDPMBase()->getTimeMax());
    variableCumulativeVolumeFlowRate_ = variableCumulativeVolumeFlowRate;
    samplingInterval_ = samplingInterval;
}

/*!
 * \brief Sets the range of particle radii that may be generated to a PSD defined by the user.
 */
/// \todo TP: Consider std::move instead of a set function. This would result in a speedup + no one has to set the PSD for insertionBoundaries again as it is set when a PSD is inserted.
void InsertionBoundary::setPSD(PSD psd)
{
    particleSizeDistribution_ = psd;
}

/*!
 * \details gets the user defined particle size distribution
 * \return a vector of the PSD::RadiusAndProbability class containing the particle size distribution
 */
PSD InsertionBoundary::getPSD()
{
    return particleSizeDistribution_;
}

/*!
 * \details sets the isManuallyInserting_ to TRUE, resulting in a top-down class-by-class insertion routine to insert PSDs
 * as accurate as possible.
 */
void InsertionBoundary::setManualInsertion(bool isManuallyInserting)
{
    isManuallyInserting_ = isManuallyInserting;
}

void InsertionBoundary::setDistribution(InsertionBoundary::Distribution distribution)
{
    distribution_ = distribution;
}

/*!
 * \brief Gets the range of particle radii that may be generated.
 */
InsertionBoundary::Distribution InsertionBoundary::getDistribution()
{
    return distribution_;
}

/*!
 * \details write distribution type variables to file.
 */
std::ostream& operator<<(std::ostream& os, InsertionBoundary::Distribution type)
{
    os << static_cast<unsigned>(type);
    return os;
}

/*!
 * \details write Distribution type variables from file.
 */
std::istream& operator>>(std::istream& is, InsertionBoundary::Distribution& type)
{
    unsigned uType;
    is >> uType;
    type = static_cast<InsertionBoundary::Distribution>(uType);
    return is;
}