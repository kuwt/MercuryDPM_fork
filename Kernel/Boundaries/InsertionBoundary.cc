//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
 * \details Default constructor, sets all data members to 0 or nullptr.
 */
InsertionBoundary::InsertionBoundary()
{   
    numberOfParticlesInserted_ = 0;
    massInserted_ = 0;
    volumeInserted_ = 0;
    particleToCopy_ = nullptr;
    maxFailed_ = 0;
    isActivated_ = true;
    volumeFlowRate_ = inf;
}

/*!
 * \details Copy constructor
 */
InsertionBoundary::InsertionBoundary(const InsertionBoundary& other) 
{
    numberOfParticlesInserted_ = other.numberOfParticlesInserted_;
    massInserted_ = other.massInserted_;
    volumeInserted_ = other.volumeInserted_;
    maxFailed_ = other.maxFailed_;
    isActivated_ = other.isActivated_;
    volumeFlowRate_ = other.volumeFlowRate_;

    if (other.particleToCopy_!=nullptr) {
        particleToCopy_ = other.particleToCopy_->copy();
    } else {
        particleToCopy_ = nullptr;
    }

}

/*!
 * \details Destructor that deletes the BaseParticle that is copied and inserted
 *          at every insertion.
 */
InsertionBoundary::~InsertionBoundary()
{
    if (particleToCopy_!=nullptr)
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
    if (particleToCopy!=nullptr)
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
    return P;
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

    /* Each time step, the InsertionBoundary attempts to fill up a region with
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
     * Otherwise, the processes terminates for this time step.
     * */

    // Keep count of how many successive times we have failed to place a new
    // particle. 
    unsigned int failed = 0;
    while (failed <= maxFailed_ && (volumeInserted_<=getVolumeFlowRate()*md->getNextTime())) // 'generating' loop
    {
        /* Generate random *intrinsic* properties for the new particle. */
        logger(VERBOSE, "about to call generateParticle\n");
        auto p0 = generateParticle(md->random);
        logger(VERBOSE, "generated a particle with intrinsics %", p0);

        while (true) // 'placing' loop
        {
            /* Generate extrinsic properties (position and velocity) for this
             * new particle. */
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
                    particle = copyDataFromParticleToMPIParticle(p0);
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
            if (md->checkParticleForInteraction(*p0))
            {
                //Note: in parallel only one of the domains will actually add the particle
                md->particleHandler.copyAndAddObject(p0);
                failed = 0;

                ++numberOfParticlesInserted_;
                massInserted_ += p0->getMass();
                volumeInserted_ += p0->getVolume();
                logger(VERBOSE, "successfully placed a particle %, with position: % after % fails.", p0, p0->getPosition(),failed);

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

double InsertionBoundary::getMassOfParticlesInserted() const
{
    return massInserted_;
}

double InsertionBoundary::getVolumeOfParticlesInserted() const
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

/*!
 * \details Sets the maximum number of times InsertionBoundary::checkBoundaryBeforeTimeStep()
 * may try to insert a particle and fail, before the insertion of particles stops.
 * \param[in] maxFailed     the maximum number of particle insertion trials
 */
void InsertionBoundary::setMaxFailed(unsigned int maxFailed)
{
    maxFailed_=maxFailed;
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
    if (particleToCopy!=nullptr)
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
    is >> dummy >> maxFailed_;
    is >> dummy >> numberOfParticlesInserted_;
    is >> dummy;
    if (particleToCopy_!=nullptr)
        delete particleToCopy_;
    particleToCopy_ = getHandler()->getDPMBase()->particleHandler.readAndCreateObject(is);

    // The .restart file records the index of the particle's species, but
    // doesn't record the pointer, i.e. the memory address of the species within
    // the speciesHandler. The latter needs to be reset now.
    particleToCopy_->setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(
                    particleToCopy_->getIndSpecies() 
                ));

    logger(VERBOSE, "maxFailed_ = %d\n", maxFailed_);
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
    os << " numberOfParticlesInserted " << numberOfParticlesInserted_;
    //os << " particleToCopy " << particleToCopy_;
}

Mdouble InsertionBoundary::getVolumeFlowRate() const {
    return volumeFlowRate_;
}

void InsertionBoundary::setVolumeFlowRate(Mdouble volumeFlowRate) {
    volumeFlowRate_ = volumeFlowRate;
}
