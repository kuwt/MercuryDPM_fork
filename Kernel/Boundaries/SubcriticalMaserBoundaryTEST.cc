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


#include "DPMBase.h"
#include "SubcriticalMaserBoundaryTEST.h"

/*!
 * \details MaserBoundary constructor, sets all scalars to non-sensible values.
 */
SubcriticalMaserBoundaryTEST::SubcriticalMaserBoundaryTEST()
{
    distanceLeft_ = std::numeric_limits<Mdouble>::quiet_NaN();
    distanceRight_ = std::numeric_limits<Mdouble>::quiet_NaN();
    maserIsActivated_ = false;
    activationTime_ = std::numeric_limits<Mdouble>::quiet_NaN();
    copyFlowParticles_ = false;
}

/*!
 * \details destructor
 */
SubcriticalMaserBoundaryTEST::~SubcriticalMaserBoundaryTEST()
{
    logger(DEBUG, "SubcriticalMaserBoundaryTEST::~SubcriticalBoundary() finished");
}

/*!
 * \details Copy method, creates a copy of the object on the heap and returns a
 * pointer to it.
 * \return      pointer to the copy
 */
SubcriticalMaserBoundaryTEST* SubcriticalMaserBoundaryTEST::copy() const
{
    return new SubcriticalMaserBoundaryTEST(*this);
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream from which the boundary must be read
 */
void SubcriticalMaserBoundaryTEST::read(std::istream& is)
{
    PeriodicBoundary::read(is);
    std::string dummy;
    is >> dummy >> maserIsActivated_
       >> dummy >> activationTime_;
    logger(INFO, "maser is activated: %", maserIsActivated_);
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream to which the boundary must be written
 */
void SubcriticalMaserBoundaryTEST::write(std::ostream& os) const
{
    PeriodicBoundary::write(os);
    os << " maserIsActivated " << maserIsActivated_
       << " activationTime " << activationTime_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'SubcriticalMaserBoundaryTEST'
 */
std::string SubcriticalMaserBoundaryTEST::getName() const
{
    return "SubcriticalMaserBoundaryTEST";
}

void SubcriticalMaserBoundaryTEST::setCopyFlowParticles(bool copyFlowParticles)
{
    copyFlowParticles_ = copyFlowParticles;
}

/*!
 * Helper function to make sure that the particles are maser particles after restarting.
 */
void SubcriticalMaserBoundaryTEST::actionsBeforeTimeLoop()
{
    if (maserIsActivated_)
    {
        for (BaseParticle* p : getHandler()->getDPMBase()->particleHandler)
        {
            if (getDistance(p->getPosition()) > 0)
            {
                p->setMaserParticle(true);
            }
        }
    }
}

/*!
 * \details If the maser is active only maser particles can create a ghost of they are close
 * to the right boundary. If this is the case a ghost will be created and added to the particleHandler
 * at the appropriate position. If the maser is not active then this acts as a normal periodic boundary
 * and so the left boundary is also allowed to create ghosts.
 * \param[in] p         Particle to be checked and possibly periodically copied
 * \param[in,out] pH    System's ParticleHandler, (1) from which the interaction radius
 *                      of its largest particle is retrieved to determine the maximum
 *                      distance from the wall at which a particle should still have
 *                      a periodic copy created, and (2) to which a possible periodic
 *                      copy of the particle will be added
 */
void SubcriticalMaserBoundaryTEST::createPeriodicParticle(BaseParticle* p, ParticleHandler& pH)
{
    ///\todo why 2.0?
    const Mdouble proximityDistance = p->getMaxInteractionRadius() + 2.0* pH.getLargestParticle()->getMaxInteractionRadius();
    if (maserIsActivated_)
    {
        if (p->isMaserParticle())
        {
            //Only if the particle is close to the right wall a ghost can be generated
            if (getDistanceFromRight(p->getPosition()) < proximityDistance)
            {
                createGhostParticle(p);
            }
        }
    }
    else
    {
        //If a particle is at any wall, the ghost can be generated
        if (getDistance(p->getPosition()) < proximityDistance)
        {
            createGhostParticle(p);
        }
    }
}

/*!
 * \details Checks whether a given particle (a) is in the Maser and (b) has
 * crossed the closest wall. If so, shifts its position so as to have it appear
 * at the other wall, and creates a 'real' equivalent in the outflow domain.
 * \param[in] p         The particle to be checked and possibly shifted and copied
 * \param pH            The ParticleHandler, which is unused in this implementation
 * \return              Returns true if particle p has interacted with the
 *                      boundary, false otherwise
 */
bool SubcriticalMaserBoundaryTEST::checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH) const
{
    // check if particle passed either of the boundary walls
    if (maserIsActivated_)
    {
        if (p->isMaserParticle() && (getDistanceFromRight(p->getPosition()) < 0))
        {
            //Do the maser stuff
            BaseParticle* pCopy = p->copy();
            pCopy->setMaserParticle(false);
            pH.addObject(pCopy);
            /*
            logger(INFO, "copying particle");
            p->write(std::cout);
            */
            logger(VERBOSE, "copying particle %", p);
            //Shift position of the particle
            shiftPosition(p);
        }
    }
    else
    {
        //note, if we already made a rough bottom on beforehand, we do not want the bottom particles to be eaten up,
        //that is, all shifted inside the domain.
        if (getDistance(p->getPosition()) < 0)
        {
            shiftPosition(p);
        }
    }
    
    return false;
}

/*!
 * \details After particles have moved, check if the maser needs to be update or not.
 * This function also updates the particles in the particleHandler based on their new position
 * \param[in] pH The particle handler, is used to loop over all particles to flag maser particles
 */
void SubcriticalMaserBoundaryTEST::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
#endif
    //Activate the maser at the given time
    if (getHandler()->getDPMBase()->getTime() > activationTime_)
    {
        if (!maserIsActivated_)
        {
            logger(INFO, "Activating the maser at time: %", activationTime_);
            activateMaser();
        }
    }
    
    //Update particles
    for (BaseParticle* p : pH)
    {
        checkBoundaryAfterParticleMoved(p, pH);
    }
#ifdef MERCURY_USE_MPI
    }

#endif
}


/*!
 * \details Activates the maser boundary by flagging all particles within as isMaserParticles.
 * These maser particle will generate real particles at the right side of the maser.
 */
void SubcriticalMaserBoundaryTEST::activateMaser()
{
    logger(INFO, "SubcriticalMaserBoundaryTEST::activateMaser ");
#ifdef MERCURY_USE_MPI
    //Clear the periodic boundaryHandler, only want real particles
    getPeriodicHandler()->clearCommunicationLists();
#endif
    //Flag that the maser is active!
    maserIsActivated_ = true;
    
    extendBottom();
    if (copyFlowParticles_)
    {
        copyExtraParticles();
    }
    
    //Loop over all particles to see if they are within the maser boundary
    ParticleHandler& pH = getHandler()->getDPMBase()->particleHandler;
    for (BaseParticle* particle : pH)
    {
        particle->setMaserParticle((getDistance(particle->getPosition()) > 0));
    }

#ifdef MERCURY_USE_MPI
    //Generate ghost particles
    getPeriodicHandler()->addNewParticles();
#endif

}

void SubcriticalMaserBoundaryTEST::deactivateMaser()
{
    if (maserIsActivated_)
    {
        for (BaseParticle* particle : getHandler()->getDPMBase()->particleHandler)
        {
            if (getDistance(particle->getPosition()) > 0)
            {
                particle->setMaserParticle(false);
            }
        }

        maserIsActivated_ = false;


        // TODO JMFT: @Marnix In activateMaser(), what do extendBottom() and
        // addNewParticles() do, and how do I undo their effects?
    }
    else
        logger(WARN, "[SubcriticalMaserBoundaryTEST::deactivateMaser()] Maser is not activated, so can't deactivate");
}

bool SubcriticalMaserBoundaryTEST::isActivated() const
{
    return maserIsActivated_; 
}

/*!
 * \details The maser is disabled by default and has to be activated. The
 * activation can be done by setting the activation time.
 * \param[in] time The time at which the maser needs to be activated
 */
void SubcriticalMaserBoundaryTEST::setActivationTime(Mdouble time)
{
    activationTime_ = time;
}

/*!
 * \details getDistance is used in a periodic boundary to measure the distance from
 * a certain position to the periodic boundary that actually mirrors are particle.
 * For a maser boundary, if the boundary is active, this is only the right boundary
 * and hence the generic getDistance function is overwritten with getDistanceFromRight()
 * \param[in] position The position of which we want to know the distance towards the nearest
 * active boundary
 */
Mdouble SubcriticalMaserBoundaryTEST::getDistance(const Vec3D& position) const
{
    if (maserIsActivated_)
    {
        return getDistanceFromRight(position);
    }
    else
    {
        return PeriodicBoundary::getDistance(position);
    }
}

/*!
 * \details The maser only requires particles at the right boundary to create ghosts
 * therefore the distance towards the right boundary is an important quantity to compute.
 * \param[in] Position from which the distance to the right wall is computd
 */
Mdouble SubcriticalMaserBoundaryTEST::getDistanceFromRight(const Vec3D& position) const
{
    Mdouble distanceFromPlaneThroughOrigin = Vec3D::dot(position, normal_);
    return distanceRight_ - distanceFromPlaneThroughOrigin;
}

/*!
 * \details In the parallel periodic boundary, everything is computed by the periodic complexity
 * of a particle. Generally it is not possible to be a real particle outside a periodic boundary,
 * but in case of a maser this is definetly possible. To accomodate this the periodic complexity
 * needs to be modified such that this particle remains to be marked as a real particle.
 * When the distance from the maser boundary is negative we give the flag 3 of the periodic
 * complexity to mark that this is actually a real particle outside the maser boundary
 * \param[in,out] complexity The periodic complexity, indicates how a position/particle is related to periodioc boundaries
 * \param[in] position The position of the given periodic complexity
 * \param[in] The index in the complexity vector that this boundary corresponds to
 */
void SubcriticalMaserBoundaryTEST::modifyPeriodicComplexity(std::vector<int>& complexity, int& totalPeriodicComplexity,
                                                            BaseParticle* particle, int i) const
{
    if (maserIsActivated_)
    {
        if (particle->isMaserParticle())
        {
            //Check if this particle is flagged periodic in this boundary
            if (complexity[i] < 0)
            {
                //Make sure that only a ghost close to the right boundary is flagged as real in this boundary
                if (getDistanceFromRight(particle->getPosition()) < 0)
                {
                    complexity[i] = 3;
                }
                
            }
            if (complexity[i] == -3)
            {
                logger(INFO, "Something went wrong in SubcriticalMaserBoundaryTEST::modifyPeriodicComplexity, "
                             "complexity[%] = -3", i);
            }
        }
        else
        {
            if (complexity[i] == 1)
            {
                totalPeriodicComplexity--;
            }
            complexity[i] = 3;
        }
    }
}

void SubcriticalMaserBoundaryTEST::modifyGhostAfterCreation(BaseParticle* particle, int i)
{
    if (maserIsActivated_)
    {
        if (!particle->isMaserParticle())
        {
            std::vector<int> periodicComplexity = particle->getPeriodicComplexity();
            periodicComplexity[i] = 3;
            particle->setPeriodicComplexity(periodicComplexity);
        }
    }
}


/*!
 * \details Before adding particles a check is made to see if the maser needs to be activated or not
 * This check is based on the give activation time, if not set by default it is NaN and the maser
 * will never be activated
 */
void SubcriticalMaserBoundaryTEST::performActionsBeforeAddingParticles()
{
    //Activate the maser at the given time
    if (getHandler()->getDPMBase()->getTime() > activationTime_)
    {
        if (!maserIsActivated_)
        {
            logger(INFO, "Activating the maser at time: %", activationTime_);
            activateMaser();
        }
    }
}

/*!
 * when activating the maser, extend the bottom periodically until the end of the domain, by copying the fixed particles
 * of the periodic part of the maser with intervals of shift_.
 */
void SubcriticalMaserBoundaryTEST::extendBottom() const
{
    logger(INFO, "extending bottom");
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    std::vector<unsigned int> numberOfParticlesPerCore(NUMBER_OF_PROCESSORS);
    std::vector<std::vector<BaseParticle*>> particlesToCores(NUMBER_OF_PROCESSORS);
    if (PROCESSOR_ID == 0)
    {
#endif
    std::vector<BaseParticle*> fixedParticles;
    std::vector<BaseParticle*> newParticles;
    for (BaseParticle* p : getHandler()->getDPMBase()->particleHandler)
    {
        if (p->isFixed() && !(p->isPeriodicGhostParticle()) && !(p->isMPIParticle()))
        {
            fixedParticles.push_back(p);
        }
    }
    
    for (BaseParticle* p : fixedParticles)
    {
        Vec3D newPosition = p->getPosition() + shift_;
        Vec3D maxDomain = getHandler()->getDPMBase()->getMax();
        Vec3D minDomain = getHandler()->getDPMBase()->getMin();
        while (newPosition.X < maxDomain.X && newPosition.Y < maxDomain.Y && newPosition.Z < maxDomain.Z
               && newPosition.X > minDomain.X && newPosition.Y > minDomain.Y && newPosition.Z > minDomain.Z)
        {
            p = p->copy();
            p->setPosition(newPosition);
            newParticles.push_back(p);
            newPosition += shift_;
        }
    }
#ifndef MERCURY_USE_MPI
    for (BaseParticle* p : newParticles)
    {
        getHandler()->getDPMBase()->particleHandler.addObject(p);
    }
#else
    if (NUMBER_OF_PROCESSORS == 1)
    {
        for (BaseParticle* p : newParticles)
        {
            getHandler()->getDPMBase()->particleHandler.addObject(p);
        }
    }
    else
    {
        //count how many particles go to each core
        for (BaseParticle* p : newParticles)
        {
            int targetGlobalIndex = getHandler()->getDPMBase()->domainHandler.getParticleDomainGlobalIndex(p);
            int targetProcessor = getHandler()->getDPMBase()->domainHandler.getParticleProcessor(
                    targetGlobalIndex);
            particlesToCores[targetProcessor].push_back(p);
        }
        for (int i = 0; i < NUMBER_OF_PROCESSORS; ++i)
        {
            numberOfParticlesPerCore[i] = (particlesToCores[i].size());
        }
    }
}
if (NUMBER_OF_PROCESSORS > 1)
{
    //broadcast numberOfParticlesPerCore to other processors
    communicator.broadcast(numberOfParticlesPerCore.data(), NUMBER_OF_PROCESSORS, 0);
    for (unsigned i = 0; i < numberOfParticlesPerCore.size(); ++i)
    {
        for (unsigned p = 0; p < numberOfParticlesPerCore[i]; ++p)
        {
            BaseParticle* particle;
            if (PROCESSOR_ID == 0)
            {
                particle = particlesToCores[i][p];
            }
            else
            {
                particle = new SphericalParticle();
            }

            getHandler()->getDPMBase()->particleHandler.addObject(0, particle);
        }
    }
}
#endif
}

/*!
 * when activating the maser, extend the bottom periodically until the end of the domain, by copying the fixed particles
 * of the periodic part of the maser with intervals of shift_.
 */
void SubcriticalMaserBoundaryTEST::copyExtraParticles() const
{
    logger(INFO, "copying flow particles");
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    std::vector<unsigned int> numberOfParticlesPerCore(NUMBER_OF_PROCESSORS);
    std::vector<std::vector<BaseParticle*>> particlesToCores(NUMBER_OF_PROCESSORS);
    if (PROCESSOR_ID == 0)
    {
#endif
    std::vector<BaseParticle*> flowParticles;
    std::vector<BaseParticle*> newParticles;
    for (BaseParticle* p : getHandler()->getDPMBase()->particleHandler)
    {
        if (!p->isFixed() && !(p->isPeriodicGhostParticle()) && !(p->isMPIParticle()))
        {
            flowParticles.push_back(p);
        }
    }
    
    for (BaseParticle* p : flowParticles)
    {
        Vec3D newPosition = p->getPosition() + shift_;
        for (unsigned i = 0; i < 4; ++i)
        {
            p = p->copy();
            p->setPosition(newPosition);
            newParticles.push_back(p);
            newPosition += shift_;
        }
    }
#ifndef MERCURY_USE_MPI
    for (BaseParticle* p : newParticles)
    {
        getHandler()->getDPMBase()->particleHandler.addObject(p);
    }
#else
    if (NUMBER_OF_PROCESSORS == 1)
    {
        for (BaseParticle* p : newParticles)
        {
            getHandler()->getDPMBase()->particleHandler.addObject(p);
        }
    }
    else
    {
        //count how many particles go to each core
        for (BaseParticle* p : newParticles)
        {
            int targetGlobalIndex = getHandler()->getDPMBase()->domainHandler.getParticleDomainGlobalIndex(p);
            int targetProcessor = getHandler()->getDPMBase()->domainHandler.getParticleProcessor(
                    targetGlobalIndex);
            particlesToCores[targetProcessor].push_back(p);
        }
        for (int i = 0; i < NUMBER_OF_PROCESSORS; ++i)
        {
            numberOfParticlesPerCore[i] = (particlesToCores[i].size());
        }
    }
}
if (NUMBER_OF_PROCESSORS > 1)
{
    //broadcast numberOfParticlesPerCore to other processors
    communicator.broadcast(numberOfParticlesPerCore.data(), NUMBER_OF_PROCESSORS, 0);
    for (unsigned i = 0; i < numberOfParticlesPerCore.size(); ++i)
    {
        for (unsigned p = 0; p < numberOfParticlesPerCore[i]; ++p)
        {
            BaseParticle* particle;
            if (PROCESSOR_ID == 0)
            {
                particle = particlesToCores[i][p];
            }
            else
            {
                particle = new SphericalParticle();
            }

            getHandler()->getDPMBase()->particleHandler.addObject(0, particle);
            //particle->setPeriodicComplexity(getPeriodicHandler()->computePeriodicComplexity(particle->getPosition()));
        }
    }
}
#endif
    
    logger(INFO, "completed copying particles");
}
