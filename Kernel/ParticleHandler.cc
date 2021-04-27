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


#include <limits>
#include <Math/Helpers.h>
#include <Particles/SphericalParticle.h>
#include <Particles/ThermalParticle.h>
#include "ParticleHandler.h"
#include "DPMBase.h"
#include "SpeciesHandler.h"
#include "Species/ParticleSpecies.h"
#include "Particles/LiquidFilmParticle.h"
#include "Particles/BaseParticle.h"
#include "Particles/SuperQuadricParticle.h"
#include "MpiContainer.h"
#include "MpiDataClass.h"
#include "BoundaryHandler.h"
#include "DomainHandler.h"
#include "Domain.h"

/*!
 * \details Constructor of the ParticleHandler class. It creates and empty 
 *          ParticleHandler.
 */
ParticleHandler::ParticleHandler()
{
    largestParticle_ = nullptr;
    smallestParticle_ = nullptr;
    logger(DEBUG, "ParticleHandler::ParticleHandler() finished");
}

/*!
 * \param[in] PH The ParticleHandler that has to be copied. 
 * \details This is not a copy constructor! It copies the DPMBase and all 
 *          BaseParticle, and sets the other variables to 0. After that, it 
 *          computes the smallest and largest particle in this handler.
 */
ParticleHandler::ParticleHandler(const ParticleHandler& PH)
        : BaseHandler<BaseParticle>()
{
    clear();
    setDPMBase(PH.getDPMBase());
    largestParticle_ = nullptr;
    smallestParticle_ = nullptr;
    copyContentsFromOtherHandler(PH);
    if (!objects_.empty())
    {
        computeLargestParticle();
        computeSmallestParticle();
    }
    logger(DEBUG, "ParticleHandler::ParticleHandler(const ParticleHandler &PH) finished");
}

/*!
 * \param[in] rhs The ParticleHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It copies the DPMBase and all 
 *          BaseParticle, and sets the other variables to 0. After that, it 
 *          computes the smallest and largest particle in this handler.
 */
ParticleHandler& ParticleHandler::operator=(const ParticleHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
        largestParticle_ = nullptr;
        smallestParticle_ = nullptr;
        copyContentsFromOtherHandler(rhs);
        if (!objects_.empty())
        {
            computeLargestParticle();
            computeSmallestParticle();
        }
    }
    logger(DEBUG, "ParticleHandler::operator = (const ParticleHandler& rhs) finished");
    return *this;
}

/*!
 * \details Set the pointers to largestParticle_ and smallestParticle_ to 
 *          nullptr, all BaseParticle are destroyed by the BaseHandler afterwards.
 */
ParticleHandler::~ParticleHandler()
{
    //First reset the pointers, such that they are not checked twice when removing particles
    largestParticle_ = nullptr;
    smallestParticle_ = nullptr;
    logger(DEBUG, "ParticleHandler::~ParticleHandler() finished");
}

/*!
 * \param[in] P A pointer to the BaseParticle that has to be added.
 * \details To add a BaseParticle to the ParticleHandler, first check if it has
 *          a species, since it is as common bug to use a BaseParticle without
 *          species, which leads to a segmentation fault. To help the user with
 *          debugging, a warning is given if a particle without species is added.
 *          After that, the actions for adding the particle to the BaseHandler
 *          are taken, which include adding it to the vector of pointers to all
 *          BaseParticle and assigning the correct id and index. Then the
 *          particle is added to the HGrid, the particle is told that this is its
 *          handler, its mass is computed and finally it is checked if this is
 *          the smallest or largest particle in this ParticleHandler.
 *          The particle exists, in other words: it has been made before.
 *          This implies it already got an id Attached to it and hence we don't want
 *          to assign a new ID to it.
 */
void ParticleHandler::addExistingObject(BaseParticle* P)
{
    if (P->getSpecies() == nullptr)
    {
        logger(WARN, "WARNING: The particle with ID % that is added in "
                     "ParticleHandler::addObject does not have a species yet."
                     "Please make sure that you have "
                     "set the species somewhere in the driver code.", P->getId());
    }
    
    //Puts the particle in the Particle list
    BaseHandler<BaseParticle>::addExistingObject(P);
    if (getDPMBase() != nullptr)
    {
        //This places the particle in this grid
        getDPMBase()->hGridInsertParticle(P);
        //This computes where the particle currently is in the grid
        getDPMBase()->hGridUpdateParticle(P);
    }
    //set the particleHandler pointer
    P->setHandler(this);
    //compute mass of the particle
    P->getSpecies()->computeMass(P);
    //Check if this particle has new extrema
    checkExtrema(P);
}


/*!
 * \param[in] P A pointer to the BaseParticle that has to be added. 
 * \details To add a BaseParticle to the ParticleHandler, first check if it has
 *          a species, since it is as common bug to use a BaseParticle without
 *          species, which leads to a segmentation fault. To help the user with
 *          debugging, a warning is given if a particle without species is added.
 *          After that, the actions for adding the particle to the BaseHandler
 *          are taken, which include adding it to the vector of pointers to all 
 *          BaseParticle and assigning the correct id and index. Then the 
 *          particle is added to the HGrid, the particle is told that this is its
 *          handler, its mass is computed and finally it is checked if this is 
 *          the smallest or largest particle in this ParticleHandler.
 */
void ParticleHandler::addObject(BaseParticle* P)
{
    if (P->getSpecies() == nullptr)
    {
        logger(WARN, "WARNING: The particle with ID % that is added in "
                     "ParticleHandler::addObject does not have a species yet."
                     "Please make sure that you have "
                     "set the species somewhere in the driver code.", P->getId());
    }
#ifdef MERCURY_USE_MPI
    bool insertParticle;
    //Check if the particle P should be added to the current domain
    if (NUMBER_OF_PROCESSORS == 1)
    {
        insertParticle = true;
    }
    else
    {
        insertParticle  = getDPMBase()->mpiInsertParticleCheck(P);
    }

    //Add the particle if it is in the mpi domain or if the domain is not yet defined
    if(insertParticle)
    {
#endif
        //Puts the particle in the Particle list
        BaseHandler<BaseParticle>::addObject(P);
        if (getDPMBase() != nullptr)
        {
            //This places the particle in this grid
            getDPMBase()->hGridInsertParticle(P);
            //This computes where the particle currently is in the grid
            getDPMBase()->hGridUpdateParticle(P);
        }
        //set the particleHandler pointer
        P->setHandler(this);
        //compute mass of the particle
        P->getSpecies()->computeMass(P);
        //Check if this particle has new extrema
        checkExtrema(P);
        if (!P->isSphericalParticle())
        {
            getDPMBase()->setRotation(true);
        }
#ifdef MERCURY_USE_MPI
        P->setPeriodicComplexity(std::vector<int>(0));
    }
    else
    {
        logger.assert(!P->isMPIParticle(),"Can't add mpi particle as it does not exist");
        logger.assert(!P->isPeriodicGhostParticle(),"Can't add mpi particle as it does not exist");
        //Somehwere a really new particle has been added, so to keep the ID's globally unique, we also update
        //the unique value on all other processors
        BaseHandler<BaseParticle>::increaseId();
    }

    //Add the particle to the ghost particle lists
    getDPMBase()->insertGhostParticle(P);
    //Check if this new particle requires an update in the mpi grid (interactionDistance).
    getDPMBase()->updateGhostGrid(P);

    //Delete the particle that was supposed to be added (but was not)
    if(!insertParticle)
    {
        delete P;
    }
#endif
}

/*!
 * \brief This function adds a particle to the simulation where the information of the particle is not available
 * by the target processor.
 * \details When a certain processor generates a particle which needs to be inserted in a domain of another processor,
 * this information needs to be transfered before the particle can be actually added
 * \param[in] fromProcessor processor containing the particle data
 * \param[in,out] particle that contains the data and receives the data
 */
void ParticleHandler::addObject(int fromProcessor, BaseParticle* p)
{
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();

    //The processor that contains the particle that needs to be copied needs to identify the target, and communicate this
    MPIParticle pInfo;
    if (communicator.getProcessorID() == fromProcessor)
    {
        pInfo.copyDataFromParticleToMPIParticle(p);
    }

    //Broadcast from processor i
    communicator.broadcast(&pInfo,MercuryMPIType::PARTICLE,fromProcessor);
    copyDataFromMPIParticleToParticle(&pInfo, p, this);

    //All processors now know have the same information and we can proceed with the collective add
    addObject(p);
#else
    logger(WARN, "Function addObject(int fromProcessor, BaseParticle* p) should not be used in serial code");
#endif
}

/*!
 * \brief This function adds a ghost particle from one processor to another processor
 * \details When a periodic ghost particle or a mpi ghost particle is created, the ID should not be updated
 * and hence the addGhostObject needs to be called. Adding a ghost keeps the ID constant.
 * Additionally this function copies the particle information given on processor fromProcessor to all
 * other processors. The target processor then adds the particle.
 * \param[in] fromProcessor processor containing the particle data
 * \param[in] toProcessor processor that needs to add the particle
 * \param[in,out] particle that contains the data and receives the data
 */
void ParticleHandler::addGhostObject(int fromProcessor, int toProcessor, BaseParticle* p)
{
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    if (fromProcessor == toProcessor)
    {
        if (communicator.getProcessorID() == fromProcessor)
        {
            addGhostObject(p);
        }

        return;
    }

    //Create communication information object
    /// \todo MX: maybe implement a direct send here, which avoids the synch step later on.
    MPIParticle pInfo;
    int tag;
    if (communicator.getProcessorID() == fromProcessor)
    {
        pInfo.copyDataFromParticleToMPIParticle(p);
        tag = fromProcessor*MAX_PROC*10 + toProcessor*10 + MercuryMPITag::PARTICLE_DATA;
        communicator.send(&pInfo, MercuryMPIType::PARTICLE, 1, toProcessor, tag);
    }

    if (communicator.getProcessorID() == toProcessor)
    {
        tag = fromProcessor*MAX_PROC*10 + toProcessor*10 + MercuryMPITag::PARTICLE_DATA;
        communicator.receive(&pInfo, MercuryMPIType::PARTICLE, 1, fromProcessor, tag);
    }

    //Sync the communication
    communicator.sync();

    //Add the data to the new particle
    if (communicator.getProcessorID() == toProcessor)
    {
        copyDataFromMPIParticleToParticle(&pInfo, p, this);
    }


    //Only toProcessor adds the particle, quietly without any check with other processors
    //IMPORTANT NOTE: When adding  a periodic particle in parallel, this is performed just before
    //finding new particles. For that reason we dont have to add a ghostparticle to the mpi communication lists
    if (communicator.getProcessorID() == toProcessor)
    {
        addGhostObject(p);
    }
#else
    logger(WARN,
           "Function addGhostObject(int fromProcessor, int toProcessor, BaseParticle* p) should not be used in serial code");
#endif
}

/*!
 * \brief Adds a BaseParticle to the ParticleHandler
 * \details This is a special function that is used in the parallel code. The functions differs from the standard
 * parallel implementation as this function does not check if the particle is an mpi particle and the mpi domain is not updated.
 * When the domain adds mpi particles to the simulation these checks are not required. It also doesnt increase the id of the particle
 * \param[in] A pointer to the BaseParticle that has to be added.
 */
void ParticleHandler::addGhostObject(BaseParticle* P)
{
#ifdef MERCURY_USE_MPI
    if (P->getSpecies() == nullptr)
    {
        logger(WARN, "[ParticleHandler::adGhostObject(BaseParticle*)] "
                "WARNING: The particle with ID % that is added in "
                "ParticleHandler::addObject does not have a species yet."
                "Please make sure that you have "
                "set the species somewhere in the driver code.", P->getId());
    }

    MPIContainer& communicator = MPIContainer::Instance();
    //Puts the particle in the Particle list
    BaseHandler<BaseParticle>::addGhostObject(P);
    if (getDPMBase() != nullptr)
    {
        //This places the particle in this grid
        getDPMBase()->hGridInsertParticle(P);
        //This computes where the particle currently is in the grid
        getDPMBase()->hGridUpdateParticle(P);
    }
    //set the particleHandler pointer
    P->setHandler(this);
    //compute mass of the particle
    P->getSpecies()->computeMass(P) ;
    //Check if this particle has new extrema
    checkExtrema(P);
#else
    logger(INFO,
           "Function ParticleHandler::mpiAddObject(BaseParticle* P) should only be called when compiling with parallel build on");
#endif
}

/*!
 * \param[in] index The index of which BaseParticle has to be removed from this 
 *                  ParticleHandler.
 * \details         The BaseParticle with index is removed and the last 
 *                  BaseParticle in the vector is moved to its position. It also 
 *                  removes the BaseParticle from the HGrid.
 *                  When running this in parallel a warning is thrown that deleting particles might cause havoc
 *                  in the communication between boundaries, as the deleted particle might still be in one of the boundary lists
 */
void ParticleHandler::removeObject(const unsigned int index)
{
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    if (communicator.getNumberOfProcessors() > 1 )
    {
        logger(WARN, "[ParticleHandler::removeObject(const unsigned int)] Using the function removeObject in parallel could lead to out-of-sync communication, Instead use deletion boundaries");
    }
#endif
#ifdef CONTACT_LIST_HGRID
    getDPMBase()->getPossibleContactList().remove_ParticlePosibleContacts(getObject(id));
#endif
    getDPMBase()->hGridRemoveParticle(getObject(index));
    BaseHandler<BaseParticle>::removeObject(index);
}

/*!
 * \param[in] index The index of which BaseParticle has to be removed from this
 *                  ParticleHandler.
 * \details         The BaseParticle with index is removed and the last
 *                  BaseParticle in the vector is moved to its position. It also
 *                  removes the BaseParticle from the HGrid.
 *                  This function is only called by the parallel routines where
 *                  we are sure the particles have been flushed out of the
 *                  communication boundaries.
 */
void ParticleHandler::removeGhostObject(const unsigned int index)
{
#ifdef MERCURY_USE_MPI
#ifdef CONTACT_LIST_HGRID
    getDPMBase()->getPossibleContactList().remove_ParticlePosibleContacts(getObject(id));
#endif
    getDPMBase()->hGridRemoveParticle(getObject(index));
    BaseHandler<BaseParticle>::removeObject(index);
#else
    logger(ERROR, "This function should only be used interally for mpi routines");
#endif
}

/*!
 * \details Function that removes the last object from this ParticleHandler. It
 *          also removes the particle from the HGrid.
 */
void ParticleHandler::removeLastObject()
{
#ifdef CONTACT_LIST_HGRID
    getDPMBase()->getPossibleContactList().remove_ParticlePosibleContacts(getLastObject());
#endif
    getDPMBase()->hGridRemoveParticle(getLastObject());
    BaseHandler<BaseParticle>::removeLastObject();
}

void ParticleHandler::computeSmallestParticle()
{
    if (getSize() == 0)
    {
        logger(DEBUG, "No particles, so cannot compute the smallest particle.");
        
        smallestParticle_ = nullptr;
        return;
    }
    Mdouble min = std::numeric_limits<Mdouble>::max();
    smallestParticle_ = nullptr;
    for (BaseParticle* const particle : objects_)
    {
        if (!(particle->isMPIParticle() || particle->isPeriodicGhostParticle()))
        {
            if (particle->getMaxInteractionRadius() < min)
            {
                min = particle->getMaxInteractionRadius();
                smallestParticle_ = particle;
            }
        }
    }
}

void ParticleHandler::computeLargestParticle()
{
    if (getSize() == 0)
    {
        logger(DEBUG, "No particles, so cannot compute the largest particle.");
        largestParticle_ = nullptr;
        return;
    }
    Mdouble max = -std::numeric_limits<Mdouble>::max();
    largestParticle_ = nullptr;
    for (BaseParticle* const particle : objects_)
    {
        if (!(particle->isMPIParticle() || particle->isPeriodicGhostParticle()))
        {
            if (particle->getMaxInteractionRadius() > max)
            {
                max = particle->getMaxInteractionRadius();
                largestParticle_ = particle;
            }
        }
    }
}

/*!
 * \return  A pointer to the to the smallest BaseParticle (by interactionRadius)
 *          in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getSmallestParticleLocal() const
{
    return smallestParticle_;
}

/*!
 * \return  A pointer to the to the smallest BaseParticle (by interactionRadius)
 *          in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getSmallestParticle() const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"getSmallestParticle should not be used in parallel; use getSmallestInteractionRadius or ParticleSpecies::getSmallestParticleMass instead");
    return nullptr;
#else
    return getSmallestParticleLocal();
#endif
}

/*!
 * \return A pointer to the largest BaseParticle (by interactionRadius) in this 
 *         ParticleHandler.
 */
BaseParticle* ParticleHandler::getLargestParticleLocal() const
{
    return largestParticle_;
}

/*!
 * \return A Pointer to the largest BaseParticle (by interactionRadius) in this ParticleHandler
 */
BaseParticle* ParticleHandler::getLargestParticle() const
{
#ifdef MERCURY_USE_MPIO
    logger(ERROR,"getLargestParticle() should not be used in parallel; use getLargestInteractionRadius instead");
    return nullptr;
#else
    return getLargestParticleLocal();
#endif
}

Mdouble ParticleHandler::getKineticEnergyLocal() const
{
    Mdouble ene = 0;
    for (auto p : *this)
    {
        if (!(p->isMPIParticle() || p->isPeriodicGhostParticle()))
        {
            ene += p->getKineticEnergy();
        }
    }
    return ene;
}

Mdouble ParticleHandler::getKineticEnergy() const
{
#ifdef MERCURY_USE_MPI
    Mdouble kineticEnergyLocal = getKineticEnergyLocal();
    Mdouble kineticEnergyGlobal = 0.0;

    //sum up over all domains
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(kineticEnergyLocal, kineticEnergyGlobal, MPI_SUM);

    return kineticEnergyGlobal;
#else
    return getKineticEnergyLocal();
#endif
}

Mdouble ParticleHandler::getRotationalEnergyLocal() const
{
    Mdouble ene = 0;
    for (auto p : *this)
    {
        if (!(p->isMPIParticle() || p->isPeriodicGhostParticle()))
        {
            ene += p->getRotationalEnergy();
        }
    }
    return ene;
}

Mdouble ParticleHandler::getRotationalEnergy() const
{
#ifdef MERCURY_USE_MPI
    Mdouble rotationalEnergyLocal = getRotationalEnergyLocal();
    Mdouble rotationalEnergyGlobal = 0.0;

    //sum up over all domains
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(rotationalEnergyLocal, rotationalEnergyGlobal, MPI_SUM);

    return rotationalEnergyGlobal;
#else
    return getRotationalEnergyLocal();
#endif
}

Mdouble ParticleHandler::getMassLocal() const
{
    Mdouble m = 0;
    for (auto p : *this)
        if (!(p->isFixed() || p->isMPIParticle() || p->isPeriodicGhostParticle()))
            m += p->getMass();
    return m;
}

Mdouble ParticleHandler::getMass() const
{
#ifdef MERCURY_USE_MPI
    Mdouble massLocal = getMassLocal();
    Mdouble massGlobal = 0.0;

    //Sum up over all domains
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(massLocal, massGlobal, MPI_SUM);

    return massGlobal;
#else
    return getMassLocal();
#endif
}

Vec3D ParticleHandler::getMassTimesPositionLocal() const
{
    Vec3D com = {0, 0, 0};
    for (auto p : *this)
        if (!(p->isFixed() || p->isMPIParticle() || p->isPeriodicGhostParticle()))
            com += p->getMass() * p->getPosition();
    return com;
}

Vec3D ParticleHandler::getMassTimesPosition() const
{
#ifdef MERCURY_USE_MPI
    Vec3D massTimesPositionLocal = getMassTimesPositionLocal();
    Vec3D massTimesPositionGlobal = {0.0, 0.0, 0.0};

    // Sum up over all domains
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(massTimesPositionLocal.X, massTimesPositionGlobal.X, MPI_SUM);
    communicator.allReduce(massTimesPositionLocal.Y, massTimesPositionGlobal.Y, MPI_SUM);
    communicator.allReduce(massTimesPositionLocal.Z, massTimesPositionGlobal.Z, MPI_SUM);

    return massTimesPositionGlobal;
#else
    return getMassTimesPositionLocal();
#endif
}

Vec3D ParticleHandler::getCentreOfMass() const
{
    Mdouble m = getMass();
    if (m == 0)
    {
        Vec3D nanvec = {constants::NaN, constants::NaN, constants::NaN};
        return nanvec;
    }
    else
        return getMassTimesPosition() / m;
}

Vec3D ParticleHandler::getMomentum() const
{
    Vec3D momentum = {0, 0, 0};
    for (auto p : *this)
        if (!(p->isFixed() || p->isMPIParticle() || p->isPeriodicGhostParticle()))
            momentum += p->getMomentum();
    return getMPISum(momentum);
}

Vec3D ParticleHandler::getAngularMomentum() const
{
    Vec3D momentum = {0, 0, 0};
    for (auto p : *this)
        if (!(p->isFixed() || p->isMPIParticle() || p->isPeriodicGhostParticle()))
            momentum += p->getAngularMomentum();
    return getMPISum(momentum);
}

/*!
 * \return A pointer to the fastest BaseParticle in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getFastestParticleLocal() const
{
    if (getSize() == 0)
    {
        logger(WARN, "No particles to set getFastestParticle()");
        return nullptr;
    }
    BaseParticle* p = nullptr;
    Mdouble maxSpeed = -std::numeric_limits<Mdouble>::max();
    for (BaseParticle* const pLoop : objects_)
    {
        if (!(pLoop->isMPIParticle() || pLoop->isPeriodicGhostParticle()))
        {
            if ((pLoop->getVelocity().getLength()) > maxSpeed)
            {
                maxSpeed = pLoop->getVelocity().getLength();
                p = pLoop;
            }
        }
    }
    return p;
}

BaseParticle* ParticleHandler::getFastestParticle() const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"This function should not be used in parallel");
#endif
    return getFastestParticleLocal();
}

/*
 * \returns The smallest interaction radius on a local domain
 */
Mdouble ParticleHandler::getSmallestInteractionRadiusLocal() const
{
    if (!(getSmallestParticleLocal() == nullptr))
    {
        return getSmallestParticleLocal()->getMaxInteractionRadius();
    }
    else
    {
        return 0.0;
    }
}

/*
 * \returns The smallest interaction radius
 */
Mdouble ParticleHandler::getSmallestInteractionRadius() const
{
#ifdef MERCURY_USE_MPI
    //Compute the local value
    Mdouble smallestInteractionRadiusLocal = getSmallestInteractionRadiusLocal();
    Mdouble smallestInteractionRadiusGlobal = 0.0;

    //Obtain the global value
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(smallestInteractionRadiusLocal, smallestInteractionRadiusGlobal, MPI_MIN);

    return smallestInteractionRadiusGlobal;

#else
    return getSmallestInteractionRadiusLocal();
#endif
}

/*
 * \returns the largest interaction radius on a local domain
 */
Mdouble ParticleHandler::getLargestInteractionRadiusLocal() const
{
    if (!(getLargestParticleLocal() == nullptr))
    {
        return getLargestParticle()->getMaxInteractionRadius();
    }
    else
    {
        return 0.0;
    }
}

/*
 * \returns the largest interaction radius
 */
Mdouble ParticleHandler::getLargestInteractionRadius() const
{
#ifdef MERCURY_USE_MPI
    Mdouble largestInteractionRadiusLocal = getLargestInteractionRadiusLocal();
    Mdouble largestInteractionRadiusGlobal = 0.0;

    //Obtain the global value
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(largestInteractionRadiusLocal, largestInteractionRadiusGlobal, MPI_MAX);

    return largestInteractionRadiusGlobal;
#else
    return getLargestInteractionRadiusLocal();
#endif
}

/*!
 * \return      The sum of all particle radii
 */
Mdouble ParticleHandler::getSumRadiusLocal() const
{
    Mdouble sumRadius = 0;
    for (BaseParticle* const p : objects_)
    {
        if (!(p->isMPIParticle() || p->isPeriodicGhostParticle()))
        {
            sumRadius += p->getRadius();
        }
    }
    return sumRadius;
}

Mdouble ParticleHandler::getMeanRadius() const
{
#ifdef MERCURY_USE_MPI
    Mdouble sumRadiusLocal = getSumRadiusLocal();
    unsigned numberOfRealParticlesLocal = getNumberOfRealObjectsLocal();

    Mdouble sumRadiusGlobal = 0.0;
    unsigned numberOfRealParticlesGlobal = 0;

    //Sum up over all domains
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(sumRadiusLocal, sumRadiusGlobal, MPI_SUM);
    communicator.allReduce(numberOfRealParticlesLocal, numberOfRealParticlesGlobal, MPI_SUM);

    return sumRadiusGlobal/numberOfRealParticlesGlobal;
#else
    return getSumRadiusLocal() / getSize();
#endif
}

/*!
 * \param[in] i Direction for which you want the particle with lowest coordinates.
 * \return      A pointer to the particle with the lowest coordinates in the
 *              given direction in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getLowestPositionComponentParticleLocal(const int i) const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"getLowestPositionComponentParticle() not implemented yet in parallel");
#endif
    if (getSize() == 0)
    {
        logger(WARN, "No getLowestPositionComponentParticle(const int i) since there are no particles.");
        return nullptr;
    }
    BaseParticle* p = nullptr;
    Mdouble min = std::numeric_limits<Mdouble>::max();
    for (BaseParticle* const pLoop : objects_)
    {
        if (!(pLoop->isMPIParticle() || pLoop->isPeriodicGhostParticle()))
        {
            if (pLoop->getPosition().getComponent(i) < min)
            {
                min = pLoop->getPosition().getComponent(i);
                p = pLoop;
            }
        }
    }
    return p;
}

BaseParticle* ParticleHandler::getLowestPositionComponentParticle(const int i) const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"This function should not be used in parallel");
#endif
    return getLowestPositionComponentParticleLocal(i);
}

/*!
 * \param[in] i Direction for which one wants the particle with highest coordinates.
 * \return      A pointer to the particle with the highest coordinates in 
 *              direction i in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getHighestPositionComponentParticleLocal(const int i) const
{
    if (getSize() == 0)
    {
        logger(WARN, "No getHighestPositionComponentParticle(const int i) since there are no particles.");
        return nullptr;
    }
    BaseParticle* p = nullptr;
    Mdouble max = -std::numeric_limits<Mdouble>::max();
    for (BaseParticle* const pLoop : objects_)
    {
        if (!(pLoop->isMPIParticle() || pLoop->isPeriodicGhostParticle()))
        {
            if (pLoop->getPosition().getComponent(i) > max)
            {
                max = pLoop->getPosition().getComponent(i);
                p = pLoop;
            }
        }
    }
    
    return p;
}

BaseParticle* ParticleHandler::getHighestPositionComponentParticle(const int i) const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"This function should not be used in parallel");
#endif
    return getHighestPositionComponentParticleLocal(i);
}

/*!
 * \param[in] i Direction for which you want the particle with lowest velocity.
 * \return      A pointer to the particle with the lowest velocity in direction 
 *              i in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getLowestVelocityComponentParticleLocal(const int i) const
{
    if (getSize() == 0)
    {
        logger(WARN, "No getLowestVelocityComponentParticle(const int i) since there are no particles");
        return nullptr;
    }
    BaseParticle* p = nullptr;
    Mdouble min = std::numeric_limits<Mdouble>::max();
    for (BaseParticle* const pLoop : objects_)
    {
        if (!(pLoop->isMPIParticle() || pLoop->isPeriodicGhostParticle()))
        {
            if (pLoop->getVelocity().getComponent(i) < min)
            {
                min = pLoop->getVelocity().getComponent(i);
                p = pLoop;
            }
        }
    }
    return p;
}

BaseParticle* ParticleHandler::getLowestVelocityComponentParticle(const int i) const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"This function should not be used in parallel");
#endif
    return getLowestVelocityComponentParticleLocal(i);
}

/*!
 * \param[in] i Direction for which you want the particle with highest velocity.
 * \return      A pointer to the particle with the highest velocity in direction
 *              i in this ParticleHandler.
 */
BaseParticle* ParticleHandler::getHighestVelocityComponentParticleLocal(const int i) const
{
    if (!getSize())
    {
        logger(WARN, "No getHighestVelocityComponentParticle(const int i) since there are no particles");
        return nullptr;
    }
    BaseParticle* p = nullptr;
    Mdouble max = -std::numeric_limits<Mdouble>::max();
    for (BaseParticle* const pLoop : objects_)
    {
        if (!(pLoop->isMPIParticle() || pLoop->isPeriodicGhostParticle()))
        {
            if (pLoop->getVelocity().getComponent(i) > max)
            {
                max = pLoop->getVelocity().getComponent(i);
                p = pLoop;
            }
        }
    }
    return p;
}

BaseParticle* ParticleHandler::getHighestVelocityComponentParticle(const int i) const
{
#ifdef MERCURY_USE_MPI
    logger(ERROR,"This function should not be used in parallel");
#endif
    return getHighestVelocityComponentParticleLocal(i);
}

/*!
 * \details Note that the pointers to smallestParticle_ and largestParticle_ are
 *          set to nullptr since these particles don't exist anymore after 
 *          calling this function.
 */
void ParticleHandler::clear()
{
    smallestParticle_ = nullptr;
    largestParticle_ = nullptr;
    BaseHandler<BaseParticle>::clear();
}

/*!
 * \brief Function returns the highest position in the x-direction
 * \details This is a prototype example of how to obtain global particle attributes in the parallel code
 * \return Returns the highest x-position value of all particles
 */
Mdouble ParticleHandler::getHighestPositionX() const
{
    //Define the attribute function
    std::function<Mdouble(BaseParticle*)> particleAttribute = [](BaseParticle* p) { return p->getPosition().X; };
    
    //Obtain the MAX attribute
    Mdouble positionXGlobal = getParticleAttribute<Mdouble>(particleAttribute, AttributeType::MAX);
    
    return positionXGlobal;
}

/*!
 *
 */
unsigned int ParticleHandler::getNumberOfFixedParticles() const
{
    return NFixedParticles_;
}

/*!
 *
 */
unsigned int ParticleHandler::getNumberOfUnfixedParticles() const
{
    return getNumberOfObjects() - NFixedParticles_;
}

/*!
 * \param[in] is The input stream from which the information is read.
 */
BaseParticle* ParticleHandler::createObject(const std::string& type)
{
    if (type == "BaseParticle") {
        //for backwards compatibility
        return new SphericalParticle;
    }
    else if (type == "SphericalParticle")
    {
        return new SphericalParticle;
    }
    else if (type == "LiquidFilmParticle")
    {
        return new LiquidFilmParticle;
    }
    else if (type == "SuperQuadricParticle")
    {
        return new SuperQuadricParticle;
    }
    else if (type == "ThermalParticle")
    {
        return new ThermalParticle;
    }
    else
    {
        logger(ERROR, "Particle type % not understood in restart file. Particle will not be read.", type);
        return nullptr;
    }
}

/*!
 * \param[in] is The input stream from which the information is read.
 */
BaseParticle* ParticleHandler::readAndCreateObject(std::istream& is)
{
    std::string type;
    BaseParticle* particle = nullptr;
    if (getStorageCapacity() > 0)
    {
        is >> type;
        logger(DEBUG, "ParticleHandler::readAndCreateObject(is): type %", type);
        particle = createObject(type);
        particle->setHandler(this);
        particle->read(is);
    }
    else //for insertion boundaries
    {
        is >> type;
        logger(DEBUG, "ParticleHandler::readAndCreateObject(is): type %", type);
        particle = createObject(type);
        particle->setHandler(this);
        particle->read(is);
    }
    particle->setSpecies(getDPMBase()->speciesHandler.getObject(particle->getIndSpecies()));
    return particle;
}

/*!
 * \param[in] is The input stream from which the information is read.
 */
void ParticleHandler::readAndAddObject(std::istream& is)
{
    BaseParticle* o = readAndCreateObject(is);
    addExistingObject(o);
}

///*!
// * \param[in] type The first value of the position.
// * \param[in] is The input stream from which the information is read.
// * \details The old objects did not have their type in the beginning of the line.
// *          Instead, the first string of the file was the position in x-direction.
// *          Since we already read the first string of the file, we need to give
// *          it to this function and convert it to the position in x-direction.
// *          The rest of the stream is then read in the usual way.
// */
//void ParticleHandler::readAndCreateOldObject(std::istream& is, const std::string& type)
//{
//    //read in next line
//    std::stringstream line;
//    helpers::getLineFromStringStream(is, line);
//    logger(VERBOSE, line.str());
//    //std::cout << line.str() << std::endl;
//
//    BaseParticle particle;
//
//    //Declare all properties of the particle
//    unsigned int indSpecies;
//    Mdouble radius, inverseMass, inverseInertia;
//    Vec3D position, velocity, euler, angularVelocity;
//
//    //Read all values
//    position.X = atof(type.c_str());
//
//    line >> position.Y >> position.Z >> velocity >> radius >> euler >> angularVelocity >> inverseMass >> inverseInertia >> indSpecies;
//
//    //Put the values in the particle
//    particle.setSpecies(getDPMBase()->speciesHandler.getObject(indSpecies));
//    particle.setPosition(position);
//    particle.setVelocity(velocity);
//    particle.setRadius(radius);
//    Quaternion q; q.setEuler(euler);
//    particle.setOrientation(q);
//    particle.setAngularVelocity(angularVelocity);
//    if (inverseMass == 0.0)
//        particle.fixParticle();
//    else
//    {
//        particle.setInverseInertia(MatrixSymmetric3D(1,0,0,1,0,1)*inverseInertia);
//    }
//
//    //Put the particle in the  handler
//    copyAndAddObject(particle);
//}

void ParticleHandler::write(std::ostream& os) const
{
#ifdef MERCURY_USE_MPI
    os << "Particles " << getNumberOfRealObjectsLocal() << std::endl;
    for (BaseParticle* it : *this)
    {
        if (!it->isPeriodicGhostParticle() && !it->isMPIParticle())
        {
            os << (*it) << '\n';
        }
    }
#else
    os << "Particles " << getSize() << '\n';
    for (BaseParticle* it : *this)
    {
        os << (*it) << '\n';
    }
#endif
    //os.flush();
}

/*!
 *  \param[in] P A pointer to the particle, which properties have to be checked 
 *               against the ParticleHandlers extrema.
 */
void ParticleHandler::checkExtrema(BaseParticle* P)
{
    if (P == largestParticle_)
    {
        //if the properties of the largest particle changes
        computeLargestParticle();
    }
    else if (!largestParticle_ || P->getMaxInteractionRadius() > largestParticle_->getMaxInteractionRadius())
    {
        largestParticle_ = P;
    }
    
    if (P == smallestParticle_)
    {
        //if the properties of the smallest particle changes
        computeSmallestParticle();
    }
    else if (!smallestParticle_ || P->getMaxInteractionRadius() < smallestParticle_->getMaxInteractionRadius())
    {
        smallestParticle_ = P;
    }
}

/*!
 * \param[in] P A pointer to the particle, which is going to get deleted.
 */
void ParticleHandler::checkExtremaOnDelete(BaseParticle* P)
{
    if (P == largestParticle_)
    {
        computeLargestParticle();
    }
    if (P == smallestParticle_)
    {
        computeSmallestParticle();
    }
}

/*!
 * \param[in] indSpecies Unsigned integer with the index of the species for which
 *                       the masses must be computed.
 */
void ParticleHandler::computeAllMasses(unsigned int indSpecies)
{
    for (BaseParticle* particle : objects_)
    {
        if (particle->getIndSpecies() == indSpecies)
        {
            particle->getSpecies()->computeMass(particle);
        }
    }
}

void ParticleHandler::computeAllMasses()
{
    for (BaseParticle* particle : objects_)
    {
        particle->getSpecies()->computeMass(particle);
    }
}

void ParticleHandler::addedFixedParticle()
{
    NFixedParticles_++;
}

void ParticleHandler::removedFixedParticle()
{
    NFixedParticles_--;
}

/*!
 * \return The string "ParticleHandler".
 */
std::string ParticleHandler::getName() const
{
    return "ParticleHandler";
}

Mdouble ParticleHandler::getVolumeLocal() const
{
    Mdouble volume = 0;
    for (auto p : *this)
    {
        if (!(p->isFixed() || p->isPeriodicGhostParticle() || p->isMPIParticle()))
            volume += p->getVolume();
    }
    return volume;
}

Mdouble ParticleHandler::getVolume() const
{
#ifdef MERCURY_USE_MPI
    Mdouble volumeLocal = getVolumeLocal();
    Mdouble volumeGlobal = 0.0;

    // sum up over all domains
    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(volumeLocal, volumeGlobal, MPI_SUM);

    return volumeGlobal;
#else
    return getVolumeLocal();
#endif
}

/*!
 * \returns the number of real particles (neglecting MPI and periodic) on a local domain
 * \todo MX: in future also add the periodic mpi particles
 */
unsigned int ParticleHandler::getNumberOfRealObjectsLocal() const
{
#ifdef MERCURY_USE_MPI
    const MPIContainer& communicator = MPIContainer::Instance();
    if (communicator.getNumberOfProcessors() > 1)
    {
        unsigned int numberOfFakeParticles = 0;
        numberOfFakeParticles += getDPMBase()->domainHandler.getCurrentDomain()->getNumberOfTrueMPIParticles();
        numberOfFakeParticles += getDPMBase()->periodicBoundaryHandler.getNumberOfPeriodicGhostParticles();
        logger.assert(numberOfFakeParticles <= getSize(), "More fake particles than getSize()");
        return (getSize() - numberOfFakeParticles);
    }
    else
    {
        return getSize();
    }
#else
    return getSize();
#endif
}

unsigned int ParticleHandler::getNumberOfRealObjects() const
{
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    unsigned int numberOfRealParticles = getNumberOfRealObjectsLocal();

    // \todo MX: use an allreduce here
    //Combine the total number of Particles into one number on processor 0
    communicator.reduce(numberOfRealParticles, MPI::SUM);

    //Broadcast new number to all the processorsi
    communicator.broadcast(numberOfRealParticles);
    return numberOfRealParticles;
#else
    return getSize();
#endif
}

/*
 * \returns the size of the ParticleHandler. In parallel code it throws an error to avoid confusion with real and fake particles
 */
unsigned int ParticleHandler::getNumberOfObjects() const
{
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    if (communicator.getNumberOfProcessors() > 1)
    {
        logger(WARN,"When compiling with MPI please do not use getNumberOfObjects(). Instead use: getNumberOfRealObjectsLocal(), getNumberOfRealObjects() or getSize()");
    }
#endif
    return getSize();
}

/*!
 * \brief Computes the number of Fixed particles on a local domain
 * \details Loops over all particles to check if the particle is fixed or not, in a local domain.
 * \returns the number of fixed particles in a local domain
 */
unsigned int ParticleHandler::getNumberOfFixedObjectsLocal() const
{
    unsigned int numberOfFixedParticles = 0;
    for (BaseParticle* particle : *this)
    {
        if (particle->isFixed())
        {
            numberOfFixedParticles++;
        }
    }
    return numberOfFixedParticles;
}

/*!
 * \brief Computes the number of fixed particles in the whole simulation
 */
unsigned int ParticleHandler::getNumberOfFixedObjects() const
{
#ifdef MERCURY_USE_MPI
    unsigned int numberOfFixedParticlesLocal = getNumberOfFixedObjectsLocal();
    unsigned int numberOfFixedParticles = 0;

    MPIContainer& communicator = MPIContainer::Instance();
    communicator.allReduce(numberOfFixedParticlesLocal, numberOfFixedParticles, MPI_SUM);
    return numberOfFixedParticles;
#else
    return getNumberOfFixedObjectsLocal();
#endif
}

void ParticleHandler::actionsAfterTimeStep()
{
    for (auto i: *this)
    {
        i->actionsAfterTimeStep();
    }
}

double ParticleHandler::getLiquidFilmVolume() const {
    double liquidVolume = 0;
    for (auto i : objects_) {
        auto j = dynamic_cast<LiquidFilmParticle*>(i);
        if (j and !j->isMPIParticle()) liquidVolume += j->getLiquidVolume();
    }
    return getMPISum(liquidVolume);
};
