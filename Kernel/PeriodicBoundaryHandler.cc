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


#include <iostream>
#include <iomanip>
#include <vector>
#include "Math/Helpers.h"

#include "PeriodicBoundaryHandler.h"
#include "Boundaries/BasePeriodicBoundary.h"
#include "DPMBase.h"
#include "ParticleHandler.h"
#include <set>
#include <Particles/LiquidFilmParticle.h>
#include "Logger.h"
#include "Particles/SphericalParticle.h"

///Constructor of the PeriodicBoundaryHandler class. It creates and empty PeriodicBoundaryHandler.
PeriodicBoundaryHandler::PeriodicBoundaryHandler()
{
    interactionDistance_ = 0;
    logger(DEBUG, "PeriodicBoundaryHandler::PeriodicBoundaryHandler() finished");
}

/*!
 * \param[in] BH The PeriodicBoundaryHandler that has to be copied.
 * \details This is not a copy constructor! It just copies all BasePeriodicBoundary from
 *          the other handler into this handler, and clears all other variables.
 */
PeriodicBoundaryHandler::PeriodicBoundaryHandler(const PeriodicBoundaryHandler& PBH)
{
    objects_.clear();
    setDPMBase(nullptr);
    interactionDistance_ = PBH.interactionDistance_;
    logger(DEBUG, "PeriodicBoundaryHandler::PeriodicBoundaryHandler(const PeriodicBoundaryHandler &BH) finished");
}

/*!
 * \param[in] rhs The PeriodicBoundaryHandler on the right hand side of the assignment.
 * \details This is not a copy assignment operator! It just copies all BaseBoundary
 *          from the other handler into this handler, and clears all other variables.
 */
PeriodicBoundaryHandler PeriodicBoundaryHandler::operator=(const PeriodicBoundaryHandler& rhs)
{
    if (this != &rhs)
    {
        objects_.clear();
        interactionDistance_ = rhs.interactionDistance_;
    }
    logger(DEBUG, "PeriodicBoundaryHandler PeriodicBoundaryHandler::operator =(const BoundaryHandler& rhs)");
    return *this;
}

///Default destructor. Note that the delete for all boundaries is done in the BaseHandler.
PeriodicBoundaryHandler::~PeriodicBoundaryHandler()
{
    /// \todo MX: Because the periodicBoundaryHandler contains a pointer to a boundary, which is also in the boundary handler
    /// currently we just need to empty the PeriodicBoundaryHandler without destroying the object
    /// Otherwise we create a segfault in the BoundaryHandler
    objects_.clear();
    
    logger(DEBUG, "PeriodicBoundaryHandler::~PeriodicBoundaryHandler() finished");
}

///\param[in] P A pointer to the BasePeriodicBoundary (or derived class) that has to be added.
///Add the object and tell the object that this is his handler.
void PeriodicBoundaryHandler::addObject(BasePeriodicBoundary* P)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
        return;
    }
    //Puts the particle in the Particle list
    objects_.push_back(P);
    //set the particleHandler pointer
    P->setPeriodicHandler(this);
    //Remove all ghost particles
    clearCommunicationLists();
    //Initialise ghost particles
    addNewParticles();
#endif
}

///The periodic boundary handler does not need to read anything, it will be reconstructed
void PeriodicBoundaryHandler::readAndAddObject(std::istream& is)
{
}

/*!
 * \details This determines the unique name of the current class
 * \return The string "PeriodicBoundaryHandler"
 */
std::string PeriodicBoundaryHandler::getName() const
{
    return "PeriodicBoundaryHandler";
}

/*!
 * \brief Sets the interaction distance.
 * \details This distance determines when a particle starts to interact with a periodic boundary
 * \param[in] interactionDistance The distance from a boundary when a particle starts to interact
 * with it.
 */
void PeriodicBoundaryHandler::setInteractionDistance(Mdouble interactionDistance)
{
    interactionDistance_ = interactionDistance;
}

/*!
 * \details This distance determines when a particle starts to interact with a periodic boundary.
 * \return The interaction distance. 
 */
Mdouble PeriodicBoundaryHandler::getInteractionDistance()
{
    return interactionDistance_;
}

/*!
 * \detail This function updates the status of periodic particles and periodic ghost particles.
 * This is done in two steps. The first step is to update the ghost particles with the position
 * of their corresponding real particles. Based on this new position their status will be updated.
 * Occasionally a particle needs to be removed and because this particle might also be listed
 * in the communication boundaries, it will be deleted after those lists have been flushed.
 * \param[in] particlesToBeDeleted A list of particles that need to be deleted, but can't yet be deleted
 * to avoid segmentation faults
 */
void PeriodicBoundaryHandler::updateStatus(std::set<BaseParticle*>& particlesToBeDeleted)
{
    //Step 1: collect position and velocity data of the real particle and send the data
    preparePositionAndVelocityUpdate();
    
    //Step 2: finalise the data transision
    finalisePositionAndVelocityUpdate();
    
    //Step 3: Update the status of the ghost and real particles
    updateParticleStatus(particlesToBeDeleted);
}

/*!
 * \detail This function shifts the position of the particle with respect to the periodic boundaries
 * based on the periodic complexity it currently has. Note that some boundaries such as the angular
 * periodic boundary do not only shift the position, but also the velocity. In that sense this
 * function name is a bit of a misnomer.
 * \param[in] particle The particle that shifts position (and possibly other properties)
 */
void PeriodicBoundaryHandler::shiftParticle(BaseParticle* particle)
{
    shiftParticle(particle, particle->getPeriodicComplexity());
}

/*!
 * \detail This function shifts the position of the particle with respect to the periodic boundaries
 * based on a given periodic complexity. Note that some boundaries such as the angular
 * periodic boundary do not only shift the position, but also the velocity. In that sense this
 * function name is a bit of a misnomer.
 * \param[in] particle The particle that shifts position (and possibly other properties)
 * \param[in] complexity The periodic complexity that determines how the particle is shifted
 */
void PeriodicBoundaryHandler::shiftParticle(BaseParticle* particle, const std::vector<int>& complexity)
{
    int boundaryIndex = 0;
    for (BasePeriodicBoundary* boundary : *this)
    {
        if (complexity[boundaryIndex] == -1)
        {
            boundary->shiftPosition(particle);
        }
        boundaryIndex++;
    }
}

/*!
 * \detail This function computes the periodic complexity, a vector of intergers that indicate how this
 * particle is related to the periodic boundaries. Every boundary has a status value following:
 * [2]  = Real particle, not in the proximity of the boundary
 * [1]  = Real particle, in the proximity of the boundary
 * [-1] = ghost particle, in the proximity of the boundary
 * [-2] = ghost particle, not in the proximity of the boundary
 * In addition to the periodic complexity it also computes the totalPeriodicComplexity that indicates with
 * how many boundaries this position is in close proximity
 * \param[in,out] periodicComplexity The periodic complexity indicating how a given position is related 
 *  to the periodic boundaries.
 * \param[in,out] totalPeriodicComplexity The number of boundaries the position is in close proximity with.
 * \param[in] position The input position for which the relation to the boundaries are determined.
 */
void
PeriodicBoundaryHandler::computePeriodicComplexity(std::vector<int>& periodicComplexity, int& totalPeriodicComplexity,
                                                   Vec3D position)
{
    //Compute the new position
    periodicComplexity.resize(this->getSize());
    totalPeriodicComplexity = 0;
    int index = 0;
    for (BasePeriodicBoundary* boundary : *this)
    {
        Mdouble distance = boundary->getDistance(position);
        if (std::abs(distance) <= interactionDistance_)
        {
            if (distance > 0)
            {
                //real particle
                periodicComplexity[index] = 1;
                totalPeriodicComplexity++;
            }
            else if (distance < 0)
            {
                //ghost particle
                periodicComplexity[index] = -1;
            }
            else //distance == 0 In extremely rare cases (test cases mostly)
            {
                //Check if on left or right side of boundary
                if (boundary->isClosestToLeftBoundary(position))
                {
                    //ghost particle on the left side
                    periodicComplexity[index] = -1;
                }
                else
                {
                    //real particle on the right side
                    periodicComplexity[index] = 1;
                    totalPeriodicComplexity++;
                }
            }
        }
        else
        {
            if (distance < 0)
            {
                periodicComplexity[index] = -2;
            }
            else
            {
                periodicComplexity[index] = 2;
            }
        }
        index++;
    }
/*
    //Modify periodic complexity in case of maser
    for (int b = 0; b < getSize(); b++)
    {
        objects_[b]->modifyPeriodicComplexity(periodicComplexity, position, b);
    }
*/
}

/*!
 * \detail This function computes the periodic complexity, a vector of intergers that indicate how this
 * particle is related to the periodic boundaries. Every boundary has a status value following:
 * [2]  = Real particle, not in the proximity of the boundary
 * [1]  = Real particle, in the proximity of the boundary
 * [-1] = ghost particle, in the proximity of the boundary
 * [-2] = ghost particle, not in the proximity of the boundary
 * \param[in] position The input position for which the relation to the boundaries are determined.
 * \return Periodic complexity of the given position
 */
std::vector<int> PeriodicBoundaryHandler::computePeriodicComplexity(Vec3D position)
{
    std::vector<int> periodicComplexity;
    int totalPeriodicComplexity;
    
    computePeriodicComplexity(periodicComplexity, totalPeriodicComplexity, position);
    
    return periodicComplexity;
}

/*!
 * \detail This function adds new particles to the periodiocBoundary lists. Particles
 * that are not yet flagged as periodic(ghost) particles are checked if they entered the interaction 
 * distances of periodic boundaries. If this is the case these will be communicated to the other
 * processors in a three step process. First the processors communicate which each other to
 * determine who is getting how many particles from who. Secondly the actual data transfer is
 * performed and thirdly the received data is processed into the correct bookkeeping lists.
 */
void PeriodicBoundaryHandler::addNewParticles()
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
        return;
    }
    //Step 0: Perform actions before finding new particles. I.e. opening maser boundary
    performActionsBeforeAddingParticles();

    //Step 1: Find new particles, such that the target domains can be found
    findNewParticles();

    //Step 3: Communicate number of particles/interactions to other dommains
    prepareNewParticleTransmission();

    //Step 4: Perform data transmission
    performNewParticleTransmission();

    //Step 5: Finalise data transmission
    finaliseNewParticleTransmission();
#endif
}

/*!
 * \detail When a user or an insertion boundary is adding a single new particle to the system,
 * this function will add it to the corresponding lists if required. If this is the case, the 
 * particle will be communicated to the receiving processor in a three step process. 
 * First the processors communicate which each other to determine who is getting how many particles 
 * from who. Secondly the actual data transfer is performed and thirdly the received data is processed 
 * into the correct bookkeeping lists.
 * \param[in] particle The particle that is (or not) being add to the periodic lists.
 */
void PeriodicBoundaryHandler::addNewParticle(BaseParticle* particle)
{
#ifdef MERCURY_USE_MPI
    if (NUMBER_OF_PROCESSORS == 1)
    {
        return;
    }
    if (getDPMBase()->getCurrentDomain()->containsParticle(particle))
    {
        //Step 1: Find if the new particle needs to be added
        findNewParticle(particle);
    }
    //Follow the standard communication routine
    prepareNewParticleTransmission();
    performNewParticleTransmission();
    finaliseNewParticleTransmission();
#endif
}

/*!
 * \details Counts the number of periodic ghost particles by simply checking how many entries there
 *  are in the periodicGhostLists.
 * \return The number of particles that have isPeriodicGhostParticle() == true
 */
unsigned int PeriodicBoundaryHandler::getNumberOfPeriodicGhostParticles()
{
    unsigned int sum = 0;
    for (auto& index : periodicGhostList_)
    {
        sum += index.size();
    }
    return sum;
}

/*!
 * \details Counts the number of true periodic ghost particles by checking how many entries
 * there are in the periodiocGhostLists, but ignores mixed ghost particles that additionally have
 * the flay isMPIParticle() == true.
 * \return The number of periodic particles that are not also MPIParticles
 */
Mdouble PeriodicBoundaryHandler::getNumberOfTruePeriodicGhostParticles()
{
    int sum = 0;
    for (auto& index : periodicGhostList_)
    {
        int numberOfMPIParticles = 0;
        for (int pIndex = 0; pIndex < index.size(); pIndex++)
        {
            if (index[pIndex]->particle->isMPIParticle() == true)
            {
                numberOfMPIParticles++;
            }
        }
        sum += (index.size() - numberOfMPIParticles);
    }
    return sum;
}

/*!
 * \details Mixed particles that are both a interacting with a periodic boundary and also in an mpi boundary are
 * special cases that need to be handled with care. As example a real particle that moves to another domain
 * will be removed by the mpi boundaries. To avoid segmentation faults the particle needs to be flushed from
 * the periodic boundaries as well. This function determines if a particle is in the MPI domain and
 * if it is a MPIParticle or not. 
 * \param[in] particle The particle for which the MPI flags are being computed
 * \param[in,out] isInMPIDomain A check to see if the particle is in the MPI domain (true) or not (false)
 * \param[in,out] isMPIParticle A bool to see if the particle is a ghost in the MPI domain (true) or not (false)
*/
void PeriodicBoundaryHandler::getMPIFlags(BaseParticle* particle, bool& isInMPIDomain, bool& isMPIParticle)
{
    //Check if the particle is actually in this domain. If not the case it is an "MG"-particle
    Domain* domain = getDPMBase()->domainHandler.getCurrentDomain();
    if (domain->containsParticle(particle))
    {
        //If the particle is in the inner domain, it is not in the mpi doamin
        if (domain->isInInnerDomain(particle))
        {
            isInMPIDomain = false;
            isMPIParticle = false;
        }
        else
        {
            isInMPIDomain = true;
            isMPIParticle = false;
        }
    }
    else
    {
        //If the particle is outside the communication zone it is not in the mpi domain
        if (domain->isInGreaterDomain(particle))
        {
            isInMPIDomain = true;
            isMPIParticle = true;
        }
        else
        {
            isInMPIDomain = false;
            isMPIParticle = true;
        }
    }
}

/*!
 * \details Mixed particles that are in periodic and mpi domains, PM-particles can create
 * ghosts of themselves PGM and PMG and additionally PGMG-particles. These last particles
 * can't be created both by the periodic routines and the paralle routines and therefore
 * the design choice is that these fully depend on the periodic routines. The update of their
 * MPI flags such as isInMPIDomain and isMPIParticle is handled by this function.
 * \param[in] particle A particle that receives the correct MPI flags
 */
void PeriodicBoundaryHandler::setMPIFlags(BaseParticle* particle)
{
    bool isInMPIDomain;
    bool isMPIParticle;
    getMPIFlags(particle, isInMPIDomain, isMPIParticle);
    particle->setInMPIDomain(isInMPIDomain);
    particle->setMPIParticle(isMPIParticle);
}

/*!
 * \details Given a real particle that interacts with periodic boundaries, ghosts need to be generated.
 * Depending on the number of boundaries this particle is interacting the number of ghosts
 * differs. A recursive function makes sure that all possibilities of the periodic complexity of the
 * real particle are generated. Note that the first periodic complexity in the list returned by
 * this function is always the periodic complexity of the real particle.
 * \param[in,out] list A list containing all possible ghost periodic complexities for a given periodic complexity
 * \param[in] periodicComplexity A periodic complexity of a real particle that is used to determine all its
 * possible ghost periodic complexities
 * \param[in,out] complexity One of the possible permutations of the ghost periodic complexity
 * \param[in] level Indicator at which level the recursive function is
 */
void
PeriodicBoundaryHandler::generateGhosts(std::vector<std::vector<int> >& list, const std::vector<int> periodicComplexity,
                                        std::vector<int>& complexity, int level)
{
    //Add the complexity to the list
    if (level == 0)
    {
        //logger(VERBOSE,"complexity: [%, %, %]",complexity[0], complexity[1], complexity[2]);
        list.push_back(complexity);
    }
    else
    {
        //Vary complexity values of level
        if (periodicComplexity[level - 1] == 1)
        {
            complexity[level - 1] = 1;
            generateGhosts(list, periodicComplexity, complexity, level - 1);
            
            complexity[level - 1] = -1;
            generateGhosts(list, periodicComplexity, complexity, level - 1);
        }
        else
        {
            generateGhosts(list, periodicComplexity, complexity, level - 1);
        }
    }
}

/*!
 * \details Sending ghosts to other processors is done in an non-blocking communication way
 * and hence the data needs to be stored somewhere. This function collects and stores the
 * data in a convenient format.
 */
void PeriodicBoundaryHandler::collectGhostParticleData()
{
    //For all new target lists
    unsigned long numberOfTargets = sendTargetList_.size();
    periodicGhostParticleSend_.resize(numberOfTargets);
    periodicGhostComplexitySend_.resize(numberOfTargets);
    for (int i = 0; i < numberOfTargets; i++)
    {
        if (sendTargetList_[i] != PROCESSOR_ID)
        {
            //For all new ghost particles copy particle data and periodic complexity
            for (int j = 0; j < numberOfNewPeriodicGhostParticlesSend_[i]; j++)
            {
                MpiPeriodicParticleID* ppid = newPeriodicParticleList_[sendTargetList_[i]][j];
                BaseParticle* particle = ppid->particle;
                periodicGhostParticleSend_[i].push_back(copyDataFromParticleToMPIParticle(particle));
                for (int k = 0; k < getSize(); k++)
                {
                    periodicGhostComplexitySend_[i].push_back(ppid->targetPeriodicComplexity[k]);
                }
            }
        }
    }
}

/*!
 * \details Not only particles need to be copied to other ghosts, also their interactions
 * with their surroundings, because these interactions contain history parameters. This
 * function collects the basic interaction data and history data into a useful
 * MPI structure
 */
void PeriodicBoundaryHandler::collectInteractionData()
{
    //For all send target lists
    InteractionHandler& iH = getDPMBase()->interactionHandler;
    unsigned long numberOfTargets = sendTargetList_.size();
    interactionDataSend_.resize(numberOfTargets);
    for (int i = 0; i < numberOfTargets; i++)
    {
        //Allocate enough space
        interactionDataSend_[i] = iH.createMPIInteractionDataArray(numberOfNewInteractionsSend_[i]);
        
        //Fill in the vector
        unsigned int indexInteraction = 0;
        for (BaseInteraction* interaction : newInteractionList_[i])
        {
            interaction->getMPIInteraction(interactionDataSend_[i], indexInteraction);
            indexInteraction++;
        }
    }
}

/*!
 * \details When processors have received new ghost particle data they have to be processed.
 * First the periodic complexity of the ghost is determined. Secondly the ghost is
 * created and given the correct flags and position and is added to the simulation.
 * Thirdly a periodicGhost ID is created for book keeping purposes.
 * A newParticle list is tracking the newly added ghost particles for a later step,
 * processing the interactions.
 * \param[in] targetIndex The index in the receiveTargetList which indicates where the data is coming from
 * \param[in] newParticles A list to which new ghost particles are added.
 */
void
PeriodicBoundaryHandler::processReceivedGhostParticleData(int targetIndex, std::vector<BaseParticle*>& newParticles)
{
    for (int j = 0; j < numberOfNewPeriodicGhostParticlesReceive_[targetIndex]; j++)
    {
        //Create the ghost particle and copy basic information
        BaseParticle* particle = MPIParticle::newParticle();
        copyDataFromMPIParticleToParticle(&(periodicGhostParticleReceive_[targetIndex][j]),
                                          particle, &(getDPMBase()->particleHandler));
        
        //Obtain real periodic complexity
        std::vector<int> realPeriodicComplexity = computePeriodicComplexity(particle->getPosition());
        
        //Obtain and set the ghost periodic complexity
        std::vector<int> ghostPeriodicComplexity(getSize());
        for (int k = 0; k < getSize(); k++)
        {
            ghostPeriodicComplexity[k] = periodicGhostComplexityReceive_[targetIndex][getSize() * j + k];
        }
        particle->setPeriodicComplexity(ghostPeriodicComplexity);
        
        //Shift the ghost particle to it's correct positions and velocities
        shiftParticle(particle, ghostPeriodicComplexity);
        
        //Add particle to simulation
        logger(VERBOSE, "Adding a ghost at position %", particle->getPosition());
        getDPMBase()->particleHandler.addGhostObject(particle);
        
        //Set the correct flags
        BaseParticle* pGhost = getDPMBase()->particleHandler.getLastObject();
        pGhost->setPeriodicGhostParticle(true);
        pGhost->setInPeriodicDomain(true);
        //Give the correct mpi flags 
        setMPIFlags(pGhost);
        
        //Create the periodic ID
        MpiPeriodicGhostParticleID* gpid = new MpiPeriodicGhostParticleID;
        gpid->realPeriodicComplexity = realPeriodicComplexity;
        gpid->particle = pGhost;
        periodicGhostList_[receiveTargetList_[targetIndex]].push_back(gpid);
        
        //Add to the newParticle list used for possible interactions
        newParticles.push_back(pGhost);
    }
}

/*!
 * \details Newly added ghost particles also need their interactions. This function adds
 * the interactions to the correct ghost particles.
 * The interaction data send over MPI only contains the ID's of the particles so a search has
 * to be performed to obtain the particle pointers. For the ghost particles we use the newParticle 
 * list that contains all newly added ghost particles. For the other particle an hgrid search is performed.
 * Finally when the two interactables are found the data is copied from the data vector.
 * \param[in] targetIndex The index in the receiveTargetList which indicates where the data is coming from
 * \param[in] newParticles A list containing the newly added ghost particles
 */
void PeriodicBoundaryHandler::processReceivedInteractionData(int targetIndex, std::vector<BaseParticle*>& newParticles)
{
    InteractionHandler& iH = getDPMBase()->interactionHandler;
    for (int l = 0; l < numberOfNewInteractionsReceive_[targetIndex]; l++)
    {
        unsigned int identificationP;
        unsigned int identificationI;
        bool isWallInteraction;
        unsigned timeStamp;
        
        //Get the general information required to setup a new interaction
        iH.getInteractionDetails(interactionDataReceive_[targetIndex], l, identificationP, identificationI,
                                 isWallInteraction, timeStamp);
        
        //Obtain the particle pointer of the ghost
        BaseParticle* pGhost;
        int idOther;
        //Check if the particle is P
        for (BaseParticle* particle : newParticles)
        {
            if (particle->getId() == identificationP)
            {
                pGhost = particle;
                idOther = identificationI;
                break;
            }
            
            if (particle->getId() == identificationI)
            {
                pGhost = particle;
                idOther = identificationP;
                break;
            }
        }
        
        //If it is a wall interaction, do stuff
        if (isWallInteraction)
        {
            BaseInteractable* I = getDPMBase()->wallHandler.getObjectById(identificationI);
            //Create interactions
            BaseInteraction* i = I->getInteractionWith(pGhost, timeStamp, &iH);
            if (i!= nullptr) i->setMPIInteraction(interactionDataReceive_[targetIndex], l, false);
        }
        else
        {
            //Obtain potential interaction particles
            std::vector<BaseParticle*> interactingParticleList;
            getDPMBase()->hGridGetInteractingParticleList(pGhost, interactingParticleList);
            
            //Find the other interacting particle
            BaseParticle* otherParticle;
            for (BaseParticle* p2 : interactingParticleList)
            {
                if (p2->getId() == idOther)
                {
                    otherParticle = p2;
                    break;
                }
            }
            
            //Add the interaction
            BaseInteraction* i = pGhost->getInteractionWith(otherParticle, timeStamp, &iH);
            if (i!= nullptr) i->setMPIInteraction(interactionDataReceive_[targetIndex], l, false);
        }
    }
}

/*!
 * \details Newly added ghost particles also need their interactions. This function adds
 * the interactions to the correct ghost particles.
 * The interaction data send over MPI only contains the ID's of the particles so a search has
 * to be performed to obtain the particle pointers. For the ghost particles we use the newParticle 
 * list that contains all newly added ghost particles. For the other particle an hgrid search is performed.
 * Finally when the two interactables are found the data is copied from the data vector.
 * Note: because these ghost particles are local, the data is actually never send.
 * \param[in] newParticles A list containing the newly added ghost particles
 */
void PeriodicBoundaryHandler::processLocalInteractionData(std::vector<BaseParticle*>& newParticles)
{
    InteractionHandler& iH = getDPMBase()->interactionHandler;
    for (int i = 0; i < sendTargetList_.size(); i++)
    {
        if (sendTargetList_[i] == PROCESSOR_ID)
        {
            for (int l = 0; l < numberOfNewInteractionsSend_[i]; l++)
            {
                unsigned int identificationP;
                unsigned int identificationI;
                bool isWallInteraction;
                unsigned timeStamp;
                
                //Get the general information required to setup a new interaction
                iH.getInteractionDetails(interactionDataSend_[i], l, identificationP, identificationI,
                                         isWallInteraction, timeStamp);
                
                //Obtain the particle pointer of the ghost
                BaseParticle* pGhost;
                int idOther;
                //Check if the particle is P
                for (BaseParticle* particle : newParticles)
                {
                    if (particle->getId() == identificationP)
                    {
                        pGhost = particle;
                        idOther = identificationI;
                        break;
                    }
                    
                    if (particle->getId() == identificationI)
                    {
                        pGhost = particle;
                        idOther = identificationP;
                        break;
                    }
                }
                
                //If it is a wall interaction, do stuff
                if (isWallInteraction)
                {
                    BaseWall* I = getDPMBase()->wallHandler.getObjectById(identificationI);
                    //Create interactions
                    BaseInteraction* j = I->getInteractionWith(pGhost, timeStamp, &iH);
                    if (j!= nullptr) j->setMPIInteraction(interactionDataSend_[i], l, false);
                    logger(VERBOSE, "Wall interaction added!");
                }
                else
                {
                    //Obtain potential interaction particles
                    std::vector<BaseParticle*> interactingParticleList;
                    getDPMBase()->hGridGetInteractingParticleList(pGhost, interactingParticleList);
                    
                    if (interactingParticleList.empty())
                    {
                        logger(VERBOSE, "Failed in creating an interaction :(");
                    }
                    
                    //Find the other interacting particle
                    BaseParticle* otherParticle = nullptr;
                    for (BaseParticle* p2 : interactingParticleList)
                    {
                        if (p2->getId() == idOther)
                        {
                            otherParticle = p2;
                            break;
                        }
                    }
                    if (otherParticle == nullptr)
                    {
                        //The interacting object can't be found
                        continue;
                    }
                    
                    //Add the interaction
                    BaseInteraction* j = pGhost->getInteractionWith(otherParticle, timeStamp, &iH);
                    if (j!= nullptr) j->setMPIInteraction(interactionDataSend_[i], l, false);
                }
            }
        }
    }
}


/*!
 * \details Periodic particles that have a copy somewhere else create an ID for book keeping
 * This function also flags the particle being in the periodic domain such that it is ignored
 * when finding new periodic particles.
 */
void PeriodicBoundaryHandler::processPeriodicParticles()
{
    for (int i : sendTargetList_)
    {
        for (int j = 0; j < newPeriodicParticleList_[i].size(); j++)
        {
            //Update the particle status
            MpiPeriodicParticleID* ppid = newPeriodicParticleList_[i][j];
            ppid->particle->setInPeriodicDomain(true);
            
            //make new entry in the list
            periodicParticleList_[i].push_back(ppid);
        }
    }
}

/*!
 * \details Ghost particles sometimes are located on the same processor as the original periodic particle.
 * These particles are processed in this function. Firstly the ghost is created with the correct position,
 * periodic complexity and flags and is added to the particleHandler. Secondly ID's are made for both
 * ghost and periodic particles to keep track of their movements. Periodic particles are flaged now
 * as being in the periodic domain.
 */
void PeriodicBoundaryHandler::processLocalGhostParticles(std::vector<BaseParticle*>& newParticles)
{
    for (int i = 0; i < sendTargetList_.size(); i++)
    {
        //Check if the processor is sending ghosts to itself
        if (sendTargetList_[i] == PROCESSOR_ID)
        {
            for (int j = 0; j < numberOfNewPeriodicGhostParticlesSend_[i]; j++)
            {
                MpiPeriodicParticleID* ppid = newPeriodicParticleList_[PROCESSOR_ID][j];
                BaseParticle* particle = ppid->particle;
                
                //Create ghost particle
                BaseParticle* pGhost = particle->copy();
                
                //Obtain and set the ghost periodic complexity
                pGhost->setPeriodicComplexity(ppid->targetPeriodicComplexity);
                
                //Shift the particle
                shiftParticle(pGhost);
                
                //Set the correct flags
                pGhost->setPeriodicGhostParticle(true);
                pGhost->setInPeriodicDomain(true);
                //Note: mpi flags dont have to be set for local particles
                //these flags are copied from the original particle
                logger(VERBOSE, "Adding a ghost with id % at position %", pGhost->getId(), pGhost->getPosition());
                
                //Finally add it to the particle handler
                getDPMBase()->particleHandler.addGhostObject(pGhost);
                
                //Do some bookkeeping
                MpiPeriodicGhostParticleID* gpid = new MpiPeriodicGhostParticleID;
                gpid->particle = pGhost;
                gpid->otherParticle = particle;
                gpid->realPeriodicComplexity = particle->getPeriodicComplexity();
                periodicGhostList_[PROCESSOR_ID].push_back(gpid);
                
                //Flag real particle and do some bookkeeping
                particle->setInPeriodicDomain(true);
                ppid->otherParticle = pGhost;
                
                //Add to the new particle list for an interaction update
                newParticles.push_back(pGhost);
            }
        }
    }
}

/*!
 * \details This function updates the position/velocity and periodic complexity of ghost particles.
 * There are two versions of update, a local and a global one. The local one can update the particles
 * through pointers, as the periodic particle is located ont he same processor as the ghost.
 * The local version has received data from other processors which is stored in the corresponding vectors
 */
void PeriodicBoundaryHandler::updateParticles()
{
    MPIContainer& communicator = MPIContainer::Instance();
    
    //For all lists that contain ghost particles
    //The variable dataIndex indicates which index in update<...>Receive_ the data is located
    int dataIndex = -1;
    for (int i = 0; i < NUMBER_OF_PROCESSORS; i++)
    {
        unsigned long numberOfParticles = periodicGhostList_[i].size();
        bool global;
        if (numberOfParticles > 0)
        {
            //Check if the update is global or from the current processor
            if (i != PROCESSOR_ID)
            {
                global = true;
                dataIndex++;
            }
            else
            {
                global = false;
            }
            
            //Update all particles
            for (int p = 0; p < numberOfParticles; p++)
            {
                MpiPeriodicGhostParticleID* pgid = periodicGhostList_[i][p];
                BaseParticle* pGhost = pgid->particle;
                
                //Depending on where the particle is located the data is stored in a different place
                if (global)
                {
                    //logger(VERBOSE,"i: %, numberOfParticles: %, dataIndex: %",i,numberOfParticles, dataIndex);
                    //Check if this particle really belongs to the data that is send
                    logger.assert(pGhost->getId() == updatePositionDataReceive_[dataIndex][p].id,
                                  "Periodic particle lists are not in syc");
                    
                    //Add the real position and velocity of the particles
                    //Note: The before updating the position is the position of the ghost
                    //Note: The received position and velocity values are of the real particle
                    //Note: It will be shifted to the correct values after the complexity is computed
                    pGhost->setPreviousPosition(pGhost->getPosition());
                    pGhost->setPosition(updatePositionDataReceive_[dataIndex][p].position);
                    pGhost->setOrientation(updatePositionDataReceive_[dataIndex][p].orientation);
                    pGhost->setVelocity(updateVelocityDataReceive_[dataIndex][p].velocity);
                    pGhost->setAngularVelocity(updateVelocityDataReceive_[dataIndex][p].angularVelocity);
                    
                }
                else
                {
                    //MpiPeriodicGhostParticleID * pgip = periodicGhostList_[i][p];
                    pGhost->setPreviousPosition(pGhost->getPosition());
                    pGhost->setPosition(pgid->otherParticle->getPosition());
                    pGhost->setOrientation(pgid->otherParticle->getOrientation());
                    pGhost->setVelocity(pgid->otherParticle->getVelocity());
                    pGhost->setAngularVelocity(pgid->otherParticle->getAngularVelocity());
                }
                
                //Move current periodic complexity to previous periodic complexity
                pGhost->setPreviousPeriodicComplexity(pGhost->getPeriodicComplexity());
                
                //Compute the new realPeriodicComplexity
                std::vector<int> previousRealPeriodicComplexity = periodicGhostList_[i][p]->realPeriodicComplexity;
                //The ghost particle has the real position at the moment, so compute the real periodic complexity here
                std::vector<int> realPeriodicComplexity = computePeriodicComplexity(pGhost->getPosition());
                std::vector<int> periodicComplexity(getSize());
                for (int b = 0; b < getSize(); b++)
                {
                    int sign = mathsFunc::sign(previousRealPeriodicComplexity[b]
                                               * realPeriodicComplexity[b]
                                               * pGhost->getPreviousPeriodicComplexity()[b]);
                    periodicComplexity[b] = sign * abs(realPeriodicComplexity[b]);
                    //The maser boundary needs this correction
                    if (periodicComplexity[b] == -3)
                    {
                        periodicComplexity[b] = 1;
                    }
                    
                }
                pGhost->setPeriodicComplexity(periodicComplexity);
                
                for (int b = 0; b < getSize(); b++)
                {
                    objects_[b]->modifyGhostAfterCreation(pGhost, b);
                }
                //Shift the particle to correct position
                //Note: If the real particle changed complexity this position will be calculated incorrectly
                //Hence the previous complexity is used.
                shiftParticle(pGhost, pGhost->getPreviousPeriodicComplexity());
                
                //Update hGrid
                Vec3D displacement = pGhost->getPreviousPosition() - pGhost->getPosition();
                getDPMBase()->hGridUpdateMove(pGhost, displacement.getLengthSquared());
                
                //Do some book keeping
                periodicGhostList_[i][p]->realPeriodicComplexity = realPeriodicComplexity;
                
            }
        }
    }
    
    //Update periodic complexity of periodic particles
    for (BaseParticle* particle : getDPMBase()->particleHandler)
    {
        if ((particle->isInPeriodicDomain()) && !(particle->isPeriodicGhostParticle()))
        {
            //Update the periodicComplexity of the real particles
            particle->setPreviousPeriodicComplexity(particle->getPeriodicComplexity());
            std::vector<int> periodicComplexity;
            int totalPeriodicComplexity;
            computePeriodicComplexity(periodicComplexity, totalPeriodicComplexity, particle->getPosition());
            //Modify periodic complexity tailored to specific boundary requirements
            for (int b = 0; b < getSize(); b++)
            {
                objects_[b]->modifyPeriodicComplexity(periodicComplexity, totalPeriodicComplexity, particle, b);
            }
            particle->setPeriodicComplexity(periodicComplexity);
        }
    }
}


/*!
 * \detail An important distinction between periodic particles is if they are real or not. This check is
 * generally required for the particle status update. A real particle only has positive periodic
 * complexity values. So the moment a negative value is found it is clear the particle is not real.
 * \param[in] complexity The periodic complexity that indicates the status of the periodic particle
 * \return When the particle is real this function returns true, otherwise it will return false.
 */
bool PeriodicBoundaryHandler::checkIsReal(const std::vector<int> complexity)
{
    bool isReal = true;
    for (int i = 0; i < getSize(); i++)
    {
        if (mathsFunc::sign(complexity[i]) == -1)
        {
            isReal = false;
        }
    }
    
    return isReal;
}

/*!
 * \detail A large part of the status update requires on if the periodic complexity of the real 
 *  particle is changed. This function will check if the previous and current periodic complexity
 *  differ. If this is the case it returns true, if they remain the same the function returns false
 * \param[in] previousPeriodicComplexity A periodioc complexity vector that is used as reference
 * \param[in] currentPeriodicComplexity A periodic complexity that is checked for difference against
 * a reference periodic complexity.
 * \return True if the periodic complexity as changed with respect to the reference periodic complexity.
 * returns false if the periodic complexity is exactly the same.
 */
bool PeriodicBoundaryHandler::checkChanged(const std::vector<int> previousComplexity,
                                           const std::vector<int> currentComplexity)
{
    bool changed = false;
    
    for (int i = 0; i < currentComplexity.size(); i++)
    {
        if (previousComplexity[i] != currentComplexity[i])
        {
            changed = true;
        }
    }
    
    return changed;
}

/*!
 * \details This function updates the status of the particles based on the periodic complexity of the particles.
 * this is beneficial as no round-off errors are made due to the shift in position. If a particle changes
 * it's periodic complexity it is either removed or turned into a real particle. The real particle will
 * be re-introduced in a later step. Particles that need to be deleted will be stored in a vector as 
 * these particles might also need to be flushed from the Domain.h lists.
 * \param[in,out] particlesToBeDeleted List of particles that will need to be removed from the simulation
 */
void PeriodicBoundaryHandler::updateParticleStatus(std::set<BaseParticle*>& particlesToBeDeleted)
{
    MPIContainer& communicator = MPIContainer::Instance();
    int processorID = communicator.getProcessorID();
    int numberOfProcessors = communicator.getNumberOfProcessors();
    std::set<MpiPeriodicParticleID*> deletePeriodicIDList;
    std::set<MpiPeriodicGhostParticleID*> deletePeriodicGhostIDList;
    std::set<BaseParticle*> specialUpdateParticleList;
    
    //For all domains
    for (int i = 0; i < numberOfProcessors; i++)
    {
        int numberOfPeriodicParticles = periodicParticleList_[i].size();
        int numberOfPeriodicGhostParticles = periodicGhostList_[i].size();
        
        //Loop over all periodic particles to see if their complexity changed 
        for (int p = 0; p < numberOfPeriodicParticles; p++)
        {
            MpiPeriodicParticleID* ppid = periodicParticleList_[i][p];
            BaseParticle* particle = ppid->particle;
            
            //Check particle status
            bool isReal = checkIsReal(particle->getPeriodicComplexity());
            bool changed = checkChanged(particle->getPreviousPeriodicComplexity(), particle->getPeriodicComplexity());
            
            //Only if the particle changed we need to undertake action
            if (changed)
            {
                if (isReal)
                {
                    //Flag this particle as normal, it will be re-introduced when finding new periodic particles
                    logger(VERBOSE, "Real particle % changed complexity at: %", particle->getId(),
                           particle->getPosition());
                    particle->setInPeriodicDomain(false);
                    //Incase of a special flag 3, perform update action
                    updateMaserParticle(particle);
                }
                else
                {
                    //Oh noes, the particle became a ghost. Kill it with balefire!!... if it is necessary
                    logger(VERBOSE, "Real particle % changed to ghost at: %", particle->getId(),
                           particle->getPosition());
                    particlesToBeDeleted.insert(particle);
                }
                
                //Delete the ID
                deletePeriodicIDList.insert(ppid);
                periodicParticleList_[i][p] = nullptr;
            }
            
            
            //If a PM particle changes from to PMG it will need to be deleted.
            //The deletion will be done by the M boundary, but it still needs to be flushed from the periodic lists
            if (particle->isInMPIDomain())
            {
                //Store old values and compute new to see if the status has changed
                bool isMPIParticleOld = particle->isMPIParticle();
                bool isMPIParticle;
                bool isInMPIDomain;
                getMPIFlags(particle, isInMPIDomain, isMPIParticle);
                
                //Particle needs to be removed from lists if it becomes an MG particle
                if (isMPIParticleOld != isMPIParticle)
                {
                    logger(VERBOSE, "PM to PMG: Flush from boundary");
                    deletePeriodicIDList.insert(ppid);
                    periodicParticleList_[i][p] = nullptr;
                }
                
                //MG particle needs to be removed from lists if it moves to another domain
                if (!isInMPIDomain)
                {
                    if (particle->isMPIParticle())
                    {
                        logger(VERBOSE, "PMG leaves domain: Flush from boundary");
                        deletePeriodicIDList.insert(ppid);
                        periodicParticleList_[i][p] = nullptr;
                    }
                }
            }
        }
        
        //Loop over all periodic ghost particles to see if their complexity changed
        for (int p = 0; p < numberOfPeriodicGhostParticles; p++)
        {
            MpiPeriodicGhostParticleID* pgid = periodicGhostList_[i][p];
            BaseParticle* pGhost = pgid->particle;
            
            //Check particle status
            bool isReal = checkIsReal(pGhost->getPeriodicComplexity());
            bool changed = checkChanged(pGhost->getPreviousPeriodicComplexity(), pGhost->getPeriodicComplexity());
            
            //Update mixed particles, particles that also are in the mpi domain also need an update
            //Note that these particles are not listed in the domain lists, they are taken care here.
            //Store old values and update status of pGhost
            bool isInMPIDomainOld = pGhost->isInMPIDomain();
            bool isMPIParticleOld = pGhost->isMPIParticle();
            setMPIFlags(pGhost);
            if (isInMPIDomainOld)
            {
                
                //Case 1: pGhost changed from real to mpi particle
                if (isMPIParticleOld != pGhost->isMPIParticle())
                {
                    //Case 1: turned from M ghost to M
                    //The correct flags have been set above already
                    
                    //Case 2: Turned from M to M ghost
                    if (pGhost->isMPIParticle())
                    {
                        logger(VERBOSE, "PGM to PGMG: Deleting particle.");
                        particlesToBeDeleted.insert(pGhost);
                        deletePeriodicGhostIDList.insert(pgid);
                        periodicGhostList_[i][p] = nullptr;
                    }
                    
                }
                
                //Case 2: pGhost left the mpi domain
                if (pGhost->isInMPIDomain() != isInMPIDomainOld)
                {
                    //Case 1: Moved inside the current domain
                    //The correct flags have been set above already
                    
                    //Case 2: Moved to a neighbour domain
                    if (pGhost->isMPIParticle())
                    {
                        logger(VERBOSE, "PGMG moved out of domain: Deleting particle.");
                        particlesToBeDeleted.insert(pGhost);
                        deletePeriodicGhostIDList.insert(pgid);
                        periodicGhostList_[i][p] = nullptr;
                    }
                }
            }
            
            //Check if the particles need to be deleted based on their periodic complexity
            if (changed)
            {
                if (isReal)
                {
                    //Check if the complexity of the particle is truely real based on it's current position
                    std::vector<int> pc = computePeriodicComplexity(pGhost->getPosition());
                    int tpc = 0;
                    for (int b = 0; b < getSize(); b++)
                    {
                        objects_[b]->modifyPeriodicComplexity(pc, tpc, pGhost, b);
                    }
                    if (!checkIsReal(pc))
                    {
                        logger(ERROR, "Round-off error detected.");
                        //logger(WARN,"Round-off error corrected, phew!");
                    }
                    
                    //There are two cases
                    //Case 1: PGMG particles turning PMG; These are still not real
                    //Case 2: PGM particles turning PG; These are truely real
                    if (pGhost->isMPIParticle())
                    {
                        logger(VERBOSE, "PGMG turned PMG: delete it");
                        particlesToBeDeleted.insert(pGhost);
                    }
                    else
                    {
                        //Turn the particle real
                        logger(VERBOSE, "Ghost particle changed to real at position: %", pGhost->getPosition());
                        pGhost->setInPeriodicDomain(false);
                        pGhost->setPeriodicGhostParticle(false);
                        //Make sure this particle can be detected now by the parallel boundaries
                        pGhost->setInMPIDomain(false);
                    }
                }
                else
                {
                    //Delete the particle
                    logger(VERBOSE, "Ghost particle changed complexity at position: %", pGhost->getPosition());
                    particlesToBeDeleted.insert(pGhost);
                }
                
                //Delete the ID
                deletePeriodicGhostIDList.insert(pgid);
                periodicGhostList_[i][p] = nullptr;
            }
        }
        
    }
    
    //Delete IDs 
    for (auto ppid_it : deletePeriodicIDList)
    {
        delete ppid_it;
    }
    
    for (auto pgid_it : deletePeriodicGhostIDList)
    {
        delete pgid_it;
    }
    unsigned int nextId = getDPMBase()->particleHandler.getNextId();
    communicator.broadcast(nextId);
    getDPMBase()->particleHandler.setNextId(nextId);
}

/*!
 * \details A given complexity can be turned into a particle position, somewhere in the simulation domain.
 * This position belongs to a specific processor and this function computes which processor this actually is.
 * Note: When an PMG particle is in the periodic domain, it needs to be copied into
 * a PG"MG" particle. A particle that is in the MG domain, but is not updated by 
 * the M boundaries. The position of this MG particle is not in the actual domain to
 * where it has to be copied, but slightly outside.
 * For that reason the position is shifted to the middle of the domain the official particle
 * is found on, and that position is used to determine the actual domain the particle
 * has to go to. Note that currently this only works for a structured grid with 
 * equal domain sizes
 * \param[in] complexity A complexity vector which determines the location of a ghost particle
 * \return The processor ID of the ghost particle location
 */
int PeriodicBoundaryHandler::findTargetProcessor(const std::vector<int>& complexity)
{
    //Find the middle of this domain (or any valid point anyway)
    Vec3D middlePosition = getDPMBase()->domainHandler.getCurrentDomain()->getMiddle();
    
    //Create the particle with a target position
    SphericalParticle particle;
    particle.setPosition(middlePosition);
    shiftParticle(&particle, complexity);
    
    //Obtain target domain
    int targetGlobalIndex = getDPMBase()->domainHandler.getParticleDomainGlobalIndex(&particle);
    return getDPMBase()->domainHandler.getParticleProcessor(targetGlobalIndex);
}

/*!
 * \details Particles that are located close to a periodic boundary, but not flagged as being in
 * the periodic domain will be added to a list by this function. If they indeed have to be added
 * a periodic ID will be created that stores some basic information required to setup a ghost particle.
 * \param[in] particle A particle that is considered to be added to the periodic boundary lists
 */
void PeriodicBoundaryHandler::findNewParticle(BaseParticle* particle)
{
    //Skip the ghost particles
    if (checkIfAddNewParticle(particle))
    {
        //If the particle was not yet flagged to be in the periodic domain, check if it is in the domain
        //if(!particle->isInPeriodicDomain())
        {
            int totalPeriodicComplexity;
            std::vector<int> periodicComplexity;
            computePeriodicComplexity(periodicComplexity, totalPeriodicComplexity, particle->getPosition());
            
            //Modify periodic complexity tailored to specific boundary requirements
            for (int b = 0; b < getSize(); b++)
            {
                objects_[b]->modifyPeriodicComplexity(periodicComplexity, totalPeriodicComplexity, particle, b);
            }
            
            //Set periodicComplexity
            particle->setPeriodicComplexity(periodicComplexity);
            
            //Create ghost particle ID's
            std::vector<std::vector<int> > list(0);
            if (totalPeriodicComplexity > 0)
            {
                periodicComplexity = particle->getPeriodicComplexity();
                
                //Generating all possible complexities.
                generateGhosts(list, periodicComplexity, periodicComplexity, getSize());
                //logger(VERBOSE,"New ghost particles: %",list.size() - 1);
                
                //Note: the first one in the list is the real particle such that we skip it
                for (int i = 1; i < list.size(); i++)
                {
                    //Compute target domain
                    //logger(VERBOSE,"Particle position: %",particle->getPosition());
                    int targetProcessor = findTargetProcessor(list[i]);
                    
                    //Create periodic particle ID
                    MpiPeriodicParticleID* ppid = new MpiPeriodicParticleID;
                    ppid->particle = particle;
                    ppid->periodicComplexity = periodicComplexity;
                    ppid->targetPeriodicComplexity = list[i];
                    ppid->targetProcessor = targetProcessor;
                    newPeriodicParticleList_[targetProcessor].push_back(ppid);
                    logger(VERBOSE, "Adding a periodic particle with id % and real particle position: %",
                           particle->getId(), particle->getPosition());
                }
            }
        }
    }
}

/*!
 * \details Loops over all base particles in the domain to see if they have to be added to the
 * periodic domain lists
 */
void PeriodicBoundaryHandler::findNewParticles()
{
    //Check for all particles if there are unassigned periodic particles
    for (auto particle_it = getDPMBase()->particleHandler.begin();
         particle_it != getDPMBase()->particleHandler.end(); particle_it++)
    {
        findNewParticle(*particle_it);
    }
    
}

/*!
 * \details Ghost particles have interactions with history parameters. These interactions also need
 * to be copied to the other processor. This function finds these interactions and stores them.
 */
void PeriodicBoundaryHandler::findNewInteractions()
{
    //Check all new particle lists
    newInteractionList_.resize(sendTargetList_.size());
    for (int i = 0; i < sendTargetList_.size(); i++)
    {
        for (int p = 0; p < newPeriodicParticleList_[sendTargetList_[i]].size(); p++)
        {
            BaseParticle* particle = newPeriodicParticleList_[sendTargetList_[i]][p]->particle;
            
            //Loop over all its interactions
            std::vector<BaseInteraction*> interactions = particle->getInteractions();
            for (BaseInteraction* interaction : interactions)
            {
                //Find out what the new particle is interacting with
                BaseParticle* particleP = dynamic_cast<BaseParticle*>(interaction->getP());
                BaseParticle* objectI = dynamic_cast<BaseParticle*>(interaction->getI());
                
                //If the P in the interaction structure is the new particle, find I
                if (particle == particleP)
                {
                    //Check if the new particle is interacting with a wall
                    if (!objectI)
                    {
                        newInteractionList_[i].push_back(interaction);
                    }
                    else //is I a particle
                    {
                        newInteractionList_[i].push_back(interaction);
                    }
                }
                else //newBoundaryParticle is I in the interaction, P can only be a particle
                {
                    newInteractionList_[i].push_back(interaction);
                }
            }
        }
        
        //Store the number of interactions
        numberOfNewInteractionsSend_.push_back(newInteractionList_[i].size());
    }
}

/*!
 * \details  When communicating ghost particles from periodic boundaries in paralell, it is not clear from the
 * start which processor is sending to which processor. First all processors check to how many domains the
 * processor is sending. A second communication step is then performed to tell all other processors which
 * domains the processor is sending to.
 * After this function has completed every processor has a list of targets that the processor sends data to
 * and a list of targets the processor receives data from.
 */
void PeriodicBoundaryHandler::communicateTargetDomains()
{
#ifdef MERCURY_USE_MPI
    //Step 1: Check to how many domains this domain has to send particles
    int numberOfTargetDomainsLocal = 0;
    for (int index = 0; index < NUMBER_OF_PROCESSORS; index++) 
    {
        if (newPeriodicParticleList_[index].size() > 0)
        {
            numberOfTargetDomainsLocal++;
            sendTargetList_.push_back(index);
            numberOfNewPeriodicGhostParticlesSend_.push_back(newPeriodicParticleList_[index].size());
        }
    }

    //Communicate with other domains to how many domains this domain is sending
    std::vector<int> numberOfTargetDomains(NUMBER_OF_PROCESSORS);
    MPIContainer::Instance().allGather(numberOfTargetDomainsLocal,1,numberOfTargetDomains,1);    

    //Step 2: Communicate to other domains what the target domains are
    std::vector<int> receiveTargetList;
    for (int index = 0; index < NUMBER_OF_PROCESSORS; index++)
    {
        if (numberOfTargetDomains[index] > 0)
        {
            if (index == PROCESSOR_ID)
            {
                MPIContainer::Instance().broadcast(sendTargetList_.data(), numberOfTargetDomainsLocal, index);
            }
            else
            {
                receiveTargetList.resize(numberOfTargetDomains[index], 0);
                MPIContainer::Instance().broadcast(receiveTargetList.data(),numberOfTargetDomains[index],index);
            }
        
            //If this domain is a target, flag it.
            for (int i = 0; i < receiveTargetList.size(); i++)
            {
                if (receiveTargetList[i] == PROCESSOR_ID)
                {
                    receiveTargetList_.push_back(index);
                }
            } 
        }
    }
#endif
}

/*!
 * \details Every receiving target that will receive ghost particle data will need to know
 * how many particles are coming their way. This function performs that communication step.
 * The same for interactions, not every ghost particle receives the same number of interactions
 * and hence that step is done seperately.
 */
void PeriodicBoundaryHandler::communicateNumberOfNewParticlesAndInteractions()
{
    MPIContainer& communicator = MPIContainer::Instance();
    //Communicate to how many particles the target domains receive
    int tagSend;
    for (int i = 0; i < sendTargetList_.size(); i++)
    {
        if (sendTargetList_[i] != PROCESSOR_ID)
        {
            tagSend = PROCESSOR_ID * MAX_PROC + sendTargetList_[i] * 10 + MercuryMPITag::PARTICLE_COUNT;
            communicator.send(numberOfNewPeriodicGhostParticlesSend_[i], sendTargetList_[i], tagSend);
            
            tagSend = PROCESSOR_ID * MAX_PROC + sendTargetList_[i] * 10 + MercuryMPITag::INTERACTION_COUNT;
            communicator.send(numberOfNewInteractionsSend_[i], sendTargetList_[i], tagSend);
        }
    }
    
    //Perform the receive routines
    numberOfNewPeriodicGhostParticlesReceive_.resize(receiveTargetList_.size());
    numberOfNewInteractionsReceive_.resize(receiveTargetList_.size());
    int tagReceive;
    for (int i = 0; i < receiveTargetList_.size(); i++)
    {
        if (receiveTargetList_[i] != PROCESSOR_ID)
        {
            tagReceive = receiveTargetList_[i] * MAX_PROC + PROCESSOR_ID * 10 + MercuryMPITag::PARTICLE_COUNT;
            communicator.receive(numberOfNewPeriodicGhostParticlesReceive_[i], receiveTargetList_[i], tagReceive);
            
            tagReceive = receiveTargetList_[i] * MAX_PROC + PROCESSOR_ID * 10 + MercuryMPITag::INTERACTION_COUNT;
            communicator.receive(numberOfNewInteractionsReceive_[i], receiveTargetList_[i], tagReceive);
        }
    }
}

/*!
 * \details The first step in adding new particles is the prepare step. First the target domains
 * are determined such that every processors knows where particles are coming from. This step includes
 * finding new particles. The second step is to find the interactions that have to be send with the
 * new particles. The receiving targets need to know how many particles and interactions they receive
 * which is the third step. Finally a communication synchronisation is performed to actually perform
 * the communications which were pending in the MPIContainer.
 */
void PeriodicBoundaryHandler::prepareNewParticleTransmission()
{
    //Communicate which domains send and which domains receive particles
    communicateTargetDomains();
    
    //Find new interactions
    findNewInteractions();
    
    //Communicate the number of particles and interactions send to the other domains
    communicateNumberOfNewParticlesAndInteractions();
    
    //Synchronise the communications
    MPIContainer::Instance().sync();
}

void PeriodicBoundaryHandler::performNewParticleTransmission()
{
    MPIContainer& communicator = MPIContainer::Instance();
    int numberOfTargetsReceive = receiveTargetList_.size();
    int numberOfTargetsSend = sendTargetList_.size();
    
    //Make sure that the receiving vectors have the correct length
    periodicGhostParticleReceive_.resize(numberOfTargetsReceive);
    periodicGhostComplexityReceive_.resize(numberOfTargetsReceive);
    interactionDataReceive_.resize(numberOfTargetsReceive);
    InteractionHandler& iH = getDPMBase()->interactionHandler;
    for (int i = 0; i < numberOfTargetsReceive; i++)
    {
        periodicGhostParticleReceive_[i].resize(numberOfNewPeriodicGhostParticlesReceive_[i]);
        periodicGhostComplexityReceive_[i].resize(numberOfNewPeriodicGhostParticlesReceive_[i] * getSize());
        interactionDataReceive_[i] = iH.createMPIInteractionDataArray(numberOfNewInteractionsReceive_[i]);
    }
    
    //Collect data for sending
    collectGhostParticleData();
    collectInteractionData();
    
    //Send particle and interaction data
    for (int i = 0; i < numberOfTargetsSend; i++)
    {
        if (sendTargetList_[i] != PROCESSOR_ID)
        {
            //Send particle data
            int sendCount = numberOfNewPeriodicGhostParticlesSend_[i];
            int processor = sendTargetList_[i];
            int tagSend = PROCESSOR_ID * MAX_PROC + sendTargetList_[i] * 10 + MercuryMPITag::PARTICLE_DATA;
            communicator.send(periodicGhostParticleSend_[i].data(), MercuryMPIType::PARTICLE, sendCount, processor,
                              tagSend);
            
            //Send complexity
            tagSend = PROCESSOR_ID * MAX_PROC + sendTargetList_[i] * 10 + MercuryMPITag::PERIODIC_COMPLEXITY;
            sendCount = numberOfNewPeriodicGhostParticlesSend_[i] * getSize();
            communicator.send(periodicGhostComplexitySend_[i].data(), sendCount, processor, tagSend);
            
            //Send interactions
            tagSend = PROCESSOR_ID * MAX_PROC + sendTargetList_[i] * 10 + MercuryMPITag::INTERACTION_DATA;
            sendCount = numberOfNewInteractionsSend_[i];
            if (sendCount > 0)
            {
                communicator.send(interactionDataSend_[i], MercuryMPIType::INTERACTION, sendCount, processor, tagSend);
            }
        }
    }
    
    //Receive data
    for (int i = 0; i < numberOfTargetsReceive; i++)
    {
        if (receiveTargetList_[i] != PROCESSOR_ID)
        {
            //Receive particle data
            int receiveCount = numberOfNewPeriodicGhostParticlesReceive_[i];
            int processor = receiveTargetList_[i];
            int tagReceive = receiveTargetList_[i] * MAX_PROC + PROCESSOR_ID * 10 + MercuryMPITag::PARTICLE_DATA;
            communicator.receive(periodicGhostParticleReceive_[i].data(), MercuryMPIType::PARTICLE,
                                 receiveCount, processor, tagReceive);
            
            //Receive complexity
            tagReceive = receiveTargetList_[i] * MAX_PROC + PROCESSOR_ID * 10 + MercuryMPITag::PERIODIC_COMPLEXITY;
            receiveCount = numberOfNewPeriodicGhostParticlesReceive_[i] * getSize();
            communicator.receive(periodicGhostComplexityReceive_[i].data(), receiveCount, processor, tagReceive);
            
            //Send interactions
            tagReceive = receiveTargetList_[i] * MAX_PROC + PROCESSOR_ID * 10 + MercuryMPITag::INTERACTION_DATA;
            receiveCount = numberOfNewInteractionsReceive_[i];
            if (receiveCount > 0)
            {
                communicator.receive(interactionDataReceive_[i], MercuryMPIType::INTERACTION, receiveCount, processor,
                                     tagReceive);
            }
        }
    }
    
    //Synchronise the communications
    communicator.sync();
}

void PeriodicBoundaryHandler::finaliseNewParticleTransmission()
{
    //Process the global received data
    for (int i = 0; i < receiveTargetList_.size(); i++)
    {
        std::vector<BaseParticle*> newGhostParticles;
        processReceivedGhostParticleData(i, newGhostParticles);
        processReceivedInteractionData(i, newGhostParticles);
    }
    
    //Process the local data
    std::vector<BaseParticle*> newGhostParticles;
    processLocalGhostParticles(newGhostParticles);
    processLocalInteractionData(newGhostParticles);
    
    //Process the periodic particles;
    processPeriodicParticles();
    
    //Clear lists
    receiveTargetList_.clear();
    sendTargetList_.clear();
    numberOfNewPeriodicGhostParticlesReceive_.clear();
    numberOfNewPeriodicGhostParticlesSend_.clear();
    numberOfNewInteractionsSend_.clear();
    numberOfNewInteractionsReceive_.clear();
    periodicGhostParticleReceive_.clear();
    periodicGhostParticleSend_.clear();
    periodicGhostComplexityReceive_.clear();
    periodicGhostComplexitySend_.clear();
    interactionDataSend_.clear();
    interactionDataReceive_.clear();
    newInteractionList_.clear();
    for (auto& i : newPeriodicParticleList_)
    {
        i.clear();
    }
    for (auto& i : newPeriodicParticleList_)
    {
        i.clear();
    }
}

void PeriodicBoundaryHandler::preparePositionAndVelocityUpdate()
{
    MPIContainer& communicator = MPIContainer::Instance();
    unsigned int numberOfProcessors = communicator.getNumberOfProcessors();
    unsigned int processorID = communicator.getProcessorID();
    
    //For all lists that contain periodic particles
    for (unsigned int i = 0; i < numberOfProcessors; i++)
    {
        unsigned int numberOfParticles = periodicParticleList_[i].size();
        if (numberOfParticles > 0 && i != processorID)
        {
            //Increase the vector size;
            updatePositionDataSend_.emplace_back(0);
            updateVelocityDataSend_.emplace_back(0);
            
            //Collect the data
            for (unsigned int p = 0; p < numberOfParticles; p++)
            {
                BaseParticle* particle = periodicParticleList_[i][p]->particle;
                updatePositionDataSend_.back().push_back(copyPositionFrom(particle));
                updateVelocityDataSend_.back().push_back(copyVelocityFrom(particle));
            }
            
            //Send position data
            unsigned int count = numberOfParticles;
            unsigned int processor = i;
            unsigned int tag = processorID * MAX_PROC + processor * 10 + MercuryMPITag::POSITION_DATA;
            communicator.send(updatePositionDataSend_.back().data(), MercuryMPIType::POSITION, count, processor, tag);
            
            //Send velocity data
            tag = processorID * MAX_PROC + processor * 10 + MercuryMPITag::VELOCITY_DATA;
            communicator.send(updateVelocityDataSend_.back().data(), MercuryMPIType::VELOCITY, count, processor, tag);
        }
    }
    
    //For all lists that contain ghost particles
    for (int i = 0; i < numberOfProcessors; i++)
    {
        unsigned int numberOfParticles = periodicGhostList_[i].size();
        if (numberOfParticles > 0 && i != processorID)
        {
            //Increase the vector size
            updatePositionDataReceive_.emplace_back(numberOfParticles);
            updateVelocityDataReceive_.emplace_back(numberOfParticles);
            
            //Receive position data
            int count = numberOfParticles;
            int processor = i;
            int tag = processor * MAX_PROC + processorID * 10 + MercuryMPITag::POSITION_DATA;
            communicator.receive(updatePositionDataReceive_.back().data(), MercuryMPIType::POSITION, count, processor,
                                 tag);
            
            //Receive velocity data
            tag = processor * MAX_PROC + processorID * 10 + MercuryMPITag::VELOCITY_DATA;
            communicator.receive(updateVelocityDataReceive_.back().data(), MercuryMPIType::POSITION, count, processor,
                                 tag);
        }
    }
    
    //Synchronise all the requests
    communicator.sync();
}

void PeriodicBoundaryHandler::finalisePositionAndVelocityUpdate()
{
    //Update the positions, velocities and periodic complexities of particles
    updateParticles();
    
    //Delete vectors
    updatePositionDataSend_.clear();
    updatePositionDataReceive_.clear();
    updateVelocityDataSend_.clear();
    updateVelocityDataReceive_.clear();
}

//TODO this is for the deletion boundary?
void PeriodicBoundaryHandler::flushParticles(std::set<BaseParticle*>& particlesToBeFlushed)
{
    /*!
     * \todo This function uses a brute force method to flush the particles, this
     * can be done in a smarter way if profiling shows that this is a bottleneck.
     */
    std::set<MpiPeriodicParticleIDBase*> toBeDeleted;
    for (auto p_it : particlesToBeFlushed)
    {
        for (int i = 0; i < NUMBER_OF_PROCESSORS; i++)
        {
            for (auto& p : periodicParticleList_[i])
            {
                if (p != nullptr)
                {
                    if (p_it == p->particle)
                    {
                        toBeDeleted.insert(p);
                        p = nullptr;
                    }
                }
            }
        }
        
        for (int i = 0; i < NUMBER_OF_PROCESSORS; i++)
        {
            for (auto& p : periodicGhostList_[i])
            {
                if (p != nullptr)
                {
                    if (p_it == p->particle)
                    {
                        toBeDeleted.insert(p);
                        p = nullptr;
                    }
                }
            }
        }
    }
    
    //Delete ID's
    for (auto id_it : toBeDeleted)
    {
        delete id_it;
    }
    
}

/*!
 * \details When a particle is deleted the ID is also removed, but to ensure everything happens in
 * the same order a nullptr is introduced instead. This function removes these nullptr's by
 * adding the last entry of the list on such empty spot and making the list shorter.
 * \param[in] list A list of MPiPeriodicParticleIDBase's that need to be cleaned from empty entries
 */
void PeriodicBoundaryHandler::cleanCommunicationList(std::vector<MpiPeriodicParticleIDBase*>& list)
{
    for (int i = 0; i < list.size(); i++)
    {
        if (list[i] == nullptr)
        {
            list[i] = list.back();
            list.pop_back();
            i--;
        }
    }
}

void PeriodicBoundaryHandler::cleanCommunicationLists()
{
    for (int i = 0; i < NUMBER_OF_PROCESSORS; i++)
    {
        cleanCommunicationList(periodicParticleList_[i]);
        cleanCommunicationList(periodicGhostList_[i]);
    }
}

/*!
 * \details On compile time it is not clear how many processors will be used when running the simulation
 * and therefore the basic periodioc lists can't be initialised /bi the constructor. This
 * function initialises the vectors to avoid segmentationfaults.
 */
//TODO optimise the target finding strategy such that only "NEW target/receives" are added
void PeriodicBoundaryHandler::initialise()
{
    newPeriodicParticleList_ = std::vector<PeriodicList>(NUMBER_OF_PROCESSORS, PeriodicList(0));
    periodicParticleList_ = std::vector<PeriodicList>(NUMBER_OF_PROCESSORS, PeriodicList(0));
    periodicGhostList_ = std::vector<PeriodicGhostList>(NUMBER_OF_PROCESSORS, PeriodicGhostList(0));
}

void PeriodicBoundaryHandler::performActionsBeforeAddingParticles()
{
    for (BasePeriodicBoundary* boundary : *this)
    {
        boundary->performActionsBeforeAddingParticles();
    }
}

/*!
 * \details This function deletes all periodic/ghost ID's and all corresponding ghosts.
 * starting with a fresh clean periodicBoundaryHandler. Useful when another boundary is
 * added - after - particles already have been added
 */
void PeriodicBoundaryHandler::clearCommunicationLists()
{
    if (NUMBER_OF_PROCESSORS > 1)
    {
        //Clear ID lists
        for (auto& i : periodicParticleList_)
        {
            //logger(INFO,"Size: %",periodicParticleList_[i].size());
            for (int j = 0; j < i.size(); j++)
            {
                delete i[j];
            }
            i.clear();
        }
        
        for (auto& i : periodicGhostList_)
        {
            for (int j = 0; j < i.size(); j++)
            {
                delete i[j];
            }
            i.clear();
        }
        
        //Collect all ghost particles and unflag
        ParticleHandler& pH = getDPMBase()->particleHandler;
        std::set<BaseParticle*> toBeDeleted; //Set because I am too lazy to implement a vector based flushParticle function
        for (BaseParticle* particle : pH)
        {
            if (particle->isPeriodicGhostParticle())
            {
                toBeDeleted.insert(particle);
            }
            
            //Unflag its periodicity
            particle->setInPeriodicDomain(false);
        }
        
        //Flush from mpi boundary and clean the lists
        getDPMBase()->getCurrentDomain()->flushParticles(toBeDeleted);
        getDPMBase()->getCurrentDomain()->cleanCommunicationLists();
        
        //Delete particles
        int index = 0;
        for (BaseParticle* particle : toBeDeleted)
        {
            pH.removeGhostObject(particle->getIndex());
            index++;
        }
    }
}

/*!
 * \details When a periodic particle changes complexity and remains real, 
 * there is a possibility the particle moved over a maser boundary. If that is indeed
 * the case then flag the particle as being not a maser particle
 * \param[in,out] particle The particle that needs to check the isMaserParticle flag
 * \bug This only works if all insertions happen on the root if there is a maser present, as the ID is increased locally
 */
void PeriodicBoundaryHandler::updateMaserParticle(BaseParticle* const particle)
{
    if (particle->isMaserParticle())
    {
        for (int b = 0; b < getSize(); b++)
        {
            if (particle->getPeriodicComplexity(b) == 3)
            {
                particle->setMaserParticle(false);
                
                logger(VERBOSE, "particle % with position % goes into outflow domain, new ID = %", particle->getId(),
                       particle->getPosition(), getDPMBase()->particleHandler.getNextId());
                const unsigned int newID = getDPMBase()->particleHandler.getNextId();
                logger(VERBOSE, "new id % position X %", newID, particle->getPosition().X);
                particle->setId(newID);
                getDPMBase()->particleHandler.increaseId();
            }
        }
    }
}

bool PeriodicBoundaryHandler::checkIfAddNewParticle(BaseParticle* particle)
{
    if (particle->isPeriodicGhostParticle())
    {
        return false;
    }
    if (particle->isInPeriodicDomain())
    {
        return false;
    }
    return true;
}











