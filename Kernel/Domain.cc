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

#include "Domain.h"
#include "DPMBase.h"
#include "Logger.h"
#include "MpiContainer.h"
#include "MpiDataClass.h"
#include "Particles/SphericalParticle.h"
#include "Walls/BaseWall.h"
#include "InteractionHandler.h"
#include "DomainHandler.h"
#include "Math/Vector.h"
#include <limits>
#include <utility>
#include <vector>
#include <set>
#include <Particles/LiquidFilmParticle.h>

/*!
 * \brief Constructs an empty domain
 * \details Constructor of the domain class. It creates a domain with
 * infinite bounds and no processor or domainHandler assigned to it
 * and initialises all communication lists
 */
Domain::Domain()
{
    constructor();
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"Domain::Domain() finished"<<std::endl;
#endif
}

/*!
 * \brief Constructs an empty domain with a globalMeshIndex
 * \details Constructor of the domain class. It creates a domain with
 * infinite bounds and no processor or domainHandler assigned to it
 * and initialises all communication lists
 */
Domain::Domain(std::vector<unsigned> globalMeshIndex) : globalMeshIndex_(std::move(globalMeshIndex))
{
    constructor();
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"Domain::Domain() finished"<<std::endl;
#endif
}

/*!
 * \brief Copies a domain
 * \param[in] PH The domain that has to be copied. 
 * \details This constructor copies all the private variables of 
 * the current domain into the input domain.
 */
Domain::Domain(const Domain& b)
        : BaseObject(b)
{
    rank_ = b.rank_;
    domainHandler_ = b.domainHandler_;
    domainMin_ = b.domainMin_;
    domainMax_ = b.domainMax_;
    globalMeshIndex_ = b.globalMeshIndex_;
    middle_ = b.middle_;
    
    //A cube has 3^3=27 neighbours
    unsigned long numberOfNeighbours = 27;
    
    //Create all lists
    localIndexToGlobalIndexTable_ = std::vector<int>(numberOfNeighbours);
    localIndexToProcessorList_ = std::vector<int>(numberOfNeighbours);
    boundaryParticleList_ = std::vector<std::vector<BaseParticle*> >(numberOfNeighbours, std::vector<BaseParticle*>(0));
    boundaryParticleListNeighbour_ = std::vector<std::vector<BaseParticle*> >(numberOfNeighbours,
                                                                              std::vector<BaseParticle*>(0));
    newBoundaryParticleList_ = std::vector<std::vector<BaseParticle*> >(numberOfNeighbours,
                                                                        std::vector<BaseParticle*>(0));
    newInteractionList_ = std::vector<std::vector<BaseInteraction*> >(numberOfNeighbours,
                                                                      std::vector<BaseInteraction*>(0));
    numberOfParticlesSend_ = std::vector<unsigned>(numberOfNeighbours);
    numberOfParticlesReceive_ = std::vector<unsigned>(numberOfNeighbours);
    numNewInteractionsSend_ = std::vector<unsigned>(numberOfNeighbours);
    numNewInteractionsReceive_ = std::vector<unsigned>(numberOfNeighbours);
    boundaryParticleDataSend_ = std::vector<std::vector<MPIParticle> >(numberOfNeighbours);
    boundaryParticleDataReceive_ = std::vector<std::vector<MPIParticle> >(numberOfNeighbours);
    updatePositionDataSend_ = std::vector<std::vector<MPIParticlePosition> >(numberOfNeighbours);
    updatePositionDataReceive_ = std::vector<std::vector<MPIParticlePosition> >(numberOfNeighbours);
    updateVelocityDataSend_ = std::vector<std::vector<MPIParticleVelocity> >(numberOfNeighbours);
    updateVelocityDataReceive_ = std::vector<std::vector<MPIParticleVelocity> >(numberOfNeighbours);
    interactionDataSend_ = std::vector<void*>(numberOfNeighbours);
    interactionDataReceive_ = std::vector<void*>(numberOfNeighbours);
    activeBoundaryList_ = std::vector<bool>(numberOfNeighbours, true);
    boundaryList_ = std::vector<int>(0);
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"Domain::Domain(const Domain &b) finished"<<std::endl;
#endif
}

/*!
 * \details Destructor
 */
Domain::~Domain()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout << "Domain::~Domain() finished"<<std::endl;
#endif
}

/*!
 * \brief constructs the domain
 * \details the boundaries are set to infinity, no processor is assigned. all lists that require a full 3D domain element to function
 * are created
 */
void Domain::constructor()
{
    rank_ = constants::intMax;
    domainHandler_ = nullptr;
    domainMin_ = {-constants::inf, -constants::inf, -constants::inf};
    domainMax_ = {constants::inf, constants::inf, constants::inf};
    
    //A cube has 3^3=27 neighbours
    int numberOfNeighbours = 27;
    
    //Create all lists
    localIndexToGlobalIndexTable_ = std::vector<int>(numberOfNeighbours);
    localIndexToProcessorList_ = std::vector<int>(numberOfNeighbours);
    boundaryParticleList_ = std::vector<std::vector<BaseParticle*> >(numberOfNeighbours, std::vector<BaseParticle*>(0));
    boundaryParticleListNeighbour_ = std::vector<std::vector<BaseParticle*> >(numberOfNeighbours,
                                                                              std::vector<BaseParticle*>(0));
    newBoundaryParticleList_ = std::vector<std::vector<BaseParticle*> >(numberOfNeighbours,
                                                                        std::vector<BaseParticle*>(0));
    newInteractionList_ = std::vector<std::vector<BaseInteraction*> >(numberOfNeighbours,
                                                                      std::vector<BaseInteraction*>(0));
    numberOfParticlesSend_ = std::vector<unsigned>(numberOfNeighbours);
    numberOfParticlesReceive_ = std::vector<unsigned>(numberOfNeighbours);
    numNewInteractionsSend_ = std::vector<unsigned>(numberOfNeighbours);
    numNewInteractionsReceive_ = std::vector<unsigned>(numberOfNeighbours);
    boundaryParticleDataSend_ = std::vector<std::vector<MPIParticle> >(numberOfNeighbours);
    boundaryParticleDataReceive_ = std::vector<std::vector<MPIParticle> >(numberOfNeighbours);
    updatePositionDataSend_ = std::vector<std::vector<MPIParticlePosition> >(numberOfNeighbours);
    updatePositionDataReceive_ = std::vector<std::vector<MPIParticlePosition> >(numberOfNeighbours);
    updateVelocityDataSend_ = std::vector<std::vector<MPIParticleVelocity> >(numberOfNeighbours);
    updateVelocityDataReceive_ = std::vector<std::vector<MPIParticleVelocity> >(numberOfNeighbours);
    interactionDataSend_ = std::vector<void*>(numberOfNeighbours);
    interactionDataReceive_ = std::vector<void*>(numberOfNeighbours);
    activeBoundaryList_ = std::vector<bool>(numberOfNeighbours, true);
    boundaryList_ = std::vector<int>(0);
}

/*!
 * \brief Copy method.
 * \details Uses copy constructor to create a copy on the heap. 
 *          Useful for polymorphism.
 * \return pointer to the domain's copy
 */
Domain* Domain::copy() const
{
    return new Domain(*this);
}

/*!
 * \brief Reads the object from a given istream
 * \details Reads the object's id_ from the given istream
 * \param[in,out] is    istream the id_ is read from
 * \todo MX: why would this function be called?
 */
void Domain::read(std::istream& is)
{
    logger(WARN, "[Domain::read] should not be called");
    //BaseObject::read(is);
}

/*!
 * \details Adds the object's id_ to the given ostream
 * \param[in] os    ostream the id_ is added to
 * \todo MX: why would this function be called?
 */
void Domain::write(std::ostream& os) const
{
    logger(WARN, "[Domain::write] should not be called");
    //BaseObject::write(os);
}

/*!
 * \brief Returns the object's class name (i.e. 'DeletionBoundary').
 * \return the object's class name
 */
std::string Domain::getName() const
{
    return "Domain";
}

/*!
 * \details Sets the domain's range in a given direction. 
 * \param[in] direction  Direction 0,1,2
 * corresponds to direction x,y,z respectively
 * \param[in] min Minimum domain bound
 * \param[in] max Maximum domain bound
 */
void Domain::setRange(Direction direction, Mdouble min, Mdouble max)
{
    if (min > max)
    {
        logger(ERROR, "[MercuryMPI ERROR]: min is larger than max. (%,%)", min, max);
    }
    
    double maxClosed;
    double minClosed;
    if (min == -constants::inf)
    {
        minClosed = getHandler()->getDPMBase()->getMin().getComponent(direction);
    }
    
    if (max == constants::inf)
    {
        maxClosed = getHandler()->getDPMBase()->getMax().getComponent(direction);
    }
    
    
    switch (direction)
    {
        case Direction::XAXIS :
            domainMin_[0] = min;
            domainMax_[0] = max;
            middle_.X = minClosed + (maxClosed - minClosed) / 2.0;
            break;
        case Direction::YAXIS :
            domainMin_[1] = min;
            domainMax_[1] = max;
            middle_.Y = minClosed + (maxClosed - minClosed) / 2.0;
            break;
        case Direction::ZAXIS :
            domainMin_[2] = min;
            domainMax_[2] = max;
            middle_.Z = minClosed + (maxClosed - minClosed) / 2.0;
            break;
        default :
            logger(ERROR, "Direction is not a valid direction. (%)", direction);
            break;
    }
}

/*!
 * \details This function sets the bound for the domain, it also computes the middle
 * of the square domain, useful for mixed periodic/mpi particles.
 * \param[in] domainMin Minimum values of the domain in x-,y-,z-direction 
 * \param[in] domainMax Maximum values of the domain in x-,y-,z-direction
 */
void Domain::setBounds(std::vector<double> domainMin, std::vector<double> domainMax, bool computeMiddle)
{
    domainMin_ = domainMin;
    domainMax_ = domainMax;
    
    //Compute the middle of the closed domain
    if (computeMiddle)
    {
        double minClosed;
        double maxClosed;
        for (int i = 0; i < 3; i++)
        {
            minClosed = domainMin_[i];
            maxClosed = domainMax_[i];
            if (domainMin_[i] == -constants::inf)
            {
                minClosed = getHandler()->getDPMBase()->getMin().getComponent(i);
            }
            
            if (domainMax_[i] == constants::inf)
            {
                maxClosed = getHandler()->getDPMBase()->getMax().getComponent(i);
            }
            
            middle_.setComponent(i, minClosed + (maxClosed - minClosed) / 2.0);
        }
    }
}

/*!
 * \details Sets the domain's rank, corresponding to a processorID. 
 * \param[in] rank  Integer corresponding to a processorID
 */
void Domain::setRank(int rank)
{
    rank_ = rank;
}

/*!
 * \details Gets the domain's minimum domain bounds
 * \return Minimum bounds of the domain
 */
std::vector<double> Domain::getDomainMin()
{
    return domainMin_;
}

/*!
 * \details Gets the domain's maximum domain bounds
 * \return Maximum bounds of the domain
 */
std::vector<double> Domain::getDomainMax()
{
    return domainMax_;
}

/*!
 * \details Gets the domain's rank, corresponding to a processorID. 
 * \return rank of the domain
 */
int Domain::getRank()
{
    return rank_;
}

/*!
 * \details Sets the pointer to the DomainHandler the domain belongs to
 * \param[in] domainHandler   pointer to the DomainHandler
 */
void Domain::setHandler(DomainHandler* domainHandler)
{
    domainHandler_ = domainHandler;
}

/*!
 *  \details Returns the global index of the domain in the whole simulation domain
 * \returns Returns the global index of the domain
 */
int Domain::getGlobalIndex()
{
    return globalIndex_;
}

/*!
 * \details returns a vector with the globa mesh indices (i,j,k)
 * \return global mesh indices i,j and k
 */
std::vector<unsigned> Domain::getGlobalMeshIndex()
{
    return globalMeshIndex_;
}

/*!
 * \brief sets the global mesh index
 * \details global mesh index has the shape of (i,j,k)
 * \param[in] globalMeshIndex the mesh index vector i,j,k
 */
void Domain::setGlobalMeshIndex(std::vector<unsigned> globalMeshIndex)
{
    globalMeshIndex_ = globalMeshIndex;
}

/*!
 * \details This function disables the communication of a boundary. Each domain
 * in the parallel code has 27 boundaries: (6 sides, 12 ribs, 8 corners) and 
 * for convenience one is associated with the domain itself (used for local indexing). 
 * Boundaries that are inactive, such as domains that extend to infinity, do not need to be considered in the parallel computation.
 * This function disables that boundary.
 * \param[in] localIndex The local index of a neighbour of the current domain
 */
void Domain::disableBoundary(unsigned localIndex)
{
    activeBoundaryList_[localIndex] = false;
}

/*!
 * \details Returns a list of boundaries that are active in mpi communication
 *  \return List of active mpi communication boundaries
 */
std::vector<bool> Domain::getActiveBoundaryList()
{
    return activeBoundaryList_;
}

/*!
 * \brief Checks if a particle is in the domain
 * \details Checks if the particle is in the domain, but with an offset. 
 * this offset is set to zero standard
 * \param[in] particle Pointer to a particle
 * \param[in] offset offset from the domain boundary
 * \return True if the particle is in the domain. False if it is outside
 */
bool Domain::containsParticle(BaseParticle* particle, Mdouble offset)
{
    for (unsigned i = 0; i < 3; i++)
    {
        if (!(((domainMin_[i] + offset) < particle->getPosition().getComponent(i))
              &&
              ((domainMax_[i] - offset) >= particle->getPosition().getComponent(i)))
                )
        {
            return false;
        }
    }
    return true;
}

/*!
 * \details Check if the particle is inside the domain + the communication zone around it
 * \param[in] particle   Pointer to a particle
 * \return True if the particle is in the domain or the whole communication zone around it. False if it is outside
 */
bool Domain::isInGreaterDomain(BaseParticle* particle)
{
    return containsParticle(particle, -domainHandler_->getInteractionDistance());
}

/*! 
 * \brief Check if the particle is in the current domain but not in the communication zone
 * \param[in] particle Pointer to a particle
 * \return Returns true if the particle is in the inner domain (i.e. not in the communication zone of the domain)
 */
bool Domain::isInInnerDomain(BaseParticle* particle)
{
    return containsParticle(particle, domainHandler_->getInteractionDistance());
}

/*!
 * \details Check if the particle is inside the communication zone
 * \param[in] particle   Pointer to a particle
 * \return True if the particle is in the communication zone, i.e. in the domain, but not in the inner domain.
 * False if it is outside the communication zone.
 */
bool Domain::isInCommunicationZone(BaseParticle* particle)
{
    //First check if the particle is actually in the domain
    //Secondly check if the particle is in the inner domain
    return containsParticle(particle) && !(isInInnerDomain(particle));
}

/*!
 * \details Function that creates a lookup table from the local index system to the global index system
 */
void Domain::createLookUpTable()
{
    int localIndex;
    int globalIndex;
    //Get the global decomposition vector, (nx,ny,nz) and compute the mesh multipliers: (1,nx,nx*ny)
    std::vector<unsigned> numberOfDomains = domainHandler_->getNumberOfDomains();
    std::vector<unsigned> globalMeshMultiplier = {1,
                                                  numberOfDomains[Direction::XAXIS],
                                                  numberOfDomains[Direction::XAXIS] *
                                                  numberOfDomains[Direction::YAXIS]};
    
    //Compute the globalIndex of this domain
    globalIndex_ = globalMeshMultiplier[Direction::XAXIS] * globalMeshIndex_[Direction::XAXIS] +
                   globalMeshMultiplier[Direction::YAXIS] * globalMeshIndex_[Direction::YAXIS] +
                   globalMeshMultiplier[Direction::ZAXIS] * globalMeshIndex_[Direction::ZAXIS];
    
    //Create the lookup table from localIndex to globalIndex of the neihbours
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            for (int k = -1; k < 2; k++)
            {
                //Get the local index  
                localIndex = i + 3 * j + 9 * k + 13;
                //Compute the global index
                globalIndex = globalMeshMultiplier[Direction::XAXIS] * (globalMeshIndex_[Direction::XAXIS] + i) +
                              globalMeshMultiplier[Direction::YAXIS] * (globalMeshIndex_[Direction::YAXIS] + j) +
                              globalMeshMultiplier[Direction::ZAXIS] * (globalMeshIndex_[Direction::ZAXIS] + k);
                //Fill in the look-up table
                localIndexToGlobalIndexTable_[localIndex] = globalIndex;
            }
        }
    }
    //Create the lookup table from localIndex to processorRank
    rank_ = globalIndex_;
    localIndexToProcessorList_ = localIndexToGlobalIndexTable_;
}

/*! 
 * \details Given i,j,k (where i=0,j=0,k=0 is the current domain), 
 * get the local index of a neighbouring domain
 * \param[in] i Index in the x-direction with i=0 corresponding to the current domain
 * \param[in] j Index in the y-direction with j=0 corresponding to the current domain
 * \param[in] k Index in the z-direction with k=0 corresponding to the current domain
 * \return The localIndex of the neighbour is returned.
 */
int Domain::getLocalIndex(const int i, const int j, const int k)
{
    return i + 3 * j + 9 * k + 13;
}

/*! 
 * \details Given i,j,k (where i=0,j=0,k=0 is the current domain), 
 * get the local index of a neighbouring domain
 * \param[in] localMeshIndex (i,j,k) where (0,0,0) corresponds to the current domain
 * \return The localIndex of the neighbour is returned.
 */
int Domain::getLocalIndex(const std::vector<int> localMeshIndex)
{
    return localMeshIndex[Direction::XAXIS] + 3 * localMeshIndex[Direction::YAXIS] +
           9 * localMeshIndex[Direction::ZAXIS] + 13;
}

/*!
 * \brief Searches for a particle with a specific id in a list of particles
 * \details When interactions are copied over MPI from an MPI particle, the unique id of the particles they
 * interact with need to be looked up. This functions searches for the correct particle
 * \todo MX: This is probably not the fastest method
 * \param[in] indentification The unique ID of a particle
 * \param[in] particleList List of particles 
 * \return Either the particle with the id identification or a nullptr when the particle has not been found
 */
BaseParticle* Domain::findParticleInList(unsigned int identification, std::vector<BaseParticle*> particleList)
{
    for (BaseParticle* particle : particleList)
    {
        if (particle->getId() == identification)
        {
            return particle;
        }
    }
    return nullptr;
}

/*!
 * \brief This function finds if a given particle is close to a given boundary
 * \details This function detects in how many boundaries the particle is located and 
 * is stored in the boundaryIndex.
 * The boundary Index in the x-direction:
 * a) -1 if the particle is at the left boundary
 * b)  0 if the particle is not at a boundary
 * c)  1 if the particle is at the right boundary
 * The same holds for the y and z direction as well
 * \param[in] particle a base particle
 * \param[out] offset the L1 distance between the domain boundary and the wanted boundary
 * \return boundaryIndex a vector indicating at which boundary the particle is
 */
std::vector<int> Domain::findNearbyBoundaries(BaseParticle* particle, Mdouble offset)
{
    std::vector<int> boundaryIndex(3);
    Mdouble interactionDistance = domainHandler_->getInteractionDistance();
    
    //for x,y,z directions, find if the particle is close to the left or right wall or not at all
    for (int d = 0; d < 3; d++)
    {
        //Check if the particle is close to Lx
        if ((particle->getPosition().getComponent(d)) < (domainMin_[d] + interactionDistance + offset))
        {
            //Update switch
            boundaryIndex[d] = -1;
        }
            //Check if particle is close to Rx
        else if ((particle->getPosition().getComponent(d)) >= (domainMax_[d] - interactionDistance - offset))
        {
            //Update switch
            boundaryIndex[d] = 1;
        }
    }
    return boundaryIndex;
}

/*!
 * \brief disables all domain boundaries that have no neighbour
 * \details This function can only disable a square mesh. If you want to make a
 * concave mesh then this function has to be adapted, especially for the case when
 * the numberOfDomains[d] = 1.
 */
void Domain::disableBoundaries()
{
    std::vector<unsigned> numberOfDomains = domainHandler_->getNumberOfDomains();
    std::vector<int> localMeshIndex(3);
    int localIndex;
    
    //disble own list
    activeBoundaryList_[13] = false;
    
    //Special case when numberOfDomains[d] = 1 -> disable all d direction boundaries
    for (unsigned d = 0; d < 3; d++) //Loop over all dimensions
    {
        //If there is only one domain in the d-direction, take action
        if (numberOfDomains[d] == 1)
        {
            //Loop over all locaIndices
            for (int i = -1; i < 2; i++)
            {
                localMeshIndex[Direction::XAXIS] = i;
                for (int j = -1; j < 2; j++)
                {
                    localMeshIndex[Direction::YAXIS] = j;
                    for (int k = -1; k < 2; k++)
                    {
                        localMeshIndex[Direction::ZAXIS] = k;
                        //Disable all boundaries that are not in the "middle" of the mesh in the d-direction
                        // i.e. a localMeshIndex[d] = -1 is a neighbour in the d-direction, however there is only one
                        // element in that direction so it can't be a real neighbour.
                        if (localMeshIndex[d] != 0)
                        {
                            //Disable the boundary
                            localIndex = getLocalIndex(localMeshIndex);
                            activeBoundaryList_[localIndex] = false;
                        }
                    }
                }
            }
        }
    }
    
    //Compute the boundaryIndex for all other cases where numberOfDomains is larger than 1
    std::vector<int> boundaryIndex(3);
    for (int d = 0; d < 3; d++)
    {
        //Check if the domain is on the simulationDomainMin boundary
        if (globalMeshIndex_[d] == 0)
        {
            //Update boundary index
            boundaryIndex[d] = -1;
        }
            //Check if the domain is on the simulationDomainMax boundary
        else if (globalMeshIndex_[d] == (numberOfDomains[d] - 1))
        {
            //Update boundary index
            boundaryIndex[d] = 1;
        }
    }
    
    //Disable the boundaries that are on the edge of the simulation domain
    for (int d = 0; d < 3; d++)
    {
        //If boundary is on the edge of the simulation domain
        if (boundaryIndex[d] != 0)
        {
            //Loop over all local indices
            for (int i = -1; i < 2; i++)
            {
                for (int j = -1; j < 2; j++)
                {
                    for (int k = -1; k < 2; k++)
                    {
                        localMeshIndex = {i, j, k};
                        //Set the localMeshIndex to the boundaryIndex in the d-direction
                        //Since we are only interested in the domains on this boundary
                        localMeshIndex[d] = boundaryIndex[d];
                        //Disable the boundary
                        localIndex = getLocalIndex(localMeshIndex);
                        activeBoundaryList_[localIndex] = false;
                    }
                }
            }
        }
    }
    
    //Create a list of the active boundaries
    localIndex = 0;
    for (bool active : activeBoundaryList_)
    {
        if (active)
        {
            boundaryList_.push_back(localIndex);
        }
        localIndex++;
    }
}

/*!
 * \brief Function that adds the particles to the approriate boundary lists
 * \details This function computes the complexity of the particle. The complexity is defined by
 * the number of boundaries the particle is close-by. Depending on this complexity value the particle
 * has to be added to different boundary lists.
 * \param[in] particle pointer to base particle
 * \param[in,out] list a list of lists of particles which the particle might be added to
 */
void Domain::addParticlesToLists(BaseParticle* particle, std::vector<std::vector<BaseParticle*> >& list)
{
    std::vector<int> boundaryIndex = findNearbyBoundaries(particle);
   
    //Compute and set complexity of the particle
    unsigned int complexity = boundaryIndex[0] + 3 * boundaryIndex[1] + 9 * boundaryIndex[2] + 13;
    unsigned int list_complexity = 0;
    for (int d = 0; d < 3; d++) //Loop over all directions
    {
        list_complexity += std::abs(boundaryIndex[d]);
    }
    particle->setCommunicationComplexity(complexity);
    //particle->setCommunicationComplexity(list_complexity);
    
    //Based on the complexity of the particle, add it to the approriate list
    switch (list_complexity)
    {
        //The particle is not close at all
        case 0:
            break;
            //The particle is close to one side
            // 1 side contribution
        case 1 :
            //Add the side contribution
            list[getLocalIndex(boundaryIndex)].push_back(particle);
            break;
            //The particle is close to two neighbouring directions
            //2 side and 1 rib contrubution
        case 2 :
        {
            //Add the two side contributions
            for (int d = 0; d < 3; d++)
            {
                std::vector<int> localMeshIndex = {0, 0, 0};
                localMeshIndex[d] = boundaryIndex[d];
                //Avoid adding the particle in the wrong direction by excluding localMeshIndex[d] = 0
                if (localMeshIndex[d] != 0)
                {
                    list[getLocalIndex(localMeshIndex)].push_back(particle);
                }
            }
        }
            //Add the rib contribution
            list[getLocalIndex(boundaryIndex)].push_back(particle);
            break;
            
            //The particle is close to three neighbouring directions
            // 3 side, 3 rib and 1 corner contribution
        case 3 :
        {
            //Add the three side contributions
            for (int d = 0; d < 3; d++)
            {
                std::vector<int> localMeshIndex = {0, 0, 0};
                //Reset index vector
                localMeshIndex[d] = boundaryIndex[d];
                list[getLocalIndex(localMeshIndex)].push_back(particle);
            }
            
            //Add the three rib contributions
            for (int d = 0; d < 3; d++)
            {
                std::vector<int> localMeshIndex = boundaryIndex;
                //All rib boundary indices are given by the boundaryIndex and setting one of the components to zero
                localMeshIndex[d] = 0;
                list[getLocalIndex(localMeshIndex)].push_back(particle);
            }
        }
            
            //Add the corner contribution
            list[getLocalIndex(boundaryIndex)].push_back(particle);
            break;
        
        default :
            logger(INFO, "boundaryIndex : %,%,% | list_complexity: %", boundaryIndex[0], boundaryIndex[1], boundaryIndex[2],
                   list_complexity);
            logger(ERROR, "Particle is in contact with the wrong number of boundaries");
            break;
    }
}

/*!
 * \brief Function that finds new particles in the particle handler that should be added to the communication lists
 * \param[in] particleHandler the container that contains all particles
 */
void Domain::findNewMPIParticles(const ParticleHandler& particleHandler)
{
    //For all particles inside the given domain, loop over all the particles to
    //see if we have to add the particles to the new particle list
    for (BaseParticle* particle : particleHandler)
    {
        findNewMPIParticle(particle);
    }
}

/*!
 * \brief Function that checks if a given particle should be added to the communication lists and adds it accordingly
 * \param[in] particle point to a particle that is being checked if it should be added to the communication lists
 */
void Domain::findNewMPIParticle(BaseParticle* particle)
{
    //If the particle is an MPI particle
    if (!particle->isInMPIDomain() && !particle->isPeriodicGhostParticle())
    {
        addParticlesToLists(particle, newBoundaryParticleList_);
    }
}

bool Domain::isInNewBoundaryParticleList(BaseParticle* objectI,int localIndex) const
{
    for (BaseParticle *q : newBoundaryParticleList_[localIndex]) {
        if (q == objectI)
            return true;
    }
    return false;
}

/*!
 * \brief Finds interactions that have to be sent over to another domain
 * \details Newly added particles to the communication lists might have history interactions with
 * the particles already there. The other domains need to be aware of these history interactions and
 * this function collects the relevant interactions.
 */
void Domain::findNewMPIInteractions()
                            {
                                //For all active boundaries
                                for (int localIndex : boundaryList_)
                                {
                                    //Check all newly added particles for interactions with already known particles
                                    for (BaseParticle* newBoundaryParticle : newBoundaryParticleList_[localIndex])
                                    {
                                        //Loop over all interactions of the new particle
                                        std::vector<BaseInteraction*> interactions = newBoundaryParticle->getInteractions();
                                        for (BaseInteraction* interaction : interactions)
                                        {
                                            //Find out what the new particle is interacting with
                                            BaseParticle* particleP = dynamic_cast<BaseParticle*>(interaction->getP());
                                            BaseParticle* objectI = dynamic_cast<BaseParticle*>(interaction->getI());

                                            //If the P in the interaction structure is the new particle, find I
                                            if (newBoundaryParticle == particleP)
                                            {
                                                //Check if the new particle is interacting with a wall
                                                if (!objectI)
                                                {
                                                    //Check if the interaction is still valid
                                                    BaseWall* wall = dynamic_cast<BaseWall*>(interaction->getI());
                                                    Mdouble distance;
                                                    Vec3D normal;
                                                    wall->getDistanceAndNormal(*particleP, distance, normal);
                                                    if (distance <= particleP->getMaxInteractionRadius())
                                                    {
                                                        //Add the interaction to the list if there are still in contact (hence the if statement)
                                                        newInteractionList_[localIndex].push_back(interaction);
                        }
                    }
                    else //is I a particle
                    {
                        //check if particle is in mpi domain OR (TODO) if it is in the newBoundaryParticleList_
                        if (objectI->isInMPIDomain() || isInNewBoundaryParticleList(objectI, localIndex) ) {
                            //Is the particle still interacting after the position update?
                            Vec3D branchVector = particleP->getPosition() - objectI->getPosition();
                            //Get the square of the distance between particle i and particle j
                            Mdouble distanceSquared = Vec3D::getLengthSquared(branchVector);
                            Mdouble sumOfInteractionRadii =
                                    objectI->getSumOfInteractionRadii(particleP);
                            if (distanceSquared < (sumOfInteractionRadii * sumOfInteractionRadii))
                            {
                                //Add the interaction to the list
                                newInteractionList_[localIndex].push_back(interaction);
                            }
                        }
                    }
                }
                else //newBoundaryParticle is I in the interaction, P can only be a particle
                {
                    //it is "I" in the interaction structure, check if its in the mpi domain
                    if (particleP->isInMPIDomain() || isInNewBoundaryParticleList(particleP, localIndex))
                    {
                        //Is the particle still interacting after the position update?
                        Vec3D branchVector = particleP->getPosition() - objectI->getPosition();
                        //Get the square of the distance between particle i and particle j
                        Mdouble distanceSquared = Vec3D::getLengthSquared(branchVector);
                        Mdouble sumOfInteractionRadii =
                                objectI->getSumOfInteractionRadii(particleP);
                        if (distanceSquared < (sumOfInteractionRadii * sumOfInteractionRadii))
                        {
                            //Add the interaction to the list
                            newInteractionList_[localIndex].push_back(interaction);
                        }
                    }
                }
            }
        }
    }
}

/*!
 * \brief collects the data of a particles that has to be communicated to other processors
 * \details Only the relevant data of a baseParticle is being send to other processors, therefore
 * this function collects the relevant data and puts it in a format compatible with data transmission
 * \param[in] localIndex the local index of a boundary
 */
void Domain::collectBoundaryParticleData(int localIndex)
{
    //For all particles
    for (BaseParticle* particle : newBoundaryParticleList_[localIndex])
    {
        //Add the data to the transmission list
        boundaryParticleDataSend_[localIndex].push_back(copyDataFromParticleToMPIParticle(particle));
    }
}

/*!
 * \brief Collects the data of an interaction that has to be communicated to other processors
 * \details Interactions have to be reformatted in an MPI-suitable format to transmit over to other processors
 * \param[in] localIndex The local index of a boundary
 */
void Domain::collectInteractionData(int localIndex)
{
    //Copy interactions to the data array
    unsigned int indexInteraction = 0;
    for (BaseInteraction* interaction : newInteractionList_[localIndex])
    {
        interaction->getMPIInteraction(interactionDataSend_[localIndex], indexInteraction);
        indexInteraction++;
    }
}

/*!
 * \brief Function that copies the mpi data format of a base particle to a real particle and adds it to the particleHandler
 * \details When adding a particle with this function it's always a ghost particle and therefor gets the MPIParticle flag
 * \param[in] index The boundary index of the boundary list currently being treated
 * \param[in,out] particleHandler Container containing the particles of the current domain, will contain the new MPI particles as well
 */
void Domain::processReceivedBoundaryParticleData(const unsigned index, std::vector<BaseParticle*>& newParticles)
{
    
    ParticleHandler& particleHandler = getHandler()->getDPMBase()->particleHandler;
    //for all MPIPparticles that have been received on the other side
    for (int i = 0; i < numberOfParticlesReceive_[index]; i++)
    {
        //copy data from data
        BaseParticle* p0 = MPIParticle::newParticle();
        copyDataFromMPIParticleToParticle(&(boundaryParticleDataReceive_[index][i]), p0, &particleHandler);
        
        //Flag that the particle is a ghost particle of a real particle in the neighbour domain
        p0->setMPIParticle(true);
        //Flag that the particle is list as a particle in the mpi domain
        p0->setInMPIDomain(true);
        //add particle to handler
        particleHandler.addGhostObject(p0);
        //Set previous position as current position
        BaseParticle* pGhost = particleHandler.getLastObject();
        pGhost->setPreviousPosition(particleHandler.getLastObject()->getPosition());
        //Add to the newParticle list
        newParticles.push_back(pGhost);
        logger(VERBOSE, "Adding mpi particle % at: %", particleHandler.getLastObject()->getId(),
               particleHandler.getLastObject()->getPosition());
        
        //add pointer to the particle to the list and set the status to idle
        boundaryParticleListNeighbour_[index].push_back(particleHandler.getLastObject());
    }
}

/*!
 * \brief Bookkeep the newly send particles
 * \details Each domain keeps a list of particles that have copies on other domains. The order of these
 * lists is exactly the same order of the neighbour particle lists of the neighbour domain.
 * \param[in] index The boundary index of the boundary list currently being treated
 */
void Domain::processSentBoundaryParticles(const unsigned index)
{
    for (BaseParticle* particle : newBoundaryParticleList_[index])
    {
        boundaryParticleList_[index].push_back(particle);
        particle->setInMPIDomain(true);
    }
}

/*!
 * \brief Processes the received interactions from newly added mpi particles
 * \details The received interaction data is scanned for interaction id's P and I
 * they are then matches to the actual particle or wall and the interaction is created.
 * The history details of the interaction are filled in after creation.
 * \bug There is a subtle and very rare bug in here, where two particles with an interaction across a periodic boundary
 * pass into the next domain at the same time step. In this case, the interactions should be copied over before the
 * ghost-particles are constructed. However, since the interaction is with the ghost particle, this cannot be done. That
 * is why there is a 'continue' if one of the particles of the interaction is not found. An illustration of this
 * situation is given in PeriodicBoundaryEnteringMpiDomainMPI2Test.
 * \param[in] localIndex The boundary index of the boundary list currently being treated
 * \param[in,out] interactionHandler Handler containing all interactions on the current domain
 */
void Domain::processReceivedInteractionData(const unsigned localIndex, std::vector<BaseParticle*>& newParticles)
{
    InteractionHandler& iH = getHandler()->getDPMBase()->interactionHandler;
    for (unsigned int l = 0; l < numNewInteractionsReceive_[localIndex]; l++)
    {
        unsigned int identificationP;
        unsigned int identificationI;
        bool isWallInteraction;
        unsigned timeStamp;
        
        //Get the general information required to setup a new interaction
        iH.getInteractionDetails(interactionDataReceive_[localIndex],
                                 l, identificationP, identificationI,
                                 isWallInteraction, timeStamp);
        
        logger(VERBOSE, "interaction details: %, %, idP %, idI %, wall %, time %", interactionDataReceive_[localIndex],
               l, identificationP, identificationI,
               isWallInteraction, timeStamp);
        
        //Find the particle in the newParticle list
        BaseParticle* pGhost = nullptr;
        int idOther;
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
        if (pGhost == nullptr)
        {
            logger(WARN, "In Domain::processReceivedInteractionData: pGhost (id %) is nullptr, the interaction data is not copied. Two particles possibly moved into domain simultaneously.", identificationP);
            continue;
        }
        
        //If it is a wall interaction, do stuff
        if (isWallInteraction)
        {
            BaseInteractable* I = getHandler()->getDPMBase()->wallHandler.getObjectById(identificationI);
            //Create interactions
            BaseInteraction* j = I->getInteractionWith(pGhost, timeStamp, &iH);
            if (j!= nullptr) j->setMPIInteraction(interactionDataReceive_[localIndex], l, false);
            
        }
        else
        {
            //Obtain potential interaction particles
            std::vector<BaseParticle*> interactingParticleList;
            getHandler()->getDPMBase()->hGridGetInteractingParticleList(pGhost, interactingParticleList);
            
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
            if (otherParticle == nullptr) {
                //search for the other particle in the newParticles list
                for (BaseParticle *particle : newParticles) {
                    if (particle->getId() == idOther) {
                        otherParticle = particle;
                        break;
                    }
                }
                if (otherParticle == nullptr) {
                    logger(WARN,
                           "In Domain::processReceivedInteractionData: otherParticle (id %) is nullptr, the interaction data with pGhost (id %) is not copied. Two particles possibly moved into domain simultaneously. nt = %",
                           identificationI, identificationP, getHandler()->getDPMBase()->getNumberOfTimeSteps());
                    continue;
                }
            }
            //Add the interaction
            BaseInteraction* j = pGhost->getInteractionWith(otherParticle, timeStamp, &iH);
            if (j!= nullptr) j->setMPIInteraction(interactionDataReceive_[localIndex], l, false);
        }
    }
}

void Domain::debugInformation() {
    std::stringstream s;
    std::stringstream m;
    for (auto p : getHandler()->getDPMBase()->particleHandler) {
        if (p->isMPIParticle())
            m << std::setw(4) << p->getId();
        else
            s << std::setw(4) << p->getId();
    }
    logger(INFO,"Particles %, MPI %",s.str(),m.str());
}


/*!
 * \brief A symmetric communication between two domains exchanging a send/recieve count
 * \details Before real data (i.e. particle data) can be transmitted to other processors, these receiving processors
 * need to know how many data objects they are receiving. This function facilitates a symmetric send mechanism between
* two domains to communicate how much objects they get from eachother
        * \param[in] tag A MercuryMPITag that indicates what type of count is being send
        * \param[out] countReceive unsigned integer containing the receiving count
        * \param[in] countSend unsigned integer containing the sending count
        * \param[in] localIndexNeighbour the local index of the domain with which the communication is being done
        */
        void Domain::sendAndReceiveCount(MercuryMPITag tag, unsigned& countReceive, unsigned& countSend,
                unsigned localIndexNeighbour)
        {
            int globalIndexNeighbour = localIndexToGlobalIndexTable_[localIndexNeighbour];
            int processor = localIndexToProcessorList_[localIndexNeighbour];

            //Create communication tags
            int tagReceive = globalIndexNeighbour * MAX_PROC + globalIndex_ * 10 + tag;
            int tagSend = globalIndex_ * MAX_PROC + globalIndexNeighbour * 10 + tag;

            logger.assert(tagSend > 0, "Send tag is wrong. Tag: %", tagSend);
            logger.assert(tagReceive > 0, "Receive tag is wrong. Tag: %", tagReceive);
            logger.assert(processor >= 0, "Processor is wrong. processor: %", processor);

            //Communicate the requests
            MPIContainer::Instance().receive(countReceive, processor, tagReceive);
            MPIContainer::Instance().send(countSend, processor, tagSend);
        }

/*!
 * \brief Prepares the MPI transmission of particle and interaction data from particles in particleHandler
 * \details This function checks all the particles in the particle handler to see if they have moved
 * in the communication domain. If that is the case they are added to a list of particles that need to be send.
 * After all particles are evaluated the number of new particles to be send is send to the appropriate processors
 * The same is being done to the interactions of the newly found particles
 * \param[in] particleHandler handles all the particles in the current domain
 */
        void Domain::prepareBoundaryDataTransmission()
        {
            //Find new particles that have entered the communication zone
            findNewMPIParticles(getHandler()->getDPMBase()->particleHandler);

            //Find interactions between new particles and other particles in the communication zone
            findNewMPIInteractions();

            //Compute number of particles to send
            //For all boundaries
            for (int localIndex : boundaryList_)
            {
                numberOfParticlesSend_[localIndex] = newBoundaryParticleList_[localIndex].size();
                numNewInteractionsSend_[localIndex] = newInteractionList_[localIndex].size();
            }

            //For all active boundaries
            for (int localIndex : boundaryList_)
            {
                //Send and receive the number of new boundary particles
                sendAndReceiveCount(MercuryMPITag::PARTICLE_COUNT, (numberOfParticlesReceive_[localIndex]),
                                    (numberOfParticlesSend_[localIndex]), localIndex);

                //Send and receive the new interactions
                sendAndReceiveCount(MercuryMPITag::INTERACTION_COUNT, (numNewInteractionsReceive_[localIndex]),
                                    (numNewInteractionsSend_[localIndex]), localIndex);
            }
        }

/*!
 * \brief Prepares the MPI transmission of a single particle and its interactions
 * \details When a single particle has to be added to the communication zone (i.e. when the user adds a particle during
 * the simulation) then all processors should be informed if they receive the particle or not. The functionality of this function
 * is practically the same as that of prepareBoundaryDataTransmission(ParticleHandler &particleHandler)
 * The same is being done to the interactions of the newly found particles
 * \param[in] particle A base particle that needs to be transmitted
 */
        void Domain::prepareBoundaryDataTransmission(BaseParticle* particle)
        {
            //Assign the particle to the correct lists if the particle belongs to this domain
            if (containsParticle(particle))
            {
                findNewMPIParticle(particle);
            }

            //Note: When inserting a single particle, it has no interactions, so no interactions to be copied either

            //Compute number of particles to send
            //For all boundaries
            for (int localIndex : boundaryList_)
            {
                numberOfParticlesSend_[localIndex] = newBoundaryParticleList_[localIndex].size();
                numNewInteractionsSend_[localIndex] = newInteractionList_[localIndex].size();
            }


            //For all active boundaries
            for (int localIndex : boundaryList_)
            {
                //Send and recieve the number of new boundary particles
                sendAndReceiveCount(MercuryMPITag::PARTICLE_COUNT, (numberOfParticlesReceive_[localIndex]),
                                    (numberOfParticlesSend_[localIndex]), localIndex);

                //Send and receive the new interactions
                sendAndReceiveCount(MercuryMPITag::INTERACTION_COUNT, (numNewInteractionsReceive_[localIndex]),
                                    (numNewInteractionsSend_[localIndex]), localIndex);
            }
        }

/*!
 * \brief Collects data to be transmitted and then performs the transmission of the data
 * \details The particles and interactions are copied into a data class that can be transmitted to other processors.
 * The data is then accordingly send to the appropriate processors
 * \param[in] interactionHandler Handler that contains all interactions, required to collect interaction data
 */
        void Domain::performBoundaryDataTransmission()
        {
            //For all active boundaries
            for (int localIndex : boundaryList_)
            {
                //make sure enough memory is reserved to receive the data
                boundaryParticleDataReceive_[localIndex].resize(numberOfParticlesReceive_[localIndex]);

                //Collect the particle data
                collectBoundaryParticleData(localIndex);

                //Send the data
                sendAndReceiveMPIData(MercuryMPITag::PARTICLE_DATA, MercuryMPIType::PARTICLE,
                                      boundaryParticleDataReceive_[localIndex].data(), numberOfParticlesReceive_[localIndex],
                                      boundaryParticleDataSend_[localIndex].data(), numberOfParticlesSend_[localIndex],
                                      localIndex);


                //create arrays for sending and receiving interaction data
                interactionDataSend_[localIndex] = getHandler()->getDPMBase()->interactionHandler.createMPIInteractionDataArray(
                        numNewInteractionsSend_[localIndex]);
                interactionDataReceive_[localIndex] = getHandler()->getDPMBase()->interactionHandler.createMPIInteractionDataArray(
                        numNewInteractionsReceive_[localIndex]);

                //Collect the interaction data
                collectInteractionData(localIndex);

                //Send the data
                sendAndReceiveMPIData(MercuryMPITag::INTERACTION_DATA, MercuryMPIType::INTERACTION,
                                      interactionDataReceive_[localIndex], numNewInteractionsReceive_[localIndex],
                                      interactionDataSend_[localIndex], numNewInteractionsSend_[localIndex], localIndex);
            }
        }

/*!
 * \brief This function processes the transmitted data
 * \details copies the received mpi data into the particles/interacions and adds them to the particleHandler/interactionHandler
 * All vectors that have been used in transmitting data are then cleaned up
 * \param[in,out] particleHandler Handles all particles in the current domain, on output contains the newly received particles
 * \param[in,out] interactionHandler Handles all interactions in the current domain, on output contains the newly received interactions
 */
        void Domain::finaliseBoundaryDataTransmission()
        {
            //For all boundaries
            for (int localIndex : boundaryList_)
            {
                //copy the received data into the particleHandler and neighbour particle list
                std::vector<BaseParticle*> newParticles;
                processReceivedBoundaryParticleData(localIndex, newParticles);
                //All particles in the current domain that have been send to the other domains need to be flagged as communicating particles
                processSentBoundaryParticles(localIndex);
                //copy the received interaction data into the standard dpm structure
                processReceivedInteractionData(localIndex, newParticles);

                //Delete all send/receive data
                boundaryParticleDataSend_[localIndex].clear();
                boundaryParticleDataReceive_[localIndex].clear();
                getHandler()->getDPMBase()->interactionHandler.deleteMPIInteractionDataArray(interactionDataSend_[localIndex]);
                getHandler()->getDPMBase()->interactionHandler.deleteMPIInteractionDataArray(
                        interactionDataReceive_[localIndex]);

                //Reset all counters
                numberOfParticlesSend_[localIndex] = 0;
                numberOfParticlesReceive_[localIndex] = 0;
                numNewInteractionsSend_[localIndex] = 0;
                numNewInteractionsReceive_[localIndex] = 0;

                //Reset all lists
                newBoundaryParticleList_[localIndex].clear();
                newInteractionList_[localIndex].clear();
            }
        }

/*!
 * \brief This step updates all communication lists and particles in the communication zone
 * \details After all the mpi particles have received a position and velocity update, a check is required to see if
 * the particle needs to be treated differently. Comments inside the code show all possible options. Note that
 * when a particle needs to be deleted, it will not be deleted straight away, but added to a list to be deleted later.
 * Additionally the location in the communication lists will not be emptied by assigning a nullptr to it. When
 * all particles are considered this nullptrs are removed from the list, to ensure that the mirror list on
 * another processor keeps the same order in particles.
 * A) The particle moves from the Neighbour domain to the current domain
 * B) The particle moves out of the communication zone, into the neighbouring domain
 * C) The particle moves out of the communication zone, into the current domain
 * D) The particle changes complexity: A special case where the particle required a re-evaluation in which boundary list it belongs
 * \param[in] particleHandler Handler that takes care of all particles in this domain
 */
        void Domain::updateParticles(std::set<BaseParticle*>& ghostParticlesToBeDeleted)
        {
            int complexityNew;
            std::vector<int> boundaryIndex;

            //For all active boundaries
            for (int localIndex : boundaryList_)
            {
                //Step 1A: Remove the particles from the boundaryParticleList_ which require re-assignment
                //or have just moved away from the region
                for (int p = 0; p < boundaryParticleList_[localIndex].size(); p++)
                {
                    BaseParticle* particle = boundaryParticleList_[localIndex][p];
                    //Check if the particle is still in the domain, but not in the communication zone
                    if (containsParticle(particle))
                    {
                        //check if the complexity has changed
                        boundaryIndex = findNearbyBoundaries(particle);
                        complexityNew = boundaryIndex[0] + 3 * boundaryIndex[1] + 9 * boundaryIndex[2] + 13;
                        if (particle->getCommunicationComplexity() != complexityNew)
                        {
                            logger(VERBOSE, "time: % | global index: % in list % | particle % | CURRENT DOMAIN - CHANGES "
                                            "COMPLEXITY", domainHandler_->getDPMBase()->getTime(), globalIndex_,
                                   localIndexToGlobalIndexTable_[localIndex], particle->getId());
                            //Flag the particle that it no longer participates in the communication layer
                            //so it can be reintroduced in the transmission step
                            particle->setMPIParticle(false);
                            particle->setInMPIDomain(false);
                            particle->setCommunicationComplexity(0);
                            boundaryParticleList_[localIndex][p] = nullptr;
                        }
                    }
                    else
                    {
                        logger(VERBOSE, "time: % | global index: % in list % | particle % | CURRENT DOMAIN - TO NEIGHBOURING "
                                        "DOMAIN", domainHandler_->getDPMBase()->getTime(), globalIndex_,
                               localIndexToGlobalIndexTable_[localIndex], particle->getId());

                        ghostParticlesToBeDeleted.insert(particle);
                        boundaryParticleList_[localIndex][p] = nullptr;
                    }
                }

                //Step 1B: Remove the particles from the boundaryParticleListNeightbour_ which require re-assignment
                //or have just moved away from the region
                for (int p = 0; p < boundaryParticleListNeighbour_[localIndex].size(); p++)
                {
                    BaseParticle* particle = boundaryParticleListNeighbour_[localIndex][p];
                    //Check if the particle moved out of the neighbour domain
                    if (!(domainHandler_->getObject(localIndexToGlobalIndexTable_[localIndex])->containsParticle(particle)))
                    {
                        //The particle has moved to this domain
                        if (containsParticle(particle))
                        {
                            logger(VERBOSE,
                                   "time: % | global index: % in list % | particle % | NEIGHBOURING DOMAIN - TO CURRENT DOMAIN",
                                   domainHandler_->getDPMBase()->getTime(), globalIndex_,
                                   localIndexToGlobalIndexTable_[localIndex], particle->getId());

                            //Flag the particle as not yet communicated
                            particle->setMPIParticle(false);
                            particle->setInMPIDomain(false);
                            particle->setCommunicationComplexity(0);
                            boundaryParticleListNeighbour_[localIndex][p] = nullptr;
                        }
                            //The particle has moved to a different domain
                        else
                        {
                            logger(VERBOSE,
                                   "time: % | global index: % in list % | particle % | NEIGBOURING DOMAIN - TO OTHER DOMAIN",
                                   domainHandler_->getDPMBase()->getTime(), globalIndex_,
                                   localIndexToGlobalIndexTable_[localIndex], particle->getId());
                            //Cruelly destroy the particle without any mercy.
                            ghostParticlesToBeDeleted.insert(particle);
                            boundaryParticleListNeighbour_[localIndex][p] = nullptr;
                        }
                    }
                    else
                    {
                        //check if the complexity has changed
                        boundaryIndex = domainHandler_->getObject(
                                localIndexToGlobalIndexTable_[localIndex])->findNearbyBoundaries(particle);
                        complexityNew = boundaryIndex[0] + 3 * boundaryIndex[1] + 9 * boundaryIndex[2] + 13;
                        if (particle->getCommunicationComplexity() != complexityNew)
                        {
                            logger(VERBOSE,
                                   "time: % | global index: % in list % | particle % | NEIGHBOURING DOMAIN - CHANGES COMPLEXITY",
                                   domainHandler_->getDPMBase()->getTime(), globalIndex_,
                                   localIndexToGlobalIndexTable_[localIndex], particle->getId());
                            //Cruelly destroy the particle without any mercy.
                            ghostParticlesToBeDeleted.insert(particle);
                            boundaryParticleListNeighbour_[localIndex][p] = nullptr;
                        }
                    }
                }
            }
        }

/*!
 * \brief Updates the position of particles which are flagged as MPIParticles
 * \details Manually updates the position of the particles, the displacement and the orientation
 * It additionally updates the position of the particle in the hGrid
 * \param[in] localIndex an index to the boundary being updated
 */
        void Domain::updateParticlePosition(int localIndex)
        {
            //process the updated information
            unsigned int index = 0;
            for (BaseParticle* particle : boundaryParticleListNeighbour_[localIndex])
            {
                logger.assert(particle->getId() == updatePositionDataReceive_[localIndex][index].id,
                              "MPI particle lists are not in sync");

                //set position
        particle->setPreviousPosition(particle->getPosition());
        particle->setPosition(updatePositionDataReceive_[localIndex][index].position);
        particle->setOrientation(updatePositionDataReceive_[localIndex][index].orientation);
        if (std::is_base_of<MPILiquidFilmParticle,MPIParticle>())
            static_cast<LiquidFilmParticle*>(particle)->setLiquidVolume(updatePositionDataReceive_[localIndex][index].liquidVolume);

        //Update hGrid
        Vec3D displacement = particle->getPreviousPosition() - particle->getPosition();
        getHandler()->getDPMBase()->hGridUpdateMove(particle, displacement.getLengthSquared());
        
        index++;
    }
}

/*!
 * \brief Updates the velocity of particles which are flagged as MPIParticles
 * \details Updates the translation but also the angular velocity
 * \param[in] localIndex an index to the boundary being updated
 */
void Domain::updateParticleVelocity(int localIndex)
{
    //process the updated information
    unsigned int index = 0;
    for (BaseParticle* particle: boundaryParticleListNeighbour_[localIndex])
    {
        //set velocity
        particle->setVelocity(updateVelocityDataReceive_[localIndex][index].velocity);
        particle->setAngularVelocity(updateVelocityDataReceive_[localIndex][index].angularVelocity);
        index++;
    }
}

/*!
 * \brief Function that sends particle position and velocity data for ghost particles
 * to other processors
 * \details When the real particle moves on it's actual domain, it's ghost particles have to get an update
 * this function sends the data to the processors that require the update
 */
void Domain::preparePositionAndVelocityUpdate()
{
    //For all active boundaries
    for (int localIndex : boundaryList_)
    {
        numberOfParticlesSend_[localIndex] = boundaryParticleList_[localIndex].size();
        numberOfParticlesReceive_[localIndex] = boundaryParticleListNeighbour_[localIndex].size();
        
        //Increase capacity for the receiving data files
        updatePositionDataReceive_[localIndex].resize(numberOfParticlesReceive_[localIndex]);
        updateVelocityDataReceive_[localIndex].resize(numberOfParticlesReceive_[localIndex]);
        
        
        //Collect data
        for (BaseParticle* particle : boundaryParticleList_[localIndex])
        {
            updatePositionDataSend_[localIndex].push_back(copyPositionFrom(particle));
            updateVelocityDataSend_[localIndex].push_back(copyVelocityFrom(particle));
        }
        
        //Send and receive the data
        sendAndReceiveMPIData(MercuryMPITag::POSITION_DATA, MercuryMPIType::POSITION,
                              updatePositionDataReceive_[localIndex].data(), numberOfParticlesReceive_[localIndex],
                              updatePositionDataSend_[localIndex].data(), numberOfParticlesSend_[localIndex],
                              localIndex);
        sendAndReceiveMPIData(MercuryMPITag::VELOCITY_DATA, MercuryMPIType::VELOCITY,
                              updateVelocityDataReceive_[localIndex].data(), numberOfParticlesReceive_[localIndex],
                              updateVelocityDataSend_[localIndex].data(), numberOfParticlesSend_[localIndex],
                              localIndex);
    }
}

/*!
 * \details After the data has been received from other processors, this function will update
 * the ghost particles with the new position and velocity data. Afterwards the status of the 
 * particles will be updated based on their new positions. Sometimes a particle needs to be removed
 * from the simulation, however since this particle might still be active in the periodic boundary 
 * it is not deleted straight away but stored in ghostParticlesToBeDeleted. 
 * Finally the communication data is removed.
 * \param[in,out] ghostParticlesToBeDeleted A vector containing particles that need to be removed from the simulation
 */
void Domain::finalisePositionAndVelocityUpdate(std::set<BaseParticle*>& ghostParticlesToBeDeleted)
{
    
    //For all active boundaries
    for (int localIndex : boundaryList_)
    {
        updateParticlePosition(localIndex);
        updateParticleVelocity(localIndex);
    }
    
    //Based on the new position, update the particle lists. 
    //Remove particles that left the communication zone, they will be re-communicated in a later step
    updateParticles(ghostParticlesToBeDeleted);
    
    //For all active boundaries clear the data lists
    for (int localIndex : boundaryList_)
    {
        updatePositionDataSend_[localIndex].clear();
        updateVelocityDataSend_[localIndex].clear();
        updatePositionDataReceive_[localIndex].clear();
        updateVelocityDataReceive_[localIndex].clear();
    }
    
}

/*!
 * \details This function is not in use
 */
void Domain::prepareVelocityUpdate()
{
    //For all active boundaries
    for (int localIndex : boundaryList_)
    {
        numberOfParticlesSend_[localIndex] = boundaryParticleList_[localIndex].size();
        numberOfParticlesReceive_[localIndex] = boundaryParticleListNeighbour_[localIndex].size();
        
        //Resize the vector to the correct size
        updateVelocityDataReceive_[localIndex].resize(numberOfParticlesReceive_[localIndex]);
        
        //Collect data
        for (BaseParticle* particle : boundaryParticleList_[localIndex])
        {
            updateVelocityDataSend_[localIndex].push_back(copyVelocityFrom(particle));
        }
        
        //Send the data
        sendAndReceiveMPIData(MercuryMPITag::VELOCITY_DATA, MercuryMPIType::VELOCITY,
                              updateVelocityDataReceive_[localIndex].data(), numberOfParticlesReceive_[localIndex],
                              updateVelocityDataSend_[localIndex].data(), numberOfParticlesSend_[localIndex],
                              localIndex);
    }
}

/*!
 * \details This function is not in use
 */
void Domain::finaliseVelocityUpdate()
{
    //For all active boundaries
    for (int localIndex : boundaryList_)
    {
        updateParticleVelocity(localIndex);
    }
    
    //For all active boundaries clear the data lists
    for (int localIndex : boundaryList_)
    {
        updateVelocityDataSend_[localIndex].clear();
        updateVelocityDataReceive_[localIndex].clear();
    }
}


/*!
 * \details Returns the pointer to the DomainHandler the domain belongs to
 * \return pointer to the handler
 */
DomainHandler* Domain::getHandler() const
{
    return domainHandler_;
}

/*!
 * \details Adds new MPI ghost particles to the simulation. This is done in a three stage process
 * where first the particles are located that generate ghosts, secondly their data is transmitted
 * to the correct processor and thirdly that data is then used to create ghost particles.
 */
void Domain::addNewParticles()
{
    //Step 1: For every MPIDomain boundary, create a list of particles that have to be transmitted
    //queue send and receive instructions for the number of particles
    prepareBoundaryDataTransmission();
    MPIContainer::Instance().sync();
    
    //Step 2: queue send and receive of data
    performBoundaryDataTransmission();
    MPIContainer::Instance().sync();
    
    //Step 3: Add the received particles to the particleHandler of the current domain
    finaliseBoundaryDataTransmission();
}

/*!
 * \brief Initialises a single particle which was added by the user after the domain creation
 * \details This function is used when a single particle is added during the simulation
 * In that case it has to be added to the mpi domains manually using this function.
 * Examples are insertion boundaries that occasionally add particles to the domain
 * \param[in] particle The single particle that has to be inserted into the domain
 */
void Domain::addParticle(BaseParticle* particle)
{
    //Step1: check if the particle has to be sent to other processors
    prepareBoundaryDataTransmission(particle);
    MPIContainer::Instance().sync();
    
    //Step2: queue send and receive data. Note for an inserted particle, no interactions should be required
    performBoundaryDataTransmission();
    MPIContainer::Instance().sync();
    
    //Step3: Add the received particles to the particleHandler of the current domain
    finaliseBoundaryDataTransmission();
    MPIContainer::Instance().sync();
}

/*!
 * \brief Updates particles that are not in the current domain and communicates newly added particles
 * \details First the newly positions of the ghost particles are communicated, based on the new positions these
 * particles will be updated accordingly. Secondly a new sweep through all particles is performed to see if we need
 * to add any particles to the list which have moved into the communication zone. Particles that need to be
 * removed from the simulation are stored in a vector and will be destroyed at a later point in the code.
 * \param[in,out] ghostParticlesToBeDeleted A list of particles which will be deleted afterwards
 */
void Domain::updateStatus(std::set<BaseParticle*>& ghostParticlesToBeDeleted)
{
    //Collect new positions and velocities and send them to the other domains
    preparePositionAndVelocityUpdate();
    MPIContainer::Instance().sync();
    
    //Receive the new positions and velocities from other domains
    //and update the mpi flagged particles accordingly. removes and switched particles in the lists
    finalisePositionAndVelocityUpdate(ghostParticlesToBeDeleted);
    MPIContainer::Instance().sync();
}

/*!
 *\details: not in use
 */
void Domain::updateVelocity()
{
    //collect new velocity data and send
    prepareVelocityUpdate();
    MPIContainer::Instance().sync();
    
    //process the received data
    finaliseVelocityUpdate();
    MPIContainer::Instance().sync();
}

/*!
 * \details All the MPIParticles are located in the neighbour lists and we only need to take a sum of these lists
 * to obtain the number of MPIParticles. This function is required to keep track of the number of real particles in the
 * domain.
 */
unsigned int Domain::getNumberOfMPIParticles()
{
    unsigned int count = 0;
    for (auto& index : boundaryParticleListNeighbour_)
    {
        count += index.size();
    }
    return count;
}

/*!
 * \details All the MPIParticles are located in the neighbour lists and we only need to take a sum of these lists
 * to obtain the number of MPIParticles. This function is required to keep track of the number of real particles in the
 * domain.
 */
unsigned int Domain::getNumberOfTrueMPIParticles()
{
    unsigned int count = 0;
    for (auto& index : boundaryParticleListNeighbour_)
    {
        for (auto& p : index)
        {
            if (!p->isPeriodicGhostParticle())
            {
                count++;
            }
        }
    }
    return count;
}

/*!
 * \param[in] toBeDeletedList A list of particles that is going to be deleted
 */
void Domain::flushParticles(std::set<BaseParticle*>& toBeFlushedList)
{
    //For all active boundaries
    for (int localIndex : boundaryList_)
    {
        flushParticlesFromList(boundaryParticleList_[localIndex], toBeFlushedList);
        flushParticlesFromList(boundaryParticleListNeighbour_[localIndex], toBeFlushedList);
    }
}

void Domain::flushParticlesFromList(std::vector<BaseParticle*>& list, std::set<BaseParticle*>& toBeFlushedList)
{
    //Firstly: turn all particles that need to be flushed into nullptrs
    for (auto& p : list)
    {
        if (p != nullptr)
        {
            BaseParticle* particle1 = p;
            for (BaseParticle* particle2 : toBeFlushedList)
            {
                //If the particle was found in the list, make a nullptr
                if (particle1 == particle2)
                {
                    logger(VERBOSE, "Removing particle from mpi domain at: %", particle1->getPosition());
                    p = nullptr;
                }
            }
        }
    }
}

/*!
 * \details Returns the middle of this square domain, required for periodic particles
 * \return middle_ The middle position of a cubic domain
 */
Vec3D Domain::getMiddle() const
{
    return middle_;
}

/*!
 * \details After particles have been updated, the communication lists contains 
 * nullptrs, remove these from the boundaryParticleList and boundaryParticleListNeighbour
 */
void Domain::cleanCommunicationLists()
{
    //For all active boundaries
    for (int i : boundaryList_)
    {
        cleanCommunicationList(boundaryParticleList_[i]);
        cleanCommunicationList(boundaryParticleListNeighbour_[i]);
    }
}

/*!
 * \details Removes nullptrs from a list of base particles efficiently
 * by replacing the last entry to an empty nullptr space.m
 * \param[in] list The list that needs to be cleansed from nullptrs
 */
void Domain::cleanCommunicationList(std::vector<BaseParticle*>& list)
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
