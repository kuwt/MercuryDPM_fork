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


/*!
 * \file Domain.h
 *
 * Class Domain
 */
#ifndef Domain_H
#define Domain_H

#include "BaseObject.h"
#include "MpiContainer.h"
#include "Math/Vector.h"
#include "ParticleHandler.h"
#include <vector>
#include <set>

class DPMBase;

class DomainHandler;

class BaseParticle;

class MPIParticle;

class MPIParticlePosition;

class MPIParticleVelocity;

class MPIParticleForce;

/*!
 * \class Domain
 * \brief The simulation can be subdivided into Domain's used in parallel code 
 * \details The domain class defines a region in the simulation domain, given by
 * a min vector and a max vector. The domain also as a processorID assigned to it.
 * It has a check to see whether a given particle is in the domain or not.
 */
class Domain final : public BaseObject
{
public:
    /*!
     * \brief Default Domain constructor
     */
    Domain();
    
    /*!
     * \brief
     */
    explicit Domain(std::vector<unsigned> globalMeshIndex);
    
    /*!
     * \brief Constructor that copies the domain range and rank from a 
     * given domain.
     */
    Domain(const Domain& d);
    
    /*!
     * \brief Destructor, destroys the domain.
     */
    ~Domain() override;
    
    /*!
     * \brief contructor of a domain
     */
    void constructor();
    
    /*!
     * \brief Function that creates a copy of this current domain, 
     * using the copy constructor.
     */
    virtual Domain* copy() const;
    
    /*!
     * \brief This function does nothing
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief This function does nothing
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*!
     * \brief Sets the domain range in a given direction
     */
    void setRange(Direction direction, Mdouble min, Mdouble max);
    
    /*!
     * \brief Sets the domain bounds
     */
    void setBounds(std::vector<double> domainLeft, std::vector<double> domainRight, bool computeMiddle);
    
    /*!
     * \brief Gets the minimum domain bounds
     */
    std::vector<double> getDomainMin();
    
    /*!
     * \brief Gets the maximum domain bounds
     */
    std::vector<double> getDomainMax();
    
    /*!
     * \brief Gets the rank associated with the assigned processorID
     */
    int getRank();
    
    /*!
     * \brief Sets the rank associated with the assigned processorID
     */
    void setRank(int rank);
    
    /*!
     * \brief Sets the domainHandler
     */
    void setHandler(DomainHandler* handler);
    
    /*!
     * \brief Gets the domainHandler
     */
    DomainHandler* getHandler() const;
    
    /*!
     * \brief Gets the global index of the domain
     */
    int getGlobalIndex();
    
    /*!
     * \brief Gets the global mesh index of the domain
     */
    std::vector<unsigned> getGlobalMeshIndex();
    
    /*! 
     * \brief Sets the global mesh index of theh domain
     */
    void setGlobalMeshIndex(std::vector<unsigned> globalMeshIndex);
    
    /*!
     * \brief Disables a boundary of the domain with a neighbouring domain
     */
    void disableBoundary(unsigned localIndex);
    
    /*!
     * \brief disables all domain boundaries that have no neighbour
     */
    void disableBoundaries();
    
    /*!
     * \brief Returns a list of boundaries that are active in mpi communication
     */
    std::vector<bool> getActiveBoundaryList();
    
    /*!
     * \brief Check to see if a given particle is within the current domain
     */
    bool containsParticle(BaseParticle* particle, Mdouble offset = 0.0);
    
    /*!
     * \brief Check to see if a given particle is in the current domain or in neighbouring
     * communication zones
     */
    bool isInGreaterDomain(BaseParticle* particle);
    
    /*! 
     * \brief Check if the particle is in the current domain but not in the communication zone
     */
    bool isInInnerDomain(BaseParticle* particle);
    
    /*!
     * \brief Check if the particle is in the communication zone of the current domain
     */
    bool isInCommunicationZone(BaseParticle* particle);
    
    /*!
     * \brief Create a look up table between local index system to global index system
     */
    void createLookUpTable();
    
    /*!
     * \brief return the local index of a domain given local mesh indices i,j and k
     */
    int getLocalIndex(int i, int j, int k);
    
    /*!
     * \brief return the local index of a doman given the localMeshIndex vector
     */
    int getLocalIndex(std::vector<int> localMeshIndex);
    
    /*! 
     * \brief Searches for a particle with a specific id in a list of particles
     */
    BaseParticle* findParticleInList(unsigned int identification, std::vector<BaseParticle*> particleList);
    
    /*!
     * \brief This function finds if a given particle is close to a given boundary
     */
    std::vector<int> findNearbyBoundaries(BaseParticle* particle, Mdouble offset = 0);
    
    /*!
     * \brief Function that adds the particles to the approriate boundary list
     */
    void addParticlesToLists(BaseParticle* particle, std::vector<std::vector<BaseParticle*> >& list);
    
    /*!
     * \brief Function that finds new particles in the particle handler that should be added to the communication lists
     */
    void findNewMPIParticles(const ParticleHandler& particleHandler);
    
    /*!
     * \brief Function that check if a given particle should be added to the communication lists
     */
    void findNewMPIParticle(BaseParticle* particle);

    /*
     * Checks whether a particle is in the list of new particles received by this thread
     */
    bool isInNewBoundaryParticleList(BaseParticle* object,int localIndex) const;
    /*!
     * \brief Finds interactions that have to be send over to another domain
     */
    void findNewMPIInteractions();
    
    /*!
     * \brief collects the data of a particle that has to be communicated to other processors
     */
    void collectBoundaryParticleData(int localIndex);
    
    /*!
     * \brief Collects the data of an interaction that has to be communicated to other processors
     */
    void collectInteractionData(int localIndex);
    
    /*!
     * \brief Function that copies the mpi data format of a base particle to a real particle and adds it to the particleHandler
     */
    void processReceivedBoundaryParticleData(unsigned index, std::vector<BaseParticle*>& newParticles);
    
    /*!
     * \brief Bookkeep the newly send particles
     */
    void processSentBoundaryParticles(unsigned index);
    
    /*!
     * \brief Processes the received interactions from newly added mpi particles    
     */
    void processReceivedInteractionData(unsigned index, std::vector<BaseParticle*>& newParticles);

    /*
     * Writes a list of all MPI particles on each thread to the logger.
     */
    void debugInformation();

    /*!
     * \brief A symmetric communication between two domains exchanging a send/recieve count
     */
    void
    sendAndReceiveCount(MercuryMPITag tag, unsigned& countReceive, unsigned& countSend, unsigned localIndexNeighbour);
    
    /*!
     * \brief Prepares the MPI transmission of particle and interaction data from particles in particleHandler
     */
    void prepareBoundaryDataTransmission();
    
    /*!
     * \brief Prepares the MPI transmission of a single particle and its interactions
     */
    void prepareBoundaryDataTransmission(BaseParticle* particle);
    
    /*!
     * \brief Collects data to be transmitted and then performs the transmission of the data
     */
    void performBoundaryDataTransmission();
    
    /*!
     * \brief This function processes the transmitted data
     */
    void finaliseBoundaryDataTransmission();
    
    /*!
     * \brief This step updates all communication lists and particles in the communication zone
     */
    void updateParticles(std::set<BaseParticle*>& ghostParticlesToBeDeleted);
    
    /*!
     * \brief Updates the position of particles which are flagged as MPIParticles
     */
    void updateParticlePosition(int localIndex);
    
    /*!
     * \brief Updates the velocity of particles which are flagged as MPIParticles
     */
    void updateParticleVelocity(int localIndex);
    
    /*!
     * \brief Function that sends particle position and velocity data for ghost particles
     * to other processors
     */
    void preparePositionAndVelocityUpdate();
    
    /*! 
     * \brief processes position and velocity data for ghost particles
     */
    void finalisePositionAndVelocityUpdate(std::set<BaseParticle*>& ghostParticlesToBeDeleted);
    
    /*!
     * \brief Function that sends particle velocity data for ghost particles
     */
    void prepareVelocityUpdate();
    
    /*!
     * \brief Processes particle velocity data for ghost particles
     */
    void finaliseVelocityUpdate();
    
    /*!
     * \brief Function that sends transmissionData/positionData/velocityData to other processors
     * \details This function is heavily templated on the type of data that can be send. All the types that
     * are located in the MercuryMPIType enum class can be used to transmit data. Examples are MPIParticle 
     * and MPIParticlePosition data types
     * \param[in] tag A MercuryMPITag that indicates what type of data is being send
     * \param[in] type A MercuryMPIType that tells the MPI what MPI_Data structure to use
     * \param[out] receiveData The data that will be received by the current domain
     * \param[in] receiveCount The number of items that are being received by the current domain
     * \param[in] sendData The data that will be send by the current domain
     * \param[in] sendCount The number of items that are being send by the current domain
     * \param[in] localIndexNeighbour the local index to the neighbouring domain
     */
    template<typename T>
    void sendAndReceiveMPIData(MercuryMPITag tag, MercuryMPIType type,
                               T* receiveData, unsigned receiveCount,
                               T* sendData, unsigned sendCount, unsigned localIndexNeighbour)
    {
        int globalIndexNeighbour = localIndexToGlobalIndexTable_[localIndexNeighbour];
        int processor = localIndexToProcessorList_[localIndexNeighbour];
        int tagReceive = globalIndexNeighbour * MAX_PROC + globalIndex_ * 10 + tag;
        int tagSend = globalIndex_ * MAX_PROC + globalIndexNeighbour * 10 + tag;
        
        //Communicate the requests
        if (receiveCount != 0)
        {
            MPIContainer::Instance().receive(receiveData, type, receiveCount, processor, tagReceive);
        }
        if (sendCount != 0)
        {
            MPIContainer::Instance().send(sendData, type, sendCount, processor, tagSend);
        }
    }
    
    /*!
     * \brief Initialises the MPIParticles by communicating newly found particles
     */
    void addNewParticles();
    
    /*!
     * \brief Initialises a single particle which is added during the simulation
     */
    void addParticle(BaseParticle* particle);
    
    /*!
     * \brief Updates particles that are not in the current domain and communicates newly added particles
     */
    void updateStatus(std::set<BaseParticle*>& ghostParticlesToBeDeleted);
    
    /*!
     * \brief Updates MPI particle velocity at the half-time step
     */
    void updateVelocity();
    
    /*!
     * \brief Obtains the number of particles in the particleHandler that are MPIParticles
     */
    unsigned int getNumberOfMPIParticles();
    
    /*!
     * \brief Obtains the number of particles in the particleHandler that are MPIParticles, but NOT periodic particles
     */
    unsigned int getNumberOfTrueMPIParticles();
    
    /*!
     * \brief Particles that are going to be deleted from the simulation are flushed out of the communication boundaries
     */
    void flushParticles(std::set<BaseParticle*>& toBeDeletedList);
    
    /*!
     * \brief Particles that are going to be deleted from the simulation are flushed out of a give communcation boundary
     */
    void flushParticlesFromList(std::vector<BaseParticle*>& list, std::set<BaseParticle*>& toBeDeletedList);
    
    /*!
     * \brief Gives the middle of the domain
     */
    Vec3D getMiddle() const ;

    /*!
     * \brief Removes nullptrs from boundaryParticleList_ and boundaryParticleListNeighbour_
     */
    void cleanCommunicationLists();
    
    /*!
     * \brief Removes nullptr's from a given particle list
    */
    void cleanCommunicationList(std::vector<BaseParticle*>& list);

private:
    /*!
     * \brief Pointer to the domain's DomainHandler container
     */
    DomainHandler* domainHandler_;
    
    /*!
     * \brief Minimum domain bounds in the x,y and z direction
     */
    std::vector<double> domainMin_;
    
    /*!
     * \brief Maximum domain bounds in the x,y and z direction
     */
    std::vector<double> domainMax_;
    
    /*!
     * \brief Middle of the closed domain
     */
    Vec3D middle_;
    
    /*!
     * \brief Global index of the domain in the mesh
     * \note for a standard decomposition this compares straight to the rank of the processor
     */
    int globalIndex_;
    
    /*!
     * \brief Vector containing the global mesh indices i,j,k
     */
    std::vector<unsigned> globalMeshIndex_;
    
    /*!
     * \brief look-up table to get the global index given a local domain index
     * \todo TW@Marnix should this be unsigned int
     */
    std::vector<int> localIndexToGlobalIndexTable_;
    
    /*!
     * \brief look-up table to get the processor of the domain given a local domain index
     * \todo TW@Marnix should this be unsigned int
     */
    std::vector<int> localIndexToProcessorList_;
    
    /*!
     * \brief A list of flags corresponding to an inactive or active boundary
     */
    std::vector<bool> activeBoundaryList_;
    
    /*!
     * \brief A list of indices of all the active boundaries
     * \todo TW@Marnix should this be unsigned int
     */
    std::vector<int> boundaryList_;
    
    /*!
     * \brief A list of boundary particles in the communication zone that are ghost particles on other domains
     */
    std::vector<std::vector<BaseParticle*> > boundaryParticleList_;
    
    /*!
     * \brief a list of ghost particles on the current domain, which are real on the neighbour domain
     */
    std::vector<std::vector<BaseParticle*> > boundaryParticleListNeighbour_;
    
    /*!
     * \brief Array that queues particles that need to be transmitted
     */
    std::vector<std::vector<BaseParticle*> > newBoundaryParticleList_;
    
    /*!
     * \brief Array that queues interactions that need to be transmitted
     */
    std::vector<std::vector<BaseInteraction*> > newInteractionList_;
    
    /*!
     * \brief Counter that keeps track of the number of particles that are being send to other domains
     */
    std::vector<unsigned> numberOfParticlesSend_;
    
    /*!
     * \brief Counter that keeps track of the number of particles that are being received by this domain
     */
    std::vector<unsigned> numberOfParticlesReceive_;
    
    /*!
     * \brief Counter that keeps track of the number of interactions that are being send to other domains
     */
    std::vector<unsigned> numNewInteractionsSend_;
    
    /*!
     * \brief Counter that keeps track of the number of interactions that are being received by this domain
     */
    std::vector<unsigned> numNewInteractionsReceive_;
    
    /*!
     * \brief Container that keeps a list of MPIParticles that are being send to other domains
     */
    std::vector<std::vector<MPIParticle> > boundaryParticleDataSend_;
    
    /*!
     * \brief Container that keeps a list of MPIParticles that are being received by this domain
     */
    std::vector<std::vector<MPIParticle> > boundaryParticleDataReceive_;
    
    /*!
     * \brief Container that keeps a list of MPIParticlePositions that are being send to other domains
     */
    std::vector<std::vector<MPIParticlePosition> > updatePositionDataSend_;
    
    /*!
     * \brief Container that keeps a list of MPIParticlePositions that are being received by this domain
     */
    std::vector<std::vector<MPIParticlePosition> > updatePositionDataReceive_;
    
    /*!
     * \brief Container that keeps a list of MPIParticleVelocities that are being send to other domains
     */
    std::vector<std::vector<MPIParticleVelocity> > updateVelocityDataSend_;
    
    /*!
     * \brief Container that keeps a list of MPIParticleVelocities that are being received by this domain
     */
    std::vector<std::vector<MPIParticleVelocity> > updateVelocityDataReceive_;
    
    /*!
     * \brief Container that keeps a void array of all the interaction data that are being send to other domains, interpretation is done by the interaction handler
     */
    std::vector<void*> interactionDataSend_;
    
    /*!
     * \brief Container that keeps a void array of all the interaction data that is being received by this domain, interpretation is done by the interaction handler
     */
    std::vector<void*> interactionDataReceive_;
    
    /*!
     * \brief Rank of the domain which identifies to which processor it belongs
     */
    int rank_;
};

#endif
