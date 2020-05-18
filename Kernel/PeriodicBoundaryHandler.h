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


#ifndef PERIODICBOUNDARYHANDLER_H
#define PERIODICBOUNDARYHANDLER_H

#include "BaseHandler.h"
#include "Boundaries/BasePeriodicBoundary.h"
#include "MpiDataClass.h"
#include <set>


typedef std::vector<MpiPeriodicParticleID*> PeriodicList;
typedef std::vector<MpiPeriodicGhostParticleID*> PeriodicGhostList;

/*!
 * \class PeriodicBoundaryHandler
 * \brief Container to store pointers to all BasePeriodicBoundary objects
 * \details The PeriodicBoundaryHandler is a container to store all BasePeriodicBoundary. It is implemented by a vector
 * of pointers to BasePeriodicBoundary.
 */
class PeriodicBoundaryHandler final : public BaseHandler<BasePeriodicBoundary>
{
public:
    /*!
     * \brief Default constructor, it creates an empty PeriodicBoundaryHandler.
     */
    PeriodicBoundaryHandler();
    
    /*!
     * \brief Constructor that copies all BasePeriodicBoundary it contains and sets the other variables to 0/nullptr.
     */
    PeriodicBoundaryHandler(const PeriodicBoundaryHandler& BH);
    
    /*!
     * \brief Assignment operator, copies only the vector of BasePeriodicBoundary and sets the other variables to 0/nullptr.
     */
    PeriodicBoundaryHandler operator=(const PeriodicBoundaryHandler& rhs);
    
    /*!
     * \brief Destructor, it destructs the PeriodicBoundaryHandler and all BasePeriodicBoundary it contains.
     */
    ~PeriodicBoundaryHandler() override;
    
    /*!
     * \brief Adds a BasePeriodicBoundary to the PeriodicBoundaryHandler. 
     */
    void addObject(BasePeriodicBoundary* P) override;
    
    ///Pure virtual function needs implementation, but it does nothing for the periodicBoudnaryHandler
    void readAndAddObject(std::istream& is) override;
    
    /*!
     *  \brief Returns the name of the handler, namely the string "PeriodicBoundaryHandler".
     */
    std::string getName() const override;
    
    /*!
     * \brief Sets the interaction distance
     */
    void setInteractionDistance(Mdouble interactionDistance);
    
    /*!
     * \brief Returns the interaction distance
     */
    Mdouble getInteractionDistance();
    
    /*!
     * \brief Updates the positions/velocity of ghost particles and accordingly the status of these particles
     */
    void updateStatus(std::set<BaseParticle*>& ghostParticlesToBeDeleted);
    
    /*!
     * \brief Shifts the position of the particle based on its current periodic complexity
     */
    void shiftParticle(BaseParticle* particle);
    
    /*!
     * \brief Shifts the position of the particle based on a given periodic complexity
     */
    void shiftParticle(BaseParticle* particle, const std::vector<int>& complexity);
    
    /*!
     * \brief Computes the periodic complexity based on a given position
     */
    std::vector<int> computePeriodicComplexity(Vec3D position);
    
    /*!
     * \brief Computes the periodic complexity and total periodic complexity based on a given position
     */
    void computePeriodicComplexity(std::vector<int>& periodicComplexity, int& totalPeriodicComplexity, Vec3D position);
    
    /*!
     * \brief Adds new particles to the periodic particle lists
     */
    void addNewParticles();
    
    /*!
     * \brief Adds a new particle to the periodic list
     */
    void addNewParticle(BaseParticle* particle);
    
    /*!
     * \brief Returns the number of particles that are flagged is periodicGhostParticle
     */
    unsigned int getNumberOfPeriodicGhostParticles();
    
    /*!
     * \brief Returns the number of particles that are flagged as periodicGhostParticles,
     * but not as MPIParticles
     */
    Mdouble getNumberOfTruePeriodicGhostParticles();
    
    /*!
     * \brief Determines if a given particle is in the MPI domain and if it is an MPI Particle
     */
    void getMPIFlags(BaseParticle* particle, bool& isInMPIDomain, bool& isMPIParticle);
    
    /*!
     * \brief Sets the MPIParticle and isMPIParticle flags of a given particle
     */
    void setMPIFlags(BaseParticle* particle);
    
    /*!
     * \brief generates a list of periodic complexities corresponding to a give real particle.
     */
    void generateGhosts(std::vector<std::vector<int> >& list, std::vector<int> periodicComplexity,
                        std::vector<int>& complexity, int level);
    
    /*!
     * \brief Collects ghost particle data that needs to be be sent to other processors
     */
    void collectGhostParticleData();
    
    /*!
     * \brief Collects interaction data into an MPI data structure
     */
    void collectInteractionData();
    
    /*!
     * \brief Processes the received ghost data, creates a ghost particle and does some book keeping
     */
    void processReceivedGhostParticleData(int targetIndex, std::vector<BaseParticle*>& newParticles);
    
    /*!
     * \brief Process the received interaction data.
     */
    void processReceivedInteractionData(int targetIndex, std::vector<BaseParticle*>& newParticles);
    
    /*!
     * \brief Process the interaction data for local ghosts.
     */
    void processLocalInteractionData(std::vector<BaseParticle*>& newParticles);
    
    /*!
     * \brief Creates a periodioc particle ID for book keeping and moves the ID to the correct list
     */
    void processPeriodicParticles();
    
    /*!
     * \brief Creates ghost particles of periodic particles that are located on the same processor
     */
    void processLocalGhostParticles(std::vector<BaseParticle*>& newParticles);
    
    /*!
     * \brief Updates position/velocity and periodic complexity of ghost particles.
     */
    void updateParticles();
    
    /*!
     * \brief checks if a periodic complexity is real
     */
    bool checkIsReal(std::vector<int> complexity);
    
    /*!
     * \brief checks of two periodic complexities differ
     */
    bool checkChanged(std::vector<int> previousComplexity, std::vector<int> complexity);
    
    /*!
     * \brief Updates the status of periodic particles and ghost particles
     */
    void updateParticleStatus(std::set<BaseParticle*>& particlesToBeDeleted);
    
    /*!
     * \brief For a given complexity this function returns the target processor
     */
    int findTargetProcessor(const std::vector<int>& complexity);
    
    /*!
     * \brief Checks if a particle is in the periodic domain, but are not flagged as being in the periodic domain
     */
    void findNewParticle(BaseParticle* particle);
    
    bool checkIfAddNewParticle(BaseParticle* particle);
    
    /*!
     * \brief Loops over all particles in the simulation to check if they need to be added to the periodic lists
     */
    void findNewParticles();
    
    /*!
     * \brief Finds interactions that accompany future ghost particles
     */
    void findNewInteractions();
    
    /*!
     * \brief Creats a list of  send and receive targets for periodic/ghost particles
     */
    void communicateTargetDomains();
    
    /*!
     * \brief Communicate the number of new particles and interactions to target processors
     */
    void communicateNumberOfNewParticlesAndInteractions();
    
    /*!
     * \brief Initial preparation work for sending ghost particles
     */
    void prepareNewParticleTransmission();
    
    /*!
     * \brief Collects and sends the ghost particle data
     */
    void performNewParticleTransmission();
    
    /*!
     * \brief creates the ghost particles and performs some bookkeeping to keep track of them
     */
    void finaliseNewParticleTransmission();
    
    /*!
     * \brief Collects the position and velocity data from periodic boundaries
     */
    void preparePositionAndVelocityUpdate();
    
    /*!
     * \brief Communicates position and velocity data from periodic boundaries and updates ghost particles
     */
    void finalisePositionAndVelocityUpdate();
    
    /*!
     * \brief Removes particles from the periodiocParticleList_ and periociGhostList_
      */
    void flushParticles(std::set<BaseParticle*>& particlesToBeFlushed);
    
    /*!
     * \brief Removes the nullptr's from a communication list
     */
    void cleanCommunicationList(std::vector<MpiPeriodicParticleIDBase*>& list);
    
    void cleanCommunicationLists();
    
    /*!
     * \brief Removes all ghost particles and bookkeeping for a fresh start
     */
    void clearCommunicationLists();
    
    /*!
     * \brief Initialises the communication list vectors as they can not be determined on compile time
     */
    void initialise();
    
    /*!
    * \brief Flushes periodioc particles that need to be deleted from the periodic lists
    */
    void flushPeriodicParticles(std::set<BaseParticle*>& particlesToBeDeleted);
    
    /*!
     * \brief Actions that boundaries perform before adding new periodic/ghost particles
     */
    void performActionsBeforeAddingParticles();
    
    /*!
     * \brief Updates the maser flag of particles leaving the maser
     */
    void updateMaserParticle(BaseParticle* particle);
    
    /*!
     * \brief Disables boundaries that need to be ignored (i.e. a non-maser particle needs to ignore the maser boundary)
     */
    void
    findBoundariesToIgnore(BaseParticle* particle, std::vector<int>& periodicComplexity, int& totalPeriodicComplexity);


private:
    /*!
     * \brief The interaction distance between a position and the boundary for which
     * particles start to participate with a boundary or not
     */
    Mdouble interactionDistance_;
    
    /*!
     * \brief A list that keeps track which target processors the current processor
     * is receiving new particles from
     */
    std::vector<int> receiveTargetList_;
    
    /*!
     * \brief A list that keeps track to which targets this processor is sending
     * new particles to
     */
    std::vector<int> sendTargetList_;
    
    /*!
     * \brief A vector that stores how many new ghost particles will be received
     * from other processors.
     */
    std::vector<int> numberOfNewPeriodicGhostParticlesReceive_;
    
    /*!
     * \brief A vector that stores how many particles are going to be send to other
     * processors.
     */
    std::vector<int> numberOfNewPeriodicGhostParticlesSend_;
    
    /*!
     * \brief Stores the number of new interactions to be send to target processor corresponding to sendTargetList_
     */
    std::vector<int> numberOfNewInteractionsSend_;
    
    /*!
     * \brief Stores the number of new interactions to be received from target processor corresponding to receiveTargetList_
    */
    std::vector<int> numberOfNewInteractionsReceive_;
    
    /*!
     * \brief Data container for particles that are being received from other processors
     */
    std::vector<std::vector<MPIParticle> > periodicGhostParticleReceive_;
    
    /*!
     * \brief Data container for particles that are being send to other processors
     */
    std::vector<std::vector<MPIParticle> > periodicGhostParticleSend_;
    
    /*!
     * \brief Data container for periodic complexity that is being received from other processors
     */
    std::vector<std::vector<int> > periodicGhostComplexityReceive_;
    
    /*!
     * \brief Data container for periodic complexity that is being send to other processors
     */
    std::vector<std::vector<int> > periodicGhostComplexitySend_;
    
    /*!
     * \brief Data container for position data that is being received from other processors
     */
    std::vector<std::vector<MPIParticlePosition> > updatePositionDataReceive_;
    
    /*!
     * \brief Data container for position data that is being send to other processors
     */
    std::vector<std::vector<MPIParticlePosition> > updatePositionDataSend_;
    
    /*!
     * \brief Data container for velocity data that is being received from other processors
     */
    std::vector<std::vector<MPIParticleVelocity> > updateVelocityDataReceive_;
    
    /*!
     * \brief Data container for velocity data that is being send to other processors
     */
    std::vector<std::vector<MPIParticleVelocity> > updateVelocityDataSend_;
    
    /*!
     * \brief Stores the interaction data that is going to be send
     */
    std::vector<void*> interactionDataSend_;
    
    /*!
     * \brief Stores the interaction data that is going to be received
     */
    std::vector<void*> interactionDataReceive_;
    
    /*!
     * A vector the size of the number of processors, each entry containing
     * a vector of newly found periodic particles that need to be send to 
     * other processors
     */
    std::vector<PeriodicList> newPeriodicParticleList_;
    
    /*!
     * A list that stores the new interactions that have to be send to target processor, corresponding to sendTargetList_
     */
    std::vector<std::vector<BaseInteraction*> > newInteractionList_;
    
    /*!
     * \brief A vector the size of the number of processors, each entry containing
     *  a vector of periodic particle ID's to keep track of periodic particles and their 
     *  corresponding ghosts.
     */
    std::vector<PeriodicList> periodicParticleList_;
    
    /*!
     * \brief A vector the size of the number of processors, each entry containing
     * a vector of ghost periodioc particle ID's to keep track of periodic ghost particles
     * and their corresponding real particles
     */
    std::vector<PeriodicGhostList> periodicGhostList_;
};

#endif

