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
#include <string>
#include <cstdlib>
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include "SpeciesHandler.h"
#include "Species/BaseSpecies.h"
#include "Interactions/BaseInteraction.h"
#include "DPMBase.h"
#include "MpiDataClass.h"
#ifdef MERCURY_USE_OMP
#include "omp.h"
#endif

/*!
 * Constructor of the ParticleHandler class. It creates and empty ParticleHandler.
 */
InteractionHandler::InteractionHandler()
{
    writeVTK_ = FileType::NO_FILE;
    logger(DEBUG, "InteractionHandler::InteractionHandler() finished");
}

/*!
 * \param[in] IH The InteractionHandler that has to be copied.
 * \details This is not a copy constructor! It only clears all variables, since 
 *          by default interactions are not copied.
 * \todo Please check if interactions indeed don't need to be copied.
 */
InteractionHandler::InteractionHandler(const InteractionHandler& IH UNUSED)
        : BaseHandler<BaseInteraction>()
{
    writeVTK_ = IH.writeVTK_;
    //By default interactions are not copied.
    logger(DEBUG, "InteractionHandler::InteractionHandler(const "
                  "InteractionHandler &IH) finished, please note that no interactions"
                  " have been copied.");
}

/*!
 * \param[in] rhs The BoundaryHandler on the right hand side of the assignment.
 */
InteractionHandler& InteractionHandler::operator=(const InteractionHandler& rhs)
{
    if (this != &rhs)
    {
        clear();
    }
    writeVTK_ = rhs.writeVTK_;
    logger(DEBUG, "InteractionHandler::operator =(const InteractionHandler& rhs) finished.");
    return *this;
}

InteractionHandler::~InteractionHandler()
{
    logger(DEBUG, "InteractionHandler::~InteractionHandler() finished");
}

/*!
 * \param[in] P A pointer to the BaseInteraction (or derived class) that has to be added.
 */
void InteractionHandler::addObject(BaseInteraction* I)
{
    //set the particleHandler pointer
    I->setHandler(this);
    //Puts the interaction in the Interaction list
    BaseHandler<BaseInteraction>::addObject(I);
    I->getP()->addInteraction(I);
    I->getI()->addInteraction(I);
}

/*!
 * \param[in] P the first BaseInteractable by which the interaction is defined.
 * \param[in] I the first BaseInteractable by which the interaction is defined.
 * \return the Interaction between the BaseInteractable's P and I, if it exists.
 */
BaseInteraction*
InteractionHandler::getExistingInteraction(const BaseInteractable* const P, const BaseInteractable* const I) const
{
    //for particle-particle collision it is assumed BaseInteractable P has a lower index then I, so we only have to check for I, not P
    {
        for (unsigned j = 0; j < P->getInteractions().size(); ++j) {
            auto i = P->getInteractions()[j];
            if (i->getI() == I) {
                return i;
            }
        }
    }
    return nullptr;
}

BaseInteraction* InteractionHandler::addInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
{
    BaseSpecies* species = getDPMBase()->speciesHandler.getMixedObject(P->getIndSpecies(), I->getIndSpecies());
    BaseInteraction* C = species->getNewInteraction(P, I, timeStamp);
    addObject(C);
    return C;
}

void InteractionHandler::resetNewObjectsOMP() {
    #ifdef MERCURY_USE_OMP
    newObjects_.clear();
    newObjects_.resize(getDPMBase()->getNumberOfOMPThreads());
    #endif
}

void InteractionHandler::addNewObjectsOMP() {
    for (auto objects : newObjects_) {
        for (auto object : objects) {
            addObject(object);
        }
    }
}

/*!
 * \details Returns a pointer to the existing Interaction, if the Interaction 
 * already exists otherwise creates a new Interaction and returns a pointer to it.
 * \param[in] P the first BaseInteractable by which the interaction is defined.
 * \param[in] I the first BaseInteractable by which the interaction is defined.
 * \param[in] timeStamp the current value of DPMBase::time_.
 * \return the Interaction between the BaseInteractable's P and I
 */
BaseInteraction*
InteractionHandler::getInteraction(BaseInteractable* const P, BaseInteractable* const I, const unsigned timeStamp)
{
    const BaseSpecies* const species = getDPMBase()->speciesHandler.getMixedObject(P->getIndSpecies(),
                                                                                   I->getIndSpecies());
    
    BaseInteraction* C = getExistingInteraction(P, I);
    if (C == nullptr) {
        C = species->getNewInteraction(P, I, timeStamp);
        #ifdef MERCURY_USE_OMP
        if (omp_get_num_threads()>1) {
            newObjects_[omp_get_thread_num()].push_back(C);
            C->setHandler(this);
        } else {
            addObject(C);
        }
        #else
        addObject(C);
        #endif
    }
    
    //set timeStamp
    C->setTimeStamp(timeStamp);
    
    ///\todo TW this can bet set earlier
    //set species of collision
    C->setSpecies(species);
    
    return C;
}

/*!
 * \brief Creates an empty interaction
 * \details Empty interactions are required when dealing with parallel code.
 * it is used a proto type interaction that can handle the MPI interaction data type
 * which is a void array
 */
BaseInteraction* InteractionHandler::createEmptyInteraction() const
{
    //NOTE: assuming the interaction type of a species is the same for all species
    BaseSpecies* species = this->getDPMBase()->speciesHandler.getObject(0);
    BaseInteraction* emptyInteraction = species->getEmptyInteraction();
    
    return emptyInteraction;
}

/*!
 * \brief Deletes an empty interaction
 * \details When handling the Interaction MPI data type, an empty interaction is created.
 * This function destroys that empty function again
 * \param[in,out] interaction Pointer to a base interaction that needs to be deleted
 */
void InteractionHandler::deleteEmptyInteraction(BaseInteraction* interaction) const
{
    BaseSpecies* species = this->getDPMBase()->speciesHandler.getObject(0);
    species->deleteEmptyInteraction(interaction);
}

/*!
 * \brief creates an empty MPIInteractionDataArray
 * \details When sending interactions to other processors it is done with a void* array, because
 * interactions come in many different sizes and shapes. The size of an interaction depends on the type f
 * interaction. By creating an emptyInteraction, the size can be determined and an empty vector can be allocated
 * \param[in] numberOfInteractions The number of interactions that are added to the MPIInteractionDataArray
 * \return Returns a void pointer to a block of memory that will contain the interaction data
 */
void* InteractionHandler::createMPIInteractionDataArray(unsigned int numberOfInteractions) const
{
    BaseInteraction* emptyInteraction = this->createEmptyInteraction();
    //create a vector based on the first interaction. This is to avoid templating the whole class
    void* historyData = emptyInteraction->createMPIInteractionDataArray(numberOfInteractions);
    this->deleteEmptyInteraction(emptyInteraction);
    
    return historyData;
}

/*!
 * \brief deletes an MPIInteractionDataArray
 * \details An empty interaction is created that can recast the data array into an array of the
 * interactions. After the conversion the data of the array is cleaned up
 * \param[in] dataArray void pointer to a data array
 */
void InteractionHandler::deleteMPIInteractionDataArray(void* dataArray)
{
    BaseInteraction* emptyInteraction = this->createEmptyInteraction();
    emptyInteraction->deleteMPIInteractionDataArray(dataArray);
    this->deleteEmptyInteraction(emptyInteraction);
}

/*!
 * \brief reads the basic interaction details from an MPIInteractionDataArray
 * \details because the interactionData array is a void pointer array, it is not possibly to directly read out the interaction data
 * this function reads the basic information such as which objects are interacting and at what time.
 * \param[in] interactionData void pointer to the mpi interaction data
 * \param[in] index The index of the interaction in the interaction data
 * \param[out] identificationP the unique id of the P object
 * \param[out] idientificationI the unique id of the I object
 * \param[out] isWallInteraction a bool that flags if the interaction is a wall interaction or not
 * \param[out] timeStamp reads the timestamp of the interaction
 */
void InteractionHandler::getInteractionDetails(void* interactionData, unsigned int index, unsigned int& identificationP,
                                               unsigned int& identificationI, bool& isWallInteraction,
                                               unsigned& timeStamp)
{
    BaseInteraction* emptyInteraction = this->createEmptyInteraction();
    emptyInteraction->getInteractionDetails(interactionData, index, identificationP, identificationI, isWallInteraction,
                                            timeStamp);
    this->deleteEmptyInteraction(emptyInteraction);
}

/*!
 * \details Returns a pointer to the existing Interaction, if the Interaction
 * already exists otherwise creates a new Interaction and returns a pointer to it.
 * \param[in] P the first BaseInteractable by which the interaction is defined.
 * \param[in] I the first BaseInteractable by which the interaction is defined.
 * \param[in] timeStamp the current value of DPMBase::time_.
 * \return the Interaction between the BaseInteractable's P and I
 */
BaseInteraction*
InteractionHandler::getInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp, const Vec3D& normal)
{
    BaseInteraction* c = nullptr;
    Mdouble dNormalMax = 0.9; //minimal agreement required
    for (const auto& i : P->getInteractions())
    {
        if (i->getI() == I && i->getTimeStamp() != timeStamp)
        {
            const Mdouble dNormal = Vec3D::dot(normal, i->getNormal());
            if (dNormal > dNormalMax)
            {
                c = i;
                dNormalMax = dNormal;
            }
        }
    }
    
    BaseSpecies* species = getDPMBase()->speciesHandler.getMixedObject(P->getIndSpecies(), I->getIndSpecies());
    
    if (c)
    {
        c->setTimeStamp(timeStamp);
    }
    else
    {
        c = species->getNewInteraction(P, I, timeStamp);
        c->setSpecies(species);
        addObject(c);
        //logger(INFO,"new interaction t=%: i=%, n=%",timeStamp,getNumberOfObjects(),normal);
    }
    
    return c;
}

/*!
 * \details Deleting the three periodic interactions between two real particles 
 * is difficult, because its interaction information has to be saved. 
 * If the two real particles interacted (which can be checked by looking at the 
 * time stamp), the interaction between the real particles is kept, and all 
 * interactions that involve ghost particles gets removed; 
 * otherwise, the interaction between the lower-indexed real particle with the 
 * ghost particle of the higher indexed particles is saved (with the ghost 
 * particle replaced by the real particle), and all other interactions removed.
 * 
 * This is what this function is intended for, and it does it in the following way:
 * When an interaction is removed the periodic particle has to be stored in the I pointer
 * So when an interaction is removed where P is normal and I is periodic, 
 * and the information is new it will be transfered when the index of P is 
 * lower than the index of the real particle of I. 
 * \image html Walls/periodicBoundary.pdf
 * \param[in] the id of the Interaction that needs to be deleted.
 */
void InteractionHandler::removeObjectKeepingPeriodics(const unsigned int id)
{
    BaseInteraction* iMain = getObject(id);
    
    BaseParticle* P = dynamic_cast<BaseParticle*>(iMain->getP());
    BaseParticle* I = dynamic_cast<BaseParticle*>(iMain->getI());
    if (P != nullptr && I != nullptr) //check that both P and I are particles (not walls)
    {
        BaseInteractable* realI = I->getPeriodicFromParticle();
        if (realI != nullptr && !P->getPeriodicFromParticle()) //check that P is a real and I is a ghost particle
        {
            if (P->getIndex() < realI->getIndex())
            {
                BaseInteraction* iOther = getExistingInteraction(P, realI);
                //You have to also check for existing interactions of the particles in reverse order since the copySwitchPointer function can revert the order of the particles
                ///\todo TW The code assumes in a few places that P->getIndex()<I->getIndex(), but the copySwitchPointer function does not obey that rule; we have to check if this is valid behaviour.
                if (iOther == nullptr)
                {
                    iOther = getExistingInteraction(realI, P);
                }
                if (iOther != nullptr) //if the interaction existed before the ghost particles were created
                {
                    //Here we decide which of the two interactions should be kept:
                    //the interaction between the two real particles (iMain), or
                    //the interaction between the real and a ghost particle (iOther).
                    //It picks the one for which a collision has happened,
                    //i.e. the one with the newer timeStamp.
                    ///\todo this function will create an error if the timeStamp is in the future! This should not happen (ever), but who knows.
                    if (iOther->getTimeStamp() <
                        iMain->getTimeStamp()) //if the interaction has been active during the last computeForce routine, make this the new (real) interaction.
                    
                    {
                        iMain->setI(realI);
                        removeObject(iOther->getIndex());
                        return;
                    }
                    else  //if the interaction has not been active during the last computeForce routine
                    {
                        BaseHandler<BaseInteraction>::removeObject(id);
                        return;
                    }
                }
                else //if the interaction has been created during the last computeForce routine, make this a new (real) interaction. 
                {
                    iMain->setI(realI);
                    return;
                }
            }
        }
    }
    //this point is reached if either P or I are walls, or P and I are both ghost particles; in these cases, the interaction gets deleted
    BaseHandler<BaseInteraction>::removeObject(id);
}

/*!
 * \details Each interaction contains a time stamp, which stores the last time 
 * that an interaction object has been called. Thus, one can see that an 
 * interaction has ended by comparing the time stamp with the last value of DPMBase::time_.
 * This function erases all interactions that have ended. 
 * \param[in] lastTimeStep the last used value of DPMBase::time_.
 */
void InteractionHandler::eraseOldInteractions(unsigned currentNTime)
{
    //std::cout<<"void InteractionHandler::eraseOldInteractions(Mdouble lastTimeStep)"<<std::endl;
    //std::cout<<"Current interactions="<<getNumberOfObjects()<<std::endl;
    //Remove_if reconstructs the vector with only elements passing the check_spring_time function
    //Erase removes the end of the vector
    ///\todo TW: this function has to be sped up with sth like this: erase(remove_if(begin(), end(), bind2nd(checkSpringTime(), lastTimeStep)), end());
    for (unsigned int id = 0; id < getNumberOfObjects(); id++)
    {
        if (getObject(id)->getTimeStamp() <= currentNTime)
        {
            getObject(id)->actionsOnErase();
            removeObject(id);
            --id;
        }
    }
}

void InteractionHandler::actionsAfterTimeStep()
{
    for (auto i: *this)
    {
        i->actionsAfterTimeStep();
    }
}

/*!
 * \return      The mean overlap of all interactions
 */
Mdouble InteractionHandler::getMeanOverlap() const
{
    Mdouble sum = 0;
    Mdouble n = 0;
    for (BaseInteraction* const p : objects_)
    {
        if (p->getOverlap() > 0)
        {
            sum += p->getOverlap();
            n++;
        }
    }
    return sum / n;
}

/*!
 * \return the string InteractionHandler
 */
std::string InteractionHandler::getName() const
{
    return "InteractionHandler";
}

/*!
 * \param[in] os The output stream where the InteractionHandler must be written
 *  to, usually a restart file.
 */
void InteractionHandler::write(std::ostream& os) const
{
#ifdef MERCURY_USE_MPI
    //note: this function only prints the particles on processor 0
    //The rest of the particles are processed in the restart file write function
    MPIContainer& communicator = MPIContainer::Instance();
    unsigned int totalNumberOfInteractions = getNumberOfObjects();
    os << "Interactions " << totalNumberOfInteractions << '\n';
    //logger(INFO,"% interactions",getNumberOfObjects());
    for (BaseInteraction* it : *this)
    {
    os << (*it) << '\n';
    }
#else
    os << "Interactions " << getNumberOfObjects() << '\n';
    for (BaseInteraction* i : objects_)
        os << (*i) << '\n';
#endif
}

/// \param[in] is The input stream from which the information is read.
void InteractionHandler::read(std::istream& is)
{
    clear();
    unsigned int N;
    std::string dummy;
    is >> dummy;
    std::stringstream line;
    helpers::getLineFromStringStream(is, line);
    line >> N;
    if (N > 1e5) logger(INFO, "Reading % %", N, dummy);
    logger(VERBOSE, "In %::read(is): reading in % objects.", getName(), N);
    setStorageCapacity(N);
    // set map
    particleById.clear();
    for (BaseParticle* p : getDPMBase()->particleHandler) {
        particleById[p->getId()] = p;
    }
    wallById.clear();
    for (BaseWall* w : getDPMBase()->wallHandler) {
        wallById[w->getId()] = w;
    }
    for (unsigned int i = 0; i < N; i++) {
        readAndAddObject(is);
    }
    particleById.clear();
    wallById.clear();
}


/*!
 * \param[in] is The input stream from which the information is read.
 */
void InteractionHandler::readAndAddObject(std::istream& is)
{
    std::string type, dummy, idType;
    unsigned int id0, id1;
    Mdouble doubleTimeStamp;
    unsigned timeStamp;
    
    /// \todo Ant This is a tmp fix as in some cases the line before has not be finished reading. This should be looked at again at a later date.
    is >> type;
    logger(DEBUG, "InteractionHandler::readAndAddObject(is): reading type %.", type);
    Mdouble timeStampDouble;
    is >> idType >> id0 >> id1 >> dummy >> timeStampDouble;
    timeStamp = timeStampDouble; //in order to read old restart files

    ///\todo TW: Change identifier in restart file from id to index; is there any reason the id should be kept after restarting, once this is done? (Note, the id is set to the old one in the particle handler because interactions store id, not indices; also note id's are slow
    BaseInteraction* C;

#ifdef MERCURY_USE_MPI
    if (idType.compare("particleIds") == 0)
    {
        std::vector<BaseParticle*> list0 = getDPMBase()->particleHandler.getObjectsById(id0);
        std::vector<BaseParticle*> list1 = getDPMBase()->particleHandler.getObjectsById(id1);
        for (int p0 = 0; p0 < list0.size(); p0++)
        {
            for (int p1  = 0; p1 < list1.size(); p1++)
            {
                C = getInteraction(list0[p0],list1[p1], timeStamp);
                if (C != nullptr)
                {
                   is >> *C;
                }
            }
        }

    }
    else
    {
        std::vector<BaseParticle*> list0 = getDPMBase()->particleHandler.getObjectsById(id0);
        for (int p0 = 0; p0 < list0.size(); p0++)
        {
            C = getInteraction(list0[p0],getDPMBase()->wallHandler.getObjectById(id1), timeStamp);
            if (C != nullptr)
            {
                is >> *C;
            }
        }
    }
#else
    //this requires particleById/wallById to be set, which is done in the read() function
    if (idType == "particleIds")
        C = getInteraction(particleById[id0],particleById[id1],timeStamp);
    else
        C = getInteraction(particleById[id0],wallById[id1], timeStamp);
    is >> (*C);
#endif
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void InteractionHandler::setWriteVTK(FileType fileType)
{
    writeVTK_ = fileType;
}

FileType InteractionHandler::getWriteVTK() const
{
    return writeVTK_;
}

double InteractionHandler::getLiquidBridgeVolume() const {
    double liquidVolume = 0;
    for (auto i : objects_) {
        auto j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
        if (j and !static_cast<BaseParticle*>(j->getP())->isMPIParticle()) liquidVolume += j->getLiquidBridgeVolume();
    }
    return getMPISum(liquidVolume);
};
