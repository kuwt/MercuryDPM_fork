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
#ifndef INTERACTIONHANDLER_H
#define INTERACTIONHANDLER_H

#include <Interactions/AdhesiveForceInteractions/LiquidBridgeWilletInteraction.h>
#include <Interactions/AdhesiveForceInteractions/LiquidMigrationWilletInteraction.h>
#include "BaseHandler.h"
#include "Interactions/BaseInteraction.h"
#include "File.h"
#include "MpiContainer.h"
#include <map>

class SpeciesHandler;

/*!
 * \brief Container to store Interaction objects.
 * \details The InteractionHandler is a container to store all Interaction objects. 
 * It is implemented as a vector of BaseInteraction pointers.
 */
class InteractionHandler final : public BaseHandler<BaseInteraction>
{
public:
    /*!
     * \brief Default constructor, it creates an empty InteractionHandler.
     */
    InteractionHandler();
    
    /*!
     * \brief Copy constructor, but since interactions must not be copied, it creates an empty InteractionHandler.
     */
    InteractionHandler(const InteractionHandler& IH UNUSED);
    
    /*!
     * \brief Assignment operator.
     */
    InteractionHandler& operator=(const InteractionHandler& rhs);
    
    /*!
     * \brief Destructor, it destructs the InteractionHandler and all BaseInteraction it contains.
     */
    ~InteractionHandler() override;
    
    /*!
     * \brief Adds an Interaction to the InteractionHandler.
     */
    void addObject(BaseInteraction* I) override;

    /// \param[in] is The input stream from which the information is read.
    void read(std::istream& is);

    /*!
     * \brief Reads an Interaction into the InteractionHandler from restart data.
     */
    void readAndAddObject(std::istream& is) override;
    
    /*!
     * \brief Returns the Interaction between the BaseInteractable's P and I if 
     * it exists, otherwise returns a null pointer.
     */
    BaseInteraction* getExistingInteraction(const BaseInteractable* P, const BaseInteractable* I) const;
    
    BaseInteraction* addInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    /*
     * This function is used in omp simulations to add new objects to the interactionHandler during computeOneTimeStep.
     * It is step 1 of the add-by-reduction approach:
     *  1. Declare an empty vector newObjects_, one element per thread.
     *  2. When a new object needs to be added to the interactionHandler, don't add it immediately.
     *     Add it to newObjects_[threadNum] instead
     *  3. Add all objects in newObjects_ to teh interactionHandler.
     * This reduction approach allows step 2 to be run in parallel, while step 3 needs to be done in sequence.
     * It is efficient, since step 3 is much quicker than step 2.
     */
    void resetNewObjectsOMP();
    
    /*
     * This function is used in omp simulations to add new objects to the interactionHandler during computeOneTimeStep.
     * It is step 3 of the add-by-reduction approach:
     *  1. Declare an empty vector newObjects_, one element per thread.
     *  2. When a new object needs to be added to the interactionHandler, don't add it immediately.
     *     Add it to newObjects_[threadNum] instead
     *  3. Add all objects in newObjects_ to teh interactionHandler.
     * This reduction approach allows step 2 to be run in parallel, while step 3 needs to be done in sequence.
     * It is efficient, since step 3 is much quicker than step 2.
     */
    void addNewObjectsOMP();
    
    /*!
     * \brief Returns the Interaction between the BaseInteractable's P and I.
     */
    BaseInteraction* getInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    /*!
     * \brief Creates an empty interaction
     */
    //BaseInteraction* getExistingInteraction(BaseInteractable* P, BaseInteractable* I) const;
    
    /*!
     * \brief Returns the Interaction between the BaseInteractable's P and I.
     */
    BaseInteraction* createEmptyInteraction() const;
    
    /*!
     * \brief Deletes an empty interaction
     */
    void deleteEmptyInteraction(BaseInteraction* interaction) const;
    
    /*!
     * \brief creates an empty MPIInteractionDataArray
     */
    void* createMPIInteractionDataArray(unsigned int numberOfInteractions) const;
    
    /*!
     * \brief deletes an MPIInteractionDataArray
     */
    void deleteMPIInteractionDataArray(void* dataArray);
    
    /*!
     * \brief reads the basic interaction details from an MPIInteractionDataArray
     */
    void getInteractionDetails(void* interactionData, unsigned int index, unsigned int& identificationP,
                               unsigned int& identificationI, bool& isWallInteraction, unsigned& timeStamp);

    /*!
     * \brief Returns the Interaction between the BaseInteractable's P and I closest to the contact point (should be used when multiple contacts are possible).
     */
    BaseInteraction* getInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp, const Vec3D& normal);
    
    /*!
     * \brief Removes interactions of periodic particles when the periodic 
     * particles get deleted (see DPMBase::removeDuplicatePeriodicParticles) 
     */
    void removeObjectKeepingPeriodics(unsigned int id);
    
    /*!
     * \brief erases interactions which have an old timestamp.
     */
    void eraseOldInteractions(unsigned);
    
    void actionsAfterTimeStep();
    
    /*!
    * \brief The mean overlap of all interactions
    */
    Mdouble getMeanOverlap() const;
    
    /*!
     * \brief Writes the InteractionHandler to an output stream, for example a restart file.
     */
    void write(std::ostream& os) const;
    
    /*!
     * \brief
     */
    void setWriteVTK(FileType f);
    
    /*!
     * \brief
     */
    FileType getWriteVTK() const;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*
     * This variable is used in omp simulations to add new objects to the interactionHandler during computeOneTimeStep.
     * Objects are added in three steps:
     *  1. Declare an empty vector newObjects_, one element per thread.
     *  2. When a new object needs to be added to the interactionHandler, don't add it immediately.
     *     Add it to newObjects_[threadNum] instead
     *  3. Add all objects in newObjects_ to teh interactionHandler.
     * This reduction approach allows step 2 to be run in parallel, while step 3 needs to be done in sequence.
     * It is efficient, since step 3 is much quicker than step 2.
     */
    std::vector<std::vector<BaseInteraction*>> newObjects_;

    double getLiquidBridgeVolume() const;

private:
    
    FileType writeVTK_;

    std::map<unsigned,BaseParticle*> particleById;
    std::map<unsigned,BaseWall*> wallById;
};

#endif

