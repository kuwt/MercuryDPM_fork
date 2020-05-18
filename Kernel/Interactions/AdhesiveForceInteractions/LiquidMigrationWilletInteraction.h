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

#ifndef LiquidMigrationWilletInteraction_H
#define LiquidMigrationWilletInteraction_H

#include "Interactions/BaseInteraction.h"
#include "Math/Vector.h"
#include "ParticleHandler.h"
#include "InteractionHandler.h"

class BaseParticle;

class LiquidMigrationWilletSpecies;

class BaseInteractable;

/*!
 * \class LiquidMigrationWilletInteraction
 * \brief Defines the liquid bridge willet interaction between two particles or walls.
 */
class LiquidMigrationWilletInteraction : public virtual BaseInteraction
{
public:
    /*!
     * \brief An alias name for LiquidMigrationWilletSpecies data type.
     */
    typedef LiquidMigrationWilletSpecies SpeciesType;
    
    /*!
     * \brief Constructor.
     */
    LiquidMigrationWilletInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    //used for mpi
    LiquidMigrationWilletInteraction();
    
    /*!
     * \brief Copy constructor.
     */
    LiquidMigrationWilletInteraction(const LiquidMigrationWilletInteraction& p);
    
    /*!
     * \brief Destructor.
     */
    ~LiquidMigrationWilletInteraction() override;
    
    void actionsOnErase() override;
    
    void actionsAfterTimeStep() override;
    
    /*!
     * \brief Computes the adhesive forces for liquid bridge Willet type of interaction.
     */
    void computeAdhesionForce();
    
    /*!
     * \brief Interaction read function, which accepts an std::istream as input.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Interaction print function, which accepts an std::ostream as input.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the amount of Elastic energy involved in an interaction. Basically
     *        used in case you want to write the elastic energy into an output file. 
     */
    Mdouble getElasticEnergy() const override;
    
    /*!
     * \brief Returns a pointer to the adhesive force species LiquidMigrationWilletSpecies.
     */
    const LiquidMigrationWilletSpecies* getSpecies() const;
    
    /*!
     * \brief Returns the name of the interaction, see Interaction.h.
     */
    std::string getBaseName() const;
    
    Mdouble getLiquidBridgeVolume() const;
    
    void setLiquidBridgeVolume(Mdouble liquidBridgeVolume);
    
    void addLiquidBridgeVolume(Mdouble liquidBridgeVolume);
    
    bool getWasInContact() const;
    
    void setWasInContact(bool wasInContact);
    
    void rupture();
    
    void form();
    
    Mdouble getRuptureDistance();
    
    int getNumberOfContacts(BaseInteractable* interactable);
    
    /**
     * writes extra information to the VTK output
     */
    unsigned getNumberOfFieldsVTK() const override;
    
    std::string getTypeVTK(unsigned i) const override;
    
    std::string getNameVTK(unsigned i) const override;
    
    std::vector<Mdouble> getFieldVTK(unsigned i) const override;
    
    static Mdouble getTotalLiquidFilmVolume(ParticleHandler&);
    
    static Mdouble getTotalLiquidBridgeVolume(InteractionHandler&);

private:
    /*!
     * \brief A history parameter to store if the particles were in contact or not. Useful
     *        to compute adhesive forces.
     */
    bool wasInContact_;
    
    Mdouble liquidBridgeVolume_;
};

#endif
