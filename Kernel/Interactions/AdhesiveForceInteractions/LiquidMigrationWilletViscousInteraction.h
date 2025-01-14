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

#ifndef LiquidMigrationWilletViscousInteraction_H
#define LiquidMigrationWilletViscousInteraction_H

#include "Interactions/BaseInteraction.h"
#include "Math/Vector.h"
#include "ParticleHandler.h"
#include "InteractionHandler.h"
#include "LiquidMigrationWilletInteraction.h"

class BaseInteractable;
class LiquidMigrationWilletViscousSpecies;


/*!
 * \class LiquidMigrationWilletViscousInteraction
 * \brief Defines the liquid bridge willetViscous interaction between two particles or walls.
 */
class LiquidMigrationWilletViscousInteraction : public LiquidMigrationWilletInteraction
{
public:
    /*!
     * \brief An alias name for LiquidMigrationWilletViscousSpecies data type.
     */
    typedef LiquidMigrationWilletViscousSpecies SpeciesType;

    /*!
     * \brief Constructor.
     */
    LiquidMigrationWilletViscousInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);

    //used for mpi
    LiquidMigrationWilletViscousInteraction();

    /*!
     * \brief Copy constructor.
     */
    LiquidMigrationWilletViscousInteraction(const LiquidMigrationWilletViscousInteraction& p);

    /*!
     * \brief Destructor.
     */
    ~LiquidMigrationWilletViscousInteraction() override;

    /*!
     * \brief Returns the name of the interaction, see Interaction.h.
     */
    std::string getBaseName() const;

    /*!
     * \brief Returns a pointer to the adhesive force species LiquidMigrationWilletViscousSpecies.
     */
    const LiquidMigrationWilletViscousSpecies* getSpecies() const;

    /*!
     * \brief Computes the adhesive forces for liquid bridge Willet and viscous type of interaction.
     */
    void computeAdhesionForce();

    /*!
     * \brief Accesses the minimum distance that the viscous liquid force is valid.
     */
    Mdouble getLimitingDistance();



};

#endif
