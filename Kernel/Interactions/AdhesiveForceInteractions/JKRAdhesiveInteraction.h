//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef LINEARJKRADHESIVEINTERACTION_H
#define LINEARJKRADHESIVEINTERACTION_H

#include "Interactions/BaseInteraction.h"
#include "Math/Vector.h"

class BaseParticle;

class JKRAdhesiveSpecies;

class BaseInteractable;

/*!
 * \class ReversibleAdheseiveInteraction
 * \brief Computes the interactions between particles for reversive adhesive contact model.
 */
class JKRAdhesiveInteraction : public virtual BaseInteraction
{
public:
    /*!
     * \brief Setting an alias name for JKRAdhesiveSpecies.
     */
    typedef JKRAdhesiveSpecies SpeciesType;

    /*!
     * \brief Constructor
     */
    JKRAdhesiveInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);

    //used for mpi
    JKRAdhesiveInteraction();

    /*!
     * \brief Copy constructor.
     */
    JKRAdhesiveInteraction(const JKRAdhesiveInteraction& p);

    /*!
     * \brief Destructor.
     */
    ~JKRAdhesiveInteraction() override;

    /*!
     * \brief Computes the adhesive forces
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
    /*!
     * \brief Returns the elastic energy stored in the adhesive spring.
     */
    Mdouble getElasticEnergy() const override;

    /*!
     * \brief A dynamic_cast of BaseSpecies pointer type to a pointer to an object of
     *        type JKRAdhesiveSpecies.
     */
    const JKRAdhesiveSpecies* getSpecies() const;

    /*!
     * \brief Returns the name of the interaction, see Interaction.h.
     */
    std::string getBaseName() const;
};

#endif
