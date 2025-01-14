//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


#include "LiquidBridgeWilletInteraction.h"
#include "Species/AdhesiveForceSpecies/LiquidBridgeWilletSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
LiquidBridgeWilletInteraction::LiquidBridgeWilletInteraction(BaseInteractable* P, BaseInteractable* I,
                                                             unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeWilletInteraction::LiquidBridgeWilletInteraction() finished"<<std::endl;
#endif
}

//used for mpi
LiquidBridgeWilletInteraction::LiquidBridgeWilletInteraction()
        : BaseInteraction()
{
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeWilletInteraction::LiquidBridgeWilletInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
LiquidBridgeWilletInteraction::LiquidBridgeWilletInteraction(const LiquidBridgeWilletInteraction& p)
        : BaseInteraction(p)
{
    wasInContact_ = p.wasInContact_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidBridgeWilletInteraction::LiquidBridgeWilletInteraction(const LiquidBridgeWilletInteraction &p finished"<<std::endl;
#endif
}

/*!
 * 
 */
LiquidBridgeWilletInteraction::~LiquidBridgeWilletInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LiquidBridgeWilletInteraction::~LiquidBridgeWilletInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void LiquidBridgeWilletInteraction::write(std::ostream& os) const
{
    os << " wasInContact " << wasInContact_;
}

/*!
 * \param[in,out] is
 */
void LiquidBridgeWilletInteraction::read(std::istream& is)
{
    // we read in wasInContact_ like this because an early version did not initialize it.
    std::string dummy;
    is >> dummy >> wasInContact_;
}

/*!
 * 
 */
void LiquidBridgeWilletInteraction::computeAdhesionForce()
{
    const LiquidBridgeWilletSpecies* species = getSpecies();
    if (getOverlap() >= 0)
    {
        wasInContact_ = true;
        Mdouble effectiveRadius = 2.0 * getEffectiveRadius();
        Mdouble fdotn = -2.0 * constants::pi * effectiveRadius * species->getSurfaceTension() *
                        std::cos(species->getContactAngle());
        addForce(getNormal() * fdotn);
    }
    else if (wasInContact_)
    {
        Mdouble effectiveRadius = 2.0 * getEffectiveRadius();
        Mdouble s_c = -getOverlap() * std::sqrt(effectiveRadius / species->getLiquidBridgeVolume());
        Mdouble fdotn = -2.0 * constants::pi * effectiveRadius * species->getSurfaceTension()
                        * std::cos(species->getContactAngle()) / (1 + (1.05 + 2.5 * s_c) * s_c);
        addForce(getNormal() * fdotn);
    }
}

/*!
 * \return Mdouble
 */
Mdouble LiquidBridgeWilletInteraction::getElasticEnergy() const
{
    ///\todo TW
    return 0.0;
}

/*!
 * \return const LiquidBridgeWilletSpecies*
 */
const LiquidBridgeWilletSpecies* LiquidBridgeWilletInteraction::getSpecies() const
{
    return static_cast<const LiquidBridgeWilletSpecies*>(getBaseSpecies()->getAdhesiveForce());
;
}

/*!
 * \return std::string
 */
std::string LiquidBridgeWilletInteraction::getBaseName() const
{
    return "LiquidBridgeWillet";
}

bool LiquidBridgeWilletInteraction::getWasInContact() const
{
    return wasInContact_;
}

void LiquidBridgeWilletInteraction::setWasInContact(bool wasInContact)
{
    wasInContact_ = wasInContact;
}
