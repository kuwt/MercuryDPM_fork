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


#include "EmptyAdhesiveInteraction.h"
#include "Species/AdhesiveForceSpecies/EmptyAdhesiveSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
EmptyAdhesiveInteraction::EmptyAdhesiveInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"EmptyAdhesiveInteraction::EmptyAdhesiveInteraction() finished"<<std::endl;
#endif
}

EmptyAdhesiveInteraction::EmptyAdhesiveInteraction()
        : BaseInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"EmptyAdhesiveInteraction::EmptyAdhesiveInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] P
 */
EmptyAdhesiveInteraction::EmptyAdhesiveInteraction(const EmptyAdhesiveInteraction& p)
        : BaseInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"EmptyAdhesiveInteraction::EmptyAdhesiveInteraction(const EmptyAdhesiveInteraction &p finished"<<std::endl;
#endif
}

/*!
 *
 */
EmptyAdhesiveInteraction::~EmptyAdhesiveInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"EmptyAdhesiveInteraction::~EmptyAdhesiveInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void EmptyAdhesiveInteraction::write(std::ostream& os UNUSED) const
{}

/*!
 * \param[in,out] is
 */
void EmptyAdhesiveInteraction::read(std::istream& is UNUSED)
{}

/*!
 *
 */
void EmptyAdhesiveInteraction::computeAdhesionForce()
{}

void EmptyAdhesiveInteraction::computeAdhesionInteraction()
{}

/*!
 * \return Mdouble
 */
Mdouble EmptyAdhesiveInteraction::getElasticEnergy() const
{
    return 0.0;
}

/*!
 * \return const EmptyAdhesiveSpecies*
 */
const EmptyAdhesiveSpecies* EmptyAdhesiveInteraction::getSpecies() const
{
    return static_cast<const EmptyAdhesiveSpecies*>(getBaseSpecies()->getAdhesiveForce());
;
}

/*!
 * \return std::string
 */
std::string EmptyAdhesiveInteraction::getBaseName() const
{
    return "";
}
