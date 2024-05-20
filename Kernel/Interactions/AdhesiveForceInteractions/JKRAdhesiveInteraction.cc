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


#include "JKRAdhesiveInteraction.h"
#include "Species/AdhesiveForceSpecies/JKRAdhesiveSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
JKRAdhesiveInteraction::JKRAdhesiveInteraction(BaseInteractable* P, BaseInteractable* I,
                                                             unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"JKRAdhesiveInteraction::JKRAdhesiveInteraction() finished"<<std::endl;
#endif
}


//used for mpi
JKRAdhesiveInteraction::JKRAdhesiveInteraction()
        : BaseInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"JKRAdhesiveInteraction::JKRAdhesiveInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
JKRAdhesiveInteraction::JKRAdhesiveInteraction(const JKRAdhesiveInteraction& p)
        : BaseInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"JKRAdhesiveInteraction::JKRAdhesiveInteraction(const JKRAdhesiveInteraction &p finished"<<std::endl;
#endif
}

/*!
 *
 */
JKRAdhesiveInteraction::~JKRAdhesiveInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"JKRAdhesiveInteraction::~JKRAdhesiveInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] os
 */
void JKRAdhesiveInteraction::write(std::ostream& os  UNUSED) const
{}

/*!
 * \param[in] is
 */
void JKRAdhesiveInteraction::read(std::istream& is  UNUSED)
{}

/*!
 *
 */
void JKRAdhesiveInteraction::computeAdhesionForce()
{
    const JKRAdhesiveSpecies* species = getSpecies();
    Mdouble JKRPullOffForce = (3.0/2.0)*constants::pi*species->getSurfaceEnergy()*getEffectiveRadius();

    if (getOverlap() >= 0)
        addForce(getNormal() * (-JKRPullOffForce));
    else if (getOverlap() >= (-JKRPullOffForce/species->getAdhesionStiffness()) && getOverlap() < 0 )
        addForce(getNormal() * (-species->getAdhesionStiffness() * getOverlap() - JKRPullOffForce));
    else
        addForce(Vec3D(0.0,0.0,0.0));
}

/*!
 * \todo
 * \return Mdouble
 */
Mdouble JKRAdhesiveInteraction::getElasticEnergy() const
{
    return 0.0;
}

/*!
 * \return a constant pointer to an instance of this class.
 */
const JKRAdhesiveSpecies* JKRAdhesiveInteraction::getSpecies() const
{
    return dynamic_cast<const JKRAdhesiveSpecies*> (getBaseSpecies()); //downcast
}

/*!
 * \return std::string
 */
std::string JKRAdhesiveInteraction::getBaseName() const
{
    return "JKRAdhesive";
}
