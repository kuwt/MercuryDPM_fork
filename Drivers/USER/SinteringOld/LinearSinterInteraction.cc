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


#include "LinearSinterInteraction.h"
#include "LinearSinterSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <DPMBase.h>

LinearSinterInteraction::LinearSinterInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
    : BaseInteraction(P, I, timeStamp),  LinearPlasticViscoelasticInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterInteraction::LinearSinterInteraction() finished"<<std::endl;
#endif
}

LinearSinterInteraction::LinearSinterInteraction(const LinearSinterInteraction &p)
    : BaseInteraction(p),  LinearPlasticViscoelasticInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterInteraction::LinearSinterInteraction(const LinearSinterInteraction &p finished"<<std::endl;
#endif
}

LinearSinterInteraction::LinearSinterInteraction()
    : BaseInteraction(),  LinearPlasticViscoelasticInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterInteraction::LinearSinterInteraction() finished"<<std::endl;
#endif
}

LinearSinterInteraction::~LinearSinterInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterInteraction::~LinearSinterInteraction() finished"<<std::endl;
#endif
}

void LinearSinterInteraction::computeForce()
{
    LinearPlasticViscoelasticInteraction::computeForce();
    setMaxOverlap(getMaxOverlap() + getSpecies()->getNeckGrowthRate() * getHandler()->getDPMBase()->getTimeStep());
}

const LinearSinterSpecies *LinearSinterInteraction::getSpecies() const
{
    return dynamic_cast<const LinearSinterSpecies *>(getBaseSpecies());
}

std::string LinearSinterInteraction::getName() const
{
    return "LinearSinterSpecies";
}
