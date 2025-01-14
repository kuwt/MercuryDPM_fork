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

#ifndef HeatFluidCoupledINTERACTION_H
#define HeatFluidCoupledINTERACTION_H

#include "ThermalInteraction.h"
//#include "Species/NormalForceSpecies/HeatFluidCoupledSpecies.h"
#include "Particles/HeatFluidCoupledParticle.h"

template<class NormalForceSpecies>
class HeatFluidCoupledSpecies;

template<class NormalForceInteraction>
class HeatFluidCoupledInteraction : public ThermalInteraction<NormalForceInteraction>
{
public:
    typedef HeatFluidCoupledSpecies<typename NormalForceInteraction::SpeciesType> SpeciesType;
    
    /*!
     * \brief Constructor.
     */
    HeatFluidCoupledInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
    : BaseInteraction(P, I, timeStamp), ThermalInteraction<NormalForceInteraction>(P, I, timeStamp)
    {}
    
    /*!
     * \brief Default Constructor.
     */
    HeatFluidCoupledInteraction()
    : BaseInteraction(), ThermalInteraction<NormalForceInteraction>()
    {}


    /*!
     * \brief Copy constructor.
     */
    HeatFluidCoupledInteraction(const HeatFluidCoupledInteraction& p)
    : BaseInteraction(p), ThermalInteraction<NormalForceInteraction>(p)
    {}

    /*!
     * \brief Destructor.
     */
    virtual ~HeatFluidCoupledInteraction()
    {}
    
    /*!
     * \brief Computes the normal forces due to linear plastic visco elastic interaction.
     */
    void computeNormalForce();
};

template<class NormalForceInteraction>
void HeatFluidCoupledInteraction<NormalForceInteraction>::computeNormalForce()
{
    NormalForceInteraction::computeNormalForce();
    Mdouble radius = 2.0 * NormalForceInteraction::getEffectiveRadius();
    Mdouble contactArea = constants::pi * radius * std::max(0.0,NormalForceInteraction::getOverlap());
    const SpeciesType* species = static_cast<const SpeciesType*>(NormalForceInteraction::getBaseSpecies()->getNormalForce());
    auto pParticle = dynamic_cast<HeatFluidCoupledParticle*>(NormalForceInteraction::getP());
    auto iParticle = dynamic_cast<HeatFluidCoupledParticle*>(NormalForceInteraction::getI());
    // if both p and i are particles
    if (pParticle && iParticle)
    {
        /* heat transfer rate Q=m*c*dT/dt */
        Mdouble heatTransfer = species->getThermalConductivity()
                               * (pParticle->getTemperature() - iParticle->getTemperature())
                               * contactArea / NormalForceInteraction::getDistance();
        /* m*dT = Q/c*dt */
        Mdouble mdT = heatTransfer / species->getHeatCapacity()
                      * NormalForceInteraction::getHandler()->getDPMBase()->getTimeStep();
        pParticle->addTemperature(-mdT * pParticle->getInvMass());
        iParticle->addTemperature(mdT * iParticle->getInvMass());
    }
}

#endif
