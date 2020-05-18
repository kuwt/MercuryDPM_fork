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

#ifndef PARHAMIMCMEEKINGSINTERSPECIES_H
#define PARHAMIMCMEEKINGSINTERSPECIES_H

#include "Species/AdhesiveForceSpecies/BaseAdhesiveForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/AdhesiveForceInteractions/ParhamiMcMeekingSinterInteraction.h"

/*!
 * \brief ParhamiMcMeekingSinterSpecies contains the parameters used to describe a linear reversible short-range force.
 * \details See ParhamiMcMeekingSinterInteraction::computeForce for a description of the force law.
 */
class ParhamiMcMeekingSinterSpecies : public BaseAdhesiveForce
{
public:
    ///\brief The correct Interaction type for this AdhesiveForceSpecies
    typedef ParhamiMcMeekingSinterInteraction InteractionType;
    
    ///\brief The default constructor.
    ParhamiMcMeekingSinterSpecies();
    
    ///\brief The default constructor.
    ParhamiMcMeekingSinterSpecies(const ParhamiMcMeekingSinterSpecies& s);
    
    ///\brief The default constructor.
    ~ParhamiMcMeekingSinterSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;
    
    ///\brief creates default values for mixed species
    void mix(ParhamiMcMeekingSinterSpecies* S, ParhamiMcMeekingSinterSpecies* T);

//adhesion-specific functions

    //setters and getters
    void set(Mdouble alpha, Mdouble beta, Mdouble atomicVolume /*Omega*/, Mdouble surfaceEnergy /*gamma_s*/,
             Mdouble thicknessDiffusion /*deltaB*D0B*/, Mdouble activationEnergy /*QB*/, Mdouble temperature /*T*/,
             Mdouble pseudoSlidingFrictionCoefficient /*\etaPart*/)
    {
        alpha_ = alpha;
        beta_ = beta;
        atomicVolume_ = atomicVolume;
        surfaceEnergy_ = surfaceEnergy;
        thicknessDiffusion_ = thicknessDiffusion;
        activationEnergy_ = activationEnergy;
        temperature_ = temperature;
        pseudoSlidingFrictionCoefficient_ = pseudoSlidingFrictionCoefficient;
        
        Mdouble boltzmannConstant /*k_B*/ = 1.38064852e-23;
        Mdouble gasConstant /*R_g*/ = 8.314459848;
        Mdouble thicknessDiffusionVacancy /*DB*/ =
                thicknessDiffusion_ * exp(-activationEnergy_ / gasConstant / temperature_);
        std::cout << thicknessDiffusionVacancy << "|" << thicknessDiffusion_ << std::endl;
        Mdouble diffusionParameter /*DeltaB*/ =
                atomicVolume_ / boltzmannConstant / temperature_ * thicknessDiffusionVacancy;
        viscosityCoefficient_ = constants::pi / (2.0 * beta * diffusionParameter);
        adhesionCoefficient_ = alpha_ / beta_ * constants::pi * surfaceEnergy_;
        slidingFrictionCoefficient_ =
                pseudoSlidingFrictionCoefficient_ * constants::pi / (2.0 * beta * diffusionParameter);
    }
    
    Mdouble getViscosityCoefficient() const
    { return viscosityCoefficient_; }
    
    Mdouble getAdhesionCoefficient() const
    { return adhesionCoefficient_; }
    
    Mdouble getSlidingFrictionCoefficient() const
    { return slidingFrictionCoefficient_; }

private:
    
    ///\brief viscous force is adhesionCoefficient_*temperature*contactRadius^4*normalRelativeVelocity
    Mdouble alpha_;
    Mdouble beta_;
    Mdouble atomicVolume_; /*Omega*/
    Mdouble surfaceEnergy_; /*gamma_s*/
    Mdouble thicknessDiffusion_; /*deltaB*D0B*/
    Mdouble activationEnergy_; /*QB*/
    Mdouble temperature_; /*T*/
    Mdouble pseudoSlidingFrictionCoefficient_; /*etaPart*/
    
    ///\brief viscous force is viscosityCoefficient_*contactRadius^4*normalRelativeVelocity
    Mdouble viscosityCoefficient_;
    
    ///\brief adhesion force is adhesionCoefficient_*radius
    Mdouble adhesionCoefficient_;
    
    ///\brief tangential force is slidingFrictionCoefficient_*contactRadius^2*radius*tangentialRelativeVelocity
    Mdouble slidingFrictionCoefficient_;
};

#endif
