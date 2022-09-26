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

#include "HertzianBSHPViscoelasticNormalSpecies.h"
#include "Interactions/BaseInteraction.h"
#include <cmath>
#include <Species/FrictionForceSpecies/MindlinSpecies.h>
#include "Interactions/NormalForceInteractions/HertzianBSHPViscoelasticInteraction.h"
#include "Logger.h"
#include "Species/BaseSpecies.h"

class BaseParticle;

class BaseInteractable;

HertzianBSHPViscoelasticNormalSpecies::HertzianBSHPViscoelasticNormalSpecies()
        : HertzianViscoelasticNormalSpecies()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianBSHPViscoelasticNormalSpecies::HertzianBSHPViscoelasticNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
HertzianBSHPViscoelasticNormalSpecies::HertzianBSHPViscoelasticNormalSpecies(const HertzianBSHPViscoelasticNormalSpecies& p)
        : HertzianViscoelasticNormalSpecies(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianBSHPViscoelasticNormalSpecies::HertzianBSHPViscoelasticNormalSpecies(const HertzianBSHPViscoelasticNormalSpecies &p) finished"<<std::endl;
#endif
}

HertzianBSHPViscoelasticNormalSpecies::~HertzianBSHPViscoelasticNormalSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"HertzianBSHPViscoelasticNormalSpecies::~HertzianBSHPViscoelasticNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string HertzianBSHPViscoelasticNormalSpecies::getBaseName() const
{
    return "HertzianBSHPViscoelastic";
}

///Allows the spring constant to be changed
void HertzianBSHPViscoelasticNormalSpecies::setEffectiveElasticModulusAndRestitutionCoefficient(Mdouble elasticModulus, Mdouble rest)
{

    logger(ERROR,
           "[HertzianBSHPViscoelasticNormalSpecies::setEffectiveElasticModulusAndRestitutionCoefficient] This function is not supported for the chosen contact model.");
}

// /**
//  * \param[in] relativeVelocity input the maximum relative velocity in your system to get the mininimum collision time
//  * \param[in] particleDiameter input the minimum particle diameter in your system to get the mininimum collision time
//  * \param[in] particleDensity input the minimum particle density in your system to get the mininimum collision time
//  */
// Mdouble HertzianBSHPViscoelasticNormalSpecies::getCollisionTime(Mdouble particleDiameter, Mdouble particleDensity,
//                                                             Mdouble relativeVelocity) const
// {
//     // Here is a very nice paper describing contact modelling
//     // http://people.ds.cam.ac.uk/jae1001/CUS/research/pfizer/Antypov_Elliott_EPL_2011.pdf
//     //caution: this function assumes the contact is elastic (no dissipation)
//     //Mdouble omega0 = 2.0/constants::sqrt_pi*sqrt(getEffectiveElasticModulus()/particleDensity)/relativeVelocity;
//     //Mdouble omega1 = sqrt(omega0*omega0-getDissipation()*getDissipation());
//     return 2.214 * pow(mathsFunc::square(particleDensity / getEffectiveElasticModulus()) / relativeVelocity, 0.2) *
//            particleDiameter;
// }
