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


#include "ParhamiMcMeekingSinterInteraction.h"
#include "Species/AdhesiveForceSpecies/ParhamiMcMeekingSinterSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
ParhamiMcMeekingSinterInteraction::ParhamiMcMeekingSinterInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                     unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ParhamiMcMeekingSinterInteraction::ParhamiMcMeekingSinterInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
ParhamiMcMeekingSinterInteraction::ParhamiMcMeekingSinterInteraction(const ParhamiMcMeekingSinterInteraction& p)
        : BaseInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ParhamiMcMeekingSinterInteraction::ParhamiMcMeekingSinterInteraction(const ParhamiMcMeekingSinterInteraction &p finished"<<std::endl;
#endif
}

ParhamiMcMeekingSinterInteraction::ParhamiMcMeekingSinterInteraction()
{
#ifdef MERCURY_USE_MPI
    logger(FATAL,"ParhamiMcMeekingSinterInteractions are currently not implemented in parallel MercuryDPM");
#endif
    
}

/*!
 *
 */
ParhamiMcMeekingSinterInteraction::~ParhamiMcMeekingSinterInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"ParhamiMcMeekingSinterInteraction::~ParhamiMcMeekingSinterInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] os
 */
void ParhamiMcMeekingSinterInteraction::write(std::ostream& os  UNUSED) const
{}

/*!
 * \param[in] is
 */
void ParhamiMcMeekingSinterInteraction::read(std::istream& is  UNUSED)
{}

/*!
 *
 */
void ParhamiMcMeekingSinterInteraction::computeAdhesionForce()
{
    //std::cout << "ParhamiMcMeekingSinterInteraction::computeAdhesionForce" << std::endl;
    const SpeciesType* species = getSpecies();
    Mdouble effectiveDiameter = 2.0 * getEffectiveRadius(); // effectiveRadius = (r1*r2)/(r1+r2)
    Mdouble contactRadiusSquared = 2.0 * effectiveDiameter * getOverlap();
    Vec3D tangentialRelativeVelocity = getRelativeVelocity() - getNormal() * getNormalRelativeVelocity();
    
    //viscous force is viscosityCoefficient_*contactRadius^4*normalRelativeVelocity
    //adhesion force is adhesionCoefficient_*radius
    Mdouble normalForce = -species->getAdhesionCoefficient() * effectiveDiameter
                          - species->getViscosityCoefficient() * contactRadiusSquared * contactRadiusSquared *
                            getNormalRelativeVelocity();
    //tangential force is slidingFrictionCoefficient_*contactRadius^2*radius*tangentialRelativeVelocity
    Vec3D tangentialForce = -species->getSlidingFrictionCoefficient() * contactRadiusSquared * effectiveDiameter *
                            tangentialRelativeVelocity;
    //std::cout << "P" << species->getAdhesionCoefficient() << species->getViscosityCoefficient() << species->getSlidingFrictionCoefficient() << std::endl;
    Mdouble attractiveForce = -species->getAdhesionCoefficient() * effectiveDiameter;
    //std::cout << effectiveDiameter << "Fs=" << species->getAdhesionCoefficient() * effectiveDiameter
    //<< " Fv/del/del/del=" << 4.0*species->getViscosityCoefficient()*effectiveDiameter*effectiveDiameter << std::endl;
    //<< species->getViscosityCoefficient() << " d" <<contactRadiusSquared<< " d" <<getNormalRelativeVelocity() << std::endl;
    addForce(getNormal() * normalForce + tangentialForce);
}

/*!
 * \todo TW 
 * \return Mdouble
 */
Mdouble ParhamiMcMeekingSinterInteraction::getElasticEnergy() const
{
    return 0.0;
}

/*!
 * \return a constant pointer to an instance of this class.
 */
const ParhamiMcMeekingSinterInteraction::SpeciesType* ParhamiMcMeekingSinterInteraction::getSpecies() const
{
    return static_cast<const SpeciesType*> (getBaseSpecies()->getAdhesiveForce()); //downcast
}

/*!
 * \return std::string
 */
std::string ParhamiMcMeekingSinterInteraction::getBaseName() const
{
    return "ParhamiMcMeekingSinter";
}
