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


#include "HertzianBSHPViscoelasticInteraction.h"
#include "Species/NormalForceSpecies/HertzianBSHPViscoelasticNormalSpecies.h"
#include "Species/BaseSpecies.h"
#include "BaseInteractable.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
HertzianBSHPViscoelasticInteraction::HertzianBSHPViscoelasticInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                 unsigned timeStamp)
        : HertzianViscoelasticInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianBSHPViscoelasticInteraction::HertzianBSHPViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
HertzianBSHPViscoelasticInteraction::HertzianBSHPViscoelasticInteraction(const HertzianBSHPViscoelasticInteraction& p)
        : HertzianViscoelasticInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianBSHPViscoelasticInteraction::HertzianBSHPViscoelasticInteraction(const HertzianBSHPViscoelasticInteraction& p) finished"<<std::endl;
#endif
}

HertzianBSHPViscoelasticInteraction::HertzianBSHPViscoelasticInteraction() = default;

/*!
 *
 */
HertzianBSHPViscoelasticInteraction::~HertzianBSHPViscoelasticInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"HertzianBSHPViscoelasticInteraction::~HertzianBSHPViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void HertzianBSHPViscoelasticInteraction::write(std::ostream& os) const
{
    HertzianViscoelasticInteraction::write(os);
}

/*!
 * \param[in,out] is
 */
void HertzianBSHPViscoelasticInteraction::read(std::istream& is)
{
    HertzianViscoelasticInteraction::read(is);
}

/*!
 * \return std::string
 */
std::string HertzianBSHPViscoelasticInteraction::getBaseName() const
{
    return "HertzianBSHPViscoelastic";
}

/*!
 * The contact model is based on the description given in
 * https://doi.org/10.1103/PhysRevE.53.5382
 * It contains a viscoelastic damping proportional to \sqrt{\delta}
 * The force is calculated with
 * kn = 4/3 E* sqrt(R*)
 * fn = kn (delta^{3/2} + 3/2 kn A delta^{1/2} \dot{\delta})
 *
 * https://doi.org/10.1103/PhysRevE.84.021302 provides a way for calculating
 * the damping constant A from the coefficient of restitution, a reference velocity
 * and the material parameters.
 */
void HertzianBSHPViscoelasticInteraction::computeNormalForce()
{
    // This function is called for all particles within interactionRadius distance.
    
    // This has to be outside the loop because it is needed for the other forces
    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(
            getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));
    
    if (getOverlap() > 0) //if contact forces
    {
        const HertzianBSHPViscoelasticNormalSpecies* species = getSpecies();
        
        Mdouble stiffness = 4. / 3. * species->getEffectiveElasticModulus() * std::sqrt(getEffectiveRadius() * getOverlap());
        
        //calculating the current normal force
        //dissipation is computed such that the restitution is constant
        Mdouble dissipationCoefficient = 3. / 2. * species->getDissipation() * stiffness;
        Mdouble normalForce = stiffness * getOverlap() - dissipationCoefficient * getNormalRelativeVelocity();
        
        // Ensuring positive normal force
        normalForce = std::max(0.0, normalForce);
        
        //setting the normal force parameter in the base interaction class so that it can be accessed
        //by other classes...        
        setAbsoluteNormalForce(std::abs(normalForce)); //used for further force calculations;
        setForce(getNormal() * normalForce);
        ///\todo check for superquadrics
        setTorque(Vec3D(0.0, 0.0, 0.0));
    }
    else
    {
        setAbsoluteNormalForce(0.0);
        setForce(Vec3D(0.0, 0.0, 0.0));
        setTorque(Vec3D(0.0, 0.0, 0.0));
    }
}

/*!
 * \return const HertzianViscoelasticNormalSpecies*
 */
const HertzianBSHPViscoelasticNormalSpecies* HertzianBSHPViscoelasticInteraction::getSpecies() const
{
    return static_cast<const HertzianBSHPViscoelasticNormalSpecies*>(getBaseSpecies()->getNormalForce());
}
