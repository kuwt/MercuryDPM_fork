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


#include "HertzianViscoelasticInteraction.h"
#include "Species/NormalForceSpecies/HertzianViscoelasticNormalSpecies.h"
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
HertzianViscoelasticInteraction::HertzianViscoelasticInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                 unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianViscoelasticInteraction::HertzianViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
HertzianViscoelasticInteraction::HertzianViscoelasticInteraction(const HertzianViscoelasticInteraction& p)
        : BaseInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianViscoelasticInteraction::HertzianViscoelasticInteraction(const HertzianViscoelasticInteraction& p) finished"<<std::endl;
#endif
}

HertzianViscoelasticInteraction::HertzianViscoelasticInteraction() = default;

/*!
 *
 */
HertzianViscoelasticInteraction::~HertzianViscoelasticInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"HertzianViscoelasticInteraction::~HertzianViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void HertzianViscoelasticInteraction::write(std::ostream& os) const
{
    BaseInteraction::write(os);
}

/*!
 * \param[in,out] is
 */
void HertzianViscoelasticInteraction::read(std::istream& is)
{
    BaseInteraction::read(is);
}

/*!
 * \return std::string
 */
std::string HertzianViscoelasticInteraction::getBaseName() const
{
    return "HertzianViscoelastic";
}

/*!
 * The contact model is based on the description given in
 * http://people.ds.cam.ac.uk/jae1001/CUS/research/pfizer/Antypov_Elliott_EPL_2011.pdf
 * (which is the same as in Yade https://answers.launchpad.net/yade/+question/235934)
 * kn = 4/3 E* sqrt(R/2 delta) = sqrt(2)/3 kn
 * fn = kn delta + gamma sqrt(m/2 kn) vn =
 *
 * Note, the constants are slightly different in the C Thornton model, otherwise it's the same:
 * www.cfd.com.au/cfd_conf12/PDFs/175CUM.pdf
 */
void HertzianViscoelasticInteraction::computeNormalForce()
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
        const HertzianViscoelasticNormalSpecies* species = getSpecies();
        
        Mdouble stiffness = 4. / 3. * species->getEffectiveElasticModulus() * std::sqrt(getEffectiveRadius() * getOverlap());
        
        //calculating the current normal force
        //dissipation is computed such that the restitution is constant
        Mdouble dissipationCoefficient = species->getDissipation() * sqrt(getEffectiveMass() * stiffness);
        Mdouble normalForce = stiffness * getOverlap() - dissipationCoefficient * getNormalRelativeVelocity();
        
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
 * kn = 4/3 E* sqrt(R/2 delta)
 * fn = kn delta
 * -> E = int fn ddelta = 2/5 4/3 E sqrt(R/2 delta) delta^2
 */
Mdouble HertzianViscoelasticInteraction::getElasticEnergy() const
{
    if (getOverlap() > 0)
    {
        return 8. / 15. * getSpecies()->getEffectiveElasticModulus() * std::sqrt(getEffectiveRadius() * getOverlap()) *
               mathsFunc::square(getOverlap());
    }
    else
    {
        return 0.0;
    }
}

/*!
 * \return const HertzianViscoelasticNormalSpecies*
 */
const HertzianViscoelasticNormalSpecies* HertzianViscoelasticInteraction::getSpecies() const
{
    return static_cast<const HertzianViscoelasticNormalSpecies*>(getBaseSpecies()->getNormalForce());
}

/*!
 * Computes elastic-adhesive energy at zero overlap, assuming the energy is zero at the equilibrium point
 * mod = 4/3 E* sqrt(R/2)
 * fn = mod delta^3/2 - fadh
 * -> E = int_deltaEq^0 fn ddelta = [2/5 mod delta^5/2 - fadh * delta]_deltaEq^0 = [2/5 fn delta - 3/5 fadh delta]_deltaEq^0 = 3/5 fadh deltaEq
 * \todo TW consider renaming to getElasticAdhesiveEnergy or getElasticAdhesiveEnergyRelativeToEquilibrium
 */
Mdouble HertzianViscoelasticInteraction::getElasticEnergyAtEquilibrium(Mdouble adhesiveForce) const
{
    const HertzianViscoelasticNormalSpecies* species = getSpecies();
    const Mdouble modulus = 4. / 3. * species->getEffectiveElasticModulus() * std::sqrt(getEffectiveRadius());
    const Mdouble equilibriumOverlap = std::cbrt(mathsFunc::square(adhesiveForce / modulus));
    return 0.6 * adhesiveForce * equilibriumOverlap;//why not 0.4?
}
