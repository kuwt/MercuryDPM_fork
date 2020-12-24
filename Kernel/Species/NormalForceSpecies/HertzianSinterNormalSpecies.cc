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


#include "HertzianSinterNormalSpecies.h"
#include "Interactions/NormalForceInteractions/HertzianSinterInteraction.h"
#include "BaseHandler.h"
#include "Species/ParticleSpecies.h"
//#include <cassert>

class BaseParticle;

class BaseInteractable;

HertzianSinterNormalSpecies::HertzianSinterNormalSpecies()
        : BaseNormalForce()
{
    loadingModulus_ = 0.0;
    unloadingModulusMax_ = 0.0;
    cohesionModulus_ = 0.0;
    penetrationDepthMax_ = 0.0;
    dissipation_ = 0.0;
    sinterRate_ = 0.0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianSinterNormalSpecies::HertzianSinterNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
HertzianSinterNormalSpecies::HertzianSinterNormalSpecies(const HertzianSinterNormalSpecies& p)
        : BaseNormalForce(p)
{
    loadingModulus_ = p.loadingModulus_;
    unloadingModulusMax_ = p.unloadingModulusMax_;
    cohesionModulus_ = p.cohesionModulus_;
    penetrationDepthMax_ = p.penetrationDepthMax_;
    dissipation_ = p.dissipation_;
    sinterRate_ = p.sinterRate_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"HertzianSinterNormalSpecies::HertzianSinterNormalSpecies(const HertzianSinterNormalSpecies &p) finished"<<std::endl;
#endif
}

HertzianSinterNormalSpecies::~HertzianSinterNormalSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"HertzianSinterNormalSpecies::~HertzianSinterNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void HertzianSinterNormalSpecies::write(std::ostream& os) const
{
    os << " loadingModulus " << loadingModulus_;
    os << " maxUnloadingModulus " << unloadingModulusMax_;
    os << " cohesionModulus " << cohesionModulus_;
    os << " maxPenetration " << penetrationDepthMax_;
    os << " dissipation " << dissipation_;
    os << " sinterRate " << sinterRate_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void HertzianSinterNormalSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> loadingModulus_;
    is >> dummy >> unloadingModulusMax_;
    is >> dummy >> cohesionModulus_;
    is >> dummy >> penetrationDepthMax_;
    is >> dummy >> dissipation_;
    is >> dummy >> sinterRate_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string HertzianSinterNormalSpecies::getBaseName() const
{
    return "HertzianSinter";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void HertzianSinterNormalSpecies::mix(HertzianSinterNormalSpecies* const S, HertzianSinterNormalSpecies* const T)
{
    loadingModulus_ = BaseSpecies::average(S->getLoadingModulus(), T->getLoadingModulus());
    unloadingModulusMax_ = BaseSpecies::average(S->getUnloadingModulusMax(), T->getUnloadingModulusMax());
    cohesionModulus_ = BaseSpecies::average(S->getCohesionModulus(), T->getCohesionModulus());
    penetrationDepthMax_ = BaseSpecies::average(S->getPenetrationDepthMax(), T->getPenetrationDepthMax());
    dissipation_ = BaseSpecies::average(S->getDissipation(), T->getDissipation());
    sinterRate_ = BaseSpecies::average(S->getSinterRate(), T->getSinterRate());
}

/*!
 * \param[in] loadingModulus      the loading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] unloadingModulusMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] cohesionModulus     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setPlasticParameters(Mdouble loadingModulus, Mdouble unloadingModulusMax,
                                                       Mdouble cohesionModulus, Mdouble penetrationDepthMax)
{
    if (loadingModulus <= 0 || unloadingModulusMax <= 1.000001 * (loadingModulus + cohesionModulus) ||
        cohesionModulus < 0 || penetrationDepthMax < 0)
    {
        std::cerr << "Error: arguments of setPlasticParameters do not make sense" << std::endl;
        exit(-1);
    }
    setLoadingModulus(loadingModulus);
    setUnloadingModulusMax(unloadingModulusMax);
    setCohesionModulus(cohesionModulus);
    setPenetrationDepthMax(penetrationDepthMax);
}

/*!
 * \return the loading stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble HertzianSinterNormalSpecies::getLoadingModulus() const
{
    return loadingModulus_;
}

/*!
 * \return the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble HertzianSinterNormalSpecies::getUnloadingModulusMax() const
{
    return unloadingModulusMax_;
}

/*!
 * \return the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble HertzianSinterNormalSpecies::getCohesionModulus() const
{
    return cohesionModulus_;
}

/*!
 * \return the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
Mdouble HertzianSinterNormalSpecies::getPenetrationDepthMax() const
{
    return penetrationDepthMax_;
}

/*!
 * \param[in] loadingModulus      the loading stiffness of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setLoadingModulus(Mdouble loadingModulus)
{
    loadingModulus_ = loadingModulus;
}

/*!
 * \param[in] unloadingModulusMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setUnloadingModulusMax(Mdouble unloadingModulusMax)
{
    unloadingModulusMax_ = unloadingModulusMax;
}

/*!
 * \param[in] cohesionModulus     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setCohesionModulus(Mdouble cohesionModulus)
{
    cohesionModulus_ = cohesionModulus;
}

/*!
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setPenetrationDepthMax(Mdouble penetrationDepthMax)
{
    penetrationDepthMax_ = penetrationDepthMax;
}

/*!
 * \details Calculates collision time for stiffest spring constant, divides by 50
 * \param[in] the optimal time step is computed to resolve a collision of two particles of this mass.
 */
Mdouble HertzianSinterNormalSpecies::computeTimeStep(Mdouble mass)
{
//    if (stiffnessMax / (.5 * mass) < mathsFunc::square(dissipation_ /mass)) {
//        std::cerr << "Dissipation too high; max. allowed " << sqrt(2.0 * stiffnessMax * mass) << std::endl;
//        return 0.02 * constants::pi / std::sqrt(2.0*stiffnessMax / mass);
//    } else {
    std::cerr << "Warning: Dissipation is not taken into account when computing the time step" << std::endl;
    ParticleSpecies* p = dynamic_cast<ParticleSpecies*>(getBaseSpecies());
    logger.assert(p,"Empty particle handler");
    Mdouble radius = cbrt(mass * 3. / (4. * constants::pi * p->getDensity()));
    return 0.02 * constants::pi / std::sqrt(2.0 * getUnloadingModulusMax() * getPenetrationDepthMax() * radius / mass);
}

/*!
 * \details should be non-negative
 * \param[in] the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setDissipation(Mdouble dissipation)
{
    if (dissipation >= 0)
    {
        dissipation_ = dissipation;
    }
    else
    {
        std::cerr << "Error in setDissipation(" << dissipation << ")" << std::endl;
        exit(-1);
    }
}

/*!
 * \details should be non-negative
 * \param[in] the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void HertzianSinterNormalSpecies::setSinterRate(Mdouble sinterRate)
{
    if (sinterRate >= 0)
    {
        sinterRate_ = sinterRate;
    }
    else
    {
        std::cerr << "Error in setSinterRate(" << sinterRate << ")" << std::endl;
        exit(-1);
    }
}

/*!
 * \return the linear dissipation coefficient
 */
Mdouble HertzianSinterNormalSpecies::getDissipation() const
{
    return dissipation_;
}

/*!
 * \return the linear dissipation coefficient
 */
Mdouble HertzianSinterNormalSpecies::getSinterRate() const
{
    return sinterRate_;
}
