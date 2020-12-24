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


#include "LinearPlasticViscoelasticNormalSpecies.h"
#include "Interactions/NormalForceInteractions/LinearPlasticViscoelasticInteraction.h"
#include "Species/ParticleSpecies.h"

class BaseParticle;

class BaseInteractable;

LinearPlasticViscoelasticNormalSpecies::LinearPlasticViscoelasticNormalSpecies()
        : BaseNormalForce()
{
    loadingStiffness_ = 0.0;
    unloadingStiffnessMax_ = 0.0;
    cohesionStiffness_ = 0.0;
    penetrationDepthMax_ = 0.0;
    dissipation_ = 0.0;
    doConstantUnloadingStiffness_ = 0.0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticNormalSpecies::LinearPlasticViscoelasticNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
LinearPlasticViscoelasticNormalSpecies::LinearPlasticViscoelasticNormalSpecies(
        const LinearPlasticViscoelasticNormalSpecies& p)
        : BaseNormalForce(p)
{
    loadingStiffness_ = p.loadingStiffness_;
    unloadingStiffnessMax_ = p.unloadingStiffnessMax_;
    cohesionStiffness_ = p.cohesionStiffness_;
    penetrationDepthMax_ = p.penetrationDepthMax_;
    dissipation_ = p.dissipation_;
    doConstantUnloadingStiffness_ = p.doConstantUnloadingStiffness_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticNormalSpecies::LinearPlasticViscoelasticNormalSpecies(const LinearPlasticViscoelasticNormalSpecies &p) finished"<<std::endl;
#endif
}

LinearPlasticViscoelasticNormalSpecies::~LinearPlasticViscoelasticNormalSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LinearPlasticViscoelasticNormalSpecies::~LinearPlasticViscoelasticNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void LinearPlasticViscoelasticNormalSpecies::write(std::ostream& os) const
{
    os << " loadingStiffness " << loadingStiffness_;
    os << " maxUnloadingStiffness " << unloadingStiffnessMax_;
    if (doConstantUnloadingStiffness_) os << " doConstantUnloadingStiffness " << doConstantUnloadingStiffness_;
    os << " cohesionStiffness " << cohesionStiffness_;
    os << " maxPenetration " << penetrationDepthMax_;
    os << " dissipation " << dissipation_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void LinearPlasticViscoelasticNormalSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> loadingStiffness_;
    is >> dummy >> unloadingStiffnessMax_;
    helpers::readOptionalVariable<bool>(is,"doConstantUnloadingStiffness",doConstantUnloadingStiffness_);
    is >> dummy >> cohesionStiffness_;
    is >> dummy >> penetrationDepthMax_;
    is >> dummy >> dissipation_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string LinearPlasticViscoelasticNormalSpecies::getBaseName() const
{
    return "LinearPlasticViscoelastic";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void LinearPlasticViscoelasticNormalSpecies::mix(LinearPlasticViscoelasticNormalSpecies* const S,
                                                 LinearPlasticViscoelasticNormalSpecies* const T)
{
    loadingStiffness_ = BaseSpecies::average(S->getLoadingStiffness(), T->getLoadingStiffness());
    unloadingStiffnessMax_ = BaseSpecies::average(S->getUnloadingStiffnessMax(), T->getUnloadingStiffnessMax());
    cohesionStiffness_ = BaseSpecies::average(S->getCohesionStiffness(), T->getCohesionStiffness());
    penetrationDepthMax_ = BaseSpecies::average(S->getPenetrationDepthMax(), T->getPenetrationDepthMax());
    dissipation_ = BaseSpecies::average(S->getDissipation(), T->getDissipation());
    logger.assert(S->doConstantUnloadingStiffness_==T->doConstantUnloadingStiffness_,"you cannot mix species where doConstantUnloadingStiffness is not the same");
    doConstantUnloadingStiffness_ = S->doConstantUnloadingStiffness_;
}

/*!
 * \param[in] loadingStiffness      the loading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] unloadingStiffnessMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] cohesionStiffness     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void
LinearPlasticViscoelasticNormalSpecies::setPlasticParameters(Mdouble loadingStiffness, Mdouble unloadingStiffnessMax,
                                                             Mdouble cohesionStiffness, Mdouble penetrationDepthMax)
{
    if (loadingStiffness <= 0 || unloadingStiffnessMax < loadingStiffness || cohesionStiffness < 0 ||
        penetrationDepthMax < 0 || penetrationDepthMax > 1)
    {
        std::cerr << "Error: arguments of setPlasticParameters do not make sense" << std::endl;
        exit(-1);
    }
    setLoadingStiffness(loadingStiffness);
    setUnloadingStiffnessMax(unloadingStiffnessMax);
    setCohesionStiffness(cohesionStiffness);
    setPenetrationDepthMax(penetrationDepthMax);
}

/*!
 * \return the loading stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::getLoadingStiffness() const
{
    return loadingStiffness_;
}

/*!
 * \return the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::getUnloadingStiffnessMax() const
{
    return unloadingStiffnessMax_;
}

/*!
 * \return the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::getCohesionStiffness() const
{
    return cohesionStiffness_;
}

/*!
 * \return the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::getPenetrationDepthMax() const
{
    return penetrationDepthMax_;
}

/*!
 * \param[in] loadingStiffness      the loading stiffness of the linear plastic-viscoelastic normal force.
 */
void LinearPlasticViscoelasticNormalSpecies::setLoadingStiffness(Mdouble loadingStiffness)
{
    loadingStiffness_ = loadingStiffness;
}

/*!
 * \param[in] unloadingStiffnessMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
void LinearPlasticViscoelasticNormalSpecies::setUnloadingStiffnessMax(Mdouble unloadingStiffnessMax)
{
    unloadingStiffnessMax_ = unloadingStiffnessMax;
}

/*!
 * \param[in] cohesionStiffness     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
void LinearPlasticViscoelasticNormalSpecies::setCohesionStiffness(Mdouble cohesionStiffness)
{
    cohesionStiffness_ = cohesionStiffness;
}

/*!
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void LinearPlasticViscoelasticNormalSpecies::setPenetrationDepthMax(Mdouble penetrationDepthMax)
{
    penetrationDepthMax_ = penetrationDepthMax;
}

/*!
 * \details Calculates collision time for stiffest spring constant, divides by 50
 * \param[in] the optimal time step is computed to resolve a collision of two particles of this mass.
 * If constant restitution is enabled, the collision time is mass-independent.
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::computeTimeStep(Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    return 0.02 * constants::pi /
           std::sqrt(unloadingStiffnessMax_ / (.5 * mass) - mathsFunc::square(dissipation_ / mass));
}

/*!
 * \details should be non-negative
 * \param[in] the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void LinearPlasticViscoelasticNormalSpecies::setDissipation(Mdouble dissipation)
{
    logger.assert(dissipation >= 0, "setDissipation(%): dissipation should be non-negative",dissipation);
    dissipation_ = dissipation;
}

/*!
 * \param[in] a helper struct containing both the loading stiffness and the disssipation coefficient.
 */
void LinearPlasticViscoelasticNormalSpecies::setLoadingStiffnessAndDissipation(helpers::KAndDisp new_)
{
    setLoadingStiffness(new_.k);
    setDissipation(new_.disp);
}

/*!
 * \return the linear dissipation coefficient
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::getDissipation() const
{
    return dissipation_;
}

/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of equal mass m
 * \param[in] tc collision time
 * \param[in] eps restitution coefficient
 * \param[in] mass harmonic mean of particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 * If constant restitution is enabled, the collision time and restitution is mass-independent.
 */
///\todo TW: check that the masses are described correctly here (m_eff or m_p?))
void
LinearPlasticViscoelasticNormalSpecies::setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    if (eps == 0.0)
    {
        loadingStiffness_ = .5 * mass * mathsFunc::square(constants::pi / tc);
        dissipation_ = std::sqrt(2.0 * mass * loadingStiffness_);
    }
    else
    {
        dissipation_ = -mass / tc * std::log(eps);
        loadingStiffness_ =
                .5 * mass * (mathsFunc::square(constants::pi / tc) + mathsFunc::square(dissipation_ / mass));
    }
    unloadingStiffnessMax_ = loadingStiffness_;
}

/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of P
 * \param[in] stiffness stiffness
 * \param[in] eps restitution coefficient
 * \param[in] mass effective particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 * If constant restitution is enabled, the restitution coefficient is mass-independent.
 */
void LinearPlasticViscoelasticNormalSpecies::setStiffnessAndRestitutionCoefficient(Mdouble stiffness, Mdouble eps, Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    loadingStiffness_ = stiffness;
    setRestitutionCoefficient(eps, mass);
}

/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of P
 * \param[in] stiffness stiffness
 * \param[in] eps restitution coefficient
 * \param[in] mass effective particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 * If constant restitution is enabled, the restitution coefficient is mass-independent.
 */
void LinearPlasticViscoelasticNormalSpecies::setRestitutionCoefficient(Mdouble eps, Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    if (eps == 0.0) {
        dissipation_ = std::sqrt(2.0 * mass * getLoadingStiffness());
    } else {
        const Mdouble logEps = log(eps);
        dissipation_ = -std::sqrt(2.0 * mass * getLoadingStiffness()
                / (constants::sqr_pi + mathsFunc::square(logEps))) * logEps;
    }
}


/** Calculates collision time for two copies of a particle of given disp, k, mass
 * If constant restitution is enabled, the collision time is mass-independent.
 * \todo should this use unloading stiffness?
 */
Mdouble LinearPlasticViscoelasticNormalSpecies::getCollisionTime(Mdouble mass) const
{
    if (getConstantRestitution()) mass = 1;

    logger.assert(mass > 0, "getCollisionTime(%): mass should be positive",mass);

    Mdouble elasticContribution = loadingStiffness_ / (.5 * mass);
    Mdouble elastoDissipativeContribution = elasticContribution - mathsFunc::square(dissipation_ / mass);

    logger.assert(elastoDissipativeContribution > -1e-8 * elasticContribution,
            "getCollisionTime(%): values for mass, stiffness and dissipation lead to an overdamped system: reduce dissipation",mass);

    return constants::pi / std::sqrt(elastoDissipativeContribution);
}

Mdouble LinearPlasticViscoelasticNormalSpecies::getRestitutionCoefficient(Mdouble mass) const
{
    if (getConstantRestitution()) mass = 1;

    return std::exp(-dissipation_ / mass * getCollisionTime(mass));
}

Mdouble LinearPlasticViscoelasticNormalSpecies::computeBondNumberMax(Mdouble harmonicMeanRadius, Mdouble gravitationalAcceleration) const {
    if (getConstantRestitution()) {
        //harmonicMeanRadius unused
        const Mdouble plasticOverlapMax = penetrationDepthMax_;
        const Mdouble overlapMaxCohesion = plasticOverlapMax / (1 + cohesionStiffness_ / unloadingStiffnessMax_);
        const Mdouble cohesionAccelerationMax = cohesionStiffness_ * overlapMaxCohesion;
        return cohesionAccelerationMax / gravitationalAcceleration;
    } else {
        const Mdouble plasticOverlapMax = penetrationDepthMax_ * harmonicMeanRadius;
        const Mdouble overlapMaxCohesion = plasticOverlapMax / (1 + cohesionStiffness_ / unloadingStiffnessMax_);
        const Mdouble cohesionForceMax = cohesionStiffness_ * overlapMaxCohesion;
        auto species = dynamic_cast<const ParticleSpecies*>(getBaseSpecies());
        logger.assert(species,"computeBondNumberMax: species needs to be a ParticleSpecies");
        const Mdouble gravitationalForce = gravitationalAcceleration * species->getMassFromRadius(harmonicMeanRadius);
        return cohesionForceMax / gravitationalForce;
    }
}
