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


#include "SinterNormalSpecies.h"
#include "Interactions/NormalForceInteractions/SinterInteraction.h"
#include "Logger.h"
#include "Species/BaseSpecies.h"

class BaseParticle;

class BaseInteractable;

SinterNormalSpecies::SinterNormalSpecies()
        : BaseNormalForce()
{
    sinterType_ = SINTERTYPE::PARHAMI_MCKEEPING;
    loadingStiffness_ = 0.0;
    unloadingStiffnessMax_ = 0.0;
    cohesionStiffness_ = 0.0;
    penetrationDepthMax_ = 0.0;
    dissipation_ = 0.0;
    sinterAdhesion_ = 0.0;
    sinterRate_ = 0.0;
    inverseSinterViscosity_ = 0.0;
    complianceZero_ = 0.0;
    surfTension_ = 0.0;
    constantC1_ = 0.0;
    separationDis_ = 0.0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SinterNormalSpecies::SinterNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
SinterNormalSpecies::SinterNormalSpecies(const SinterNormalSpecies& p)
        : BaseNormalForce(p)
{
    sinterType_ = p.sinterType_;
    loadingStiffness_ = p.loadingStiffness_;
    unloadingStiffnessMax_ = p.unloadingStiffnessMax_;
    cohesionStiffness_ = p.cohesionStiffness_;
    penetrationDepthMax_ = p.penetrationDepthMax_;
    dissipation_ = p.dissipation_;
    sinterAdhesion_ = p.sinterAdhesion_;
    inverseSinterViscosity_ = p.inverseSinterViscosity_;
    sinterRate_ = p.sinterRate_;
    complianceZero_ = p.complianceZero_;
    surfTension_ = p.surfTension_;
    constantC1_ = p.constantC1_;
    separationDis_ = p.separationDis_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SinterNormalSpecies::SinterNormalSpecies(const SinterNormalSpecies &p) finished"<<std::endl;
#endif
}

SinterNormalSpecies::~SinterNormalSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SinterNormalSpecies::~SinterNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void SinterNormalSpecies::write(std::ostream& os) const
{
    os << " loadingStiffness " << loadingStiffness_;
    os << " maxUnloadingStiffness " << unloadingStiffnessMax_;
    os << " cohesionStiffness " << cohesionStiffness_;
    os << " maxPenetration " << penetrationDepthMax_;
    os << " dissipation " << dissipation_;
    os << " sinterAdhesion " << sinterAdhesion_;
    os << " inverseSinterViscosity " << inverseSinterViscosity_;
    os << " sinterRate " << sinterRate_;
    os << " complianceZero " << complianceZero_;
    os << " surfTension " << surfTension_;
    os << " constantC1 " << constantC1_;
    os << " separationDis " << separationDis_;
    os << " sinterType " << (unsigned) sinterType_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void SinterNormalSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> loadingStiffness_;
    is >> dummy >> unloadingStiffnessMax_;
    is >> dummy >> cohesionStiffness_;
    is >> dummy >> penetrationDepthMax_;
    is >> dummy >> dissipation_;
    is >> dummy >> sinterAdhesion_;
    is >> dummy >> inverseSinterViscosity_;
    is >> dummy >> sinterRate_;
    is >> dummy >> complianceZero_;
    is >> dummy >> surfTension_;
    is >> dummy >> constantC1_;
    is >> dummy >> separationDis_;
    unsigned type;
    is >> dummy >> type;
    sinterType_ = (SINTERTYPE) type;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string SinterNormalSpecies::getBaseName() const
{
    return "Sinter";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void SinterNormalSpecies::mix(SinterNormalSpecies* const S, SinterNormalSpecies* const T)
{
    loadingStiffness_ = BaseSpecies::average(S->getLoadingStiffness(), T->getLoadingStiffness());
    unloadingStiffnessMax_ = BaseSpecies::average(S->getUnloadingStiffnessMax(), T->getUnloadingStiffnessMax());
    cohesionStiffness_ = BaseSpecies::average(S->getCohesionStiffness(), T->getCohesionStiffness());
    penetrationDepthMax_ = BaseSpecies::average(S->getPenetrationDepthMax(), T->getPenetrationDepthMax());
    dissipation_ = BaseSpecies::average(S->getDissipation(), T->getDissipation());
    inverseSinterViscosity_ = BaseSpecies::average(S->getInverseSinterViscosity(), T->getInverseSinterViscosity());
    sinterAdhesion_ = BaseSpecies::average(S->getSinterAdhesion(), T->getSinterAdhesion());
    complianceZero_ = BaseSpecies::average(S->getComplianceZero(), T->getComplianceZero());
    surfTension_ = BaseSpecies::average(S->getSurfTension(), T->getSurfTension());
    constantC1_= BaseSpecies::average(S->getConstantC1(), T->getConstantC1());
    separationDis_ = BaseSpecies::average(S->getSeparationDis(), T->getSeparationDis());

}

/*!
 * \param[in] loadingStiffness      the loading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] unloadingStiffnessMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] cohesionStiffness     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setPlasticParameters(Mdouble loadingStiffness, Mdouble unloadingStiffnessMax,
                                               Mdouble cohesionStiffness, Mdouble penetrationDepthMax)
{
    if (loadingStiffness <= 0 || unloadingStiffnessMax < loadingStiffness || cohesionStiffness < 0 ||
        penetrationDepthMax < 0)
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
Mdouble SinterNormalSpecies::getLoadingStiffness() const
{
    return loadingStiffness_;
}

/*!
 * \return the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble SinterNormalSpecies::getUnloadingStiffnessMax() const
{
    return unloadingStiffnessMax_;
}

/*!
 * \return the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble SinterNormalSpecies::getCohesionStiffness() const
{
    return cohesionStiffness_;
}

/*!
 * \return the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
Mdouble SinterNormalSpecies::getPenetrationDepthMax() const
{
    return penetrationDepthMax_;
}

/*!
 * \param[in] loadingStiffness      the loading stiffness of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setLoadingStiffness(Mdouble loadingStiffness)
{
    loadingStiffness_ = loadingStiffness;
}

/*!
 * \param[in] unloadingStiffnessMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setUnloadingStiffnessMax(Mdouble unloadingStiffnessMax)
{
    unloadingStiffnessMax_ = unloadingStiffnessMax;
}

/*!
 * \param[in] cohesionStiffness     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setCohesionStiffness(Mdouble cohesionStiffness)
{
    cohesionStiffness_ = cohesionStiffness;
}

/*!
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setPenetrationDepthMax(Mdouble penetrationDepthMax)
{
    penetrationDepthMax_ = penetrationDepthMax;
}

/*!
 * \details Calculates collision time for stiffest spring constant, divides by 50
 * \param[in] mass the optimal time step is computed to resolve a collision of two particles of this mass.
 */
Mdouble SinterNormalSpecies::computeTimeStep(Mdouble mass)
{
    if (unloadingStiffnessMax_ / (.5 * mass) < mathsFunc::square(dissipation_ / mass))
    {
        std::cerr << "Dissipation too high; max. allowed " << sqrt(2.0 * unloadingStiffnessMax_ * mass) << std::endl;
        exit(-1);
    }
    return 0.02 * constants::pi /
           std::sqrt(unloadingStiffnessMax_ / (.5 * mass) - mathsFunc::square(dissipation_ / mass));
}

/*!
 * \details should be non-negative
 * \param[in] dissipation the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setDissipation(Mdouble dissipation)
{
    if (dissipation < 0)
    {
        logger(ERROR, "setDissipation(%)", dissipation);
        exit(-1);
    }
    dissipation_ = dissipation;
}

/*!
 * \details should be non-negative
 * \param[in] sinterAdhesion_
 */
void SinterNormalSpecies::setSinterAdhesion(Mdouble sinterAdhesion)
{
    if (sinterAdhesion < 0)
    {
        logger(ERROR, "setSinterAdhesion(%)", sinterAdhesion);
        exit(-1);
    }
    sinterAdhesion_ = sinterAdhesion;
}

/*!
 * \details should be non-negative
 * \param[in] sinterAdhesion_
 */
void SinterNormalSpecies::setSinterRate(Mdouble sinterRate)
{
    if (sinterRate < 0)
    {
        logger(ERROR, "setSinterRate(%)", sinterRate);
    }
    sinterRate_ = sinterRate;
}

void SinterNormalSpecies::setComplianceZero(Mdouble complianceZero)
{
    if (complianceZero < 0)
    {
        logger(ERROR, "complianceZero(%)", complianceZero);
    }
    complianceZero_ = complianceZero;
}

void SinterNormalSpecies::setSurfTension(Mdouble surfTension)
{
    if (surfTension < 0)
    {
        logger(ERROR, "SurfTension(%)", surfTension);
    }
    surfTension_ = surfTension;
}

void SinterNormalSpecies::setConstantC1(Mdouble complianceOne)
{
    if (complianceOne < 0)
    {
        logger(ERROR, "ComplianceOne(%)", complianceOne);
    }
    constantC1_ = complianceOne;
}

void SinterNormalSpecies::setSeparationDis(Mdouble separationDis)
{
    if (separationDis < 0)
    {
        logger(ERROR, "SeparationDis(%)", separationDis);
    }
    separationDis_ = separationDis;
}




/*!
 * \details should be non-negative
 * \param[in] sinterAdhesion_
 */
void SinterNormalSpecies::setSinterType(SINTERTYPE sinterType)
{
    sinterType_ = sinterType;
    if (sinterType_ == SINTERTYPE::CONSTANT_RATE)
        logger(INFO, "Sintertype set to CONSTANT_RATE");
    else if (sinterType_ == SINTERTYPE::PARHAMI_MCKEEPING)
        logger(INFO, "Sintertype set to PARHAMI_MCKEEPING");
    else if (sinterType_ == SINTERTYPE::TEMPERATURE_DEPENDENT_FRENKEL)
        logger(INFO, "Sintertype set to TEMPERATURE_DEPENDENT_FRENKEL");
    else if (sinterType_ == SINTERTYPE::REGIME_SINTERING)
        logger(INFO, "Sintertype set to REGIME_SINTERING");
    else
        logger(ERROR, "Sintertype not understood");
}

/*!
 * \details should be non-negative
 * \param[in] the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setInverseSinterViscosity(Mdouble inverseSinterViscosity)
{
    //assertOrDie(inverseSinterViscosity < 0, "inverseSinterViscosity should be non-negative!");
    if (inverseSinterViscosity < 0)
    {
        logger(ERROR, "setInverseSinterViscosity(%)", inverseSinterViscosity);
    }
    setSinterType(SINTERTYPE::PARHAMI_MCKEEPING);
    inverseSinterViscosity_ = inverseSinterViscosity;
}

void SinterNormalSpecies::setSinterForceAndTime(Mdouble adhesionForce, Mdouble sinterTime, Mdouble radius)
{
    sinterType_ = SINTERTYPE::PARHAMI_MCKEEPING;
    ///\todo fix
    setSinterAdhesion(adhesionForce / radius);
    setInverseSinterViscosity(1.0 / (93.75 * adhesionForce * sinterTime / std::pow(radius, 5)));
    logger(INFO, "Set sintering parameters: adhesion force f_a=%*r,  sinter viscosity is nu = contactRadius^4/%",
           getSinterAdhesion(), getInverseSinterViscosity());
}

/*!
 * \details should be non-negative
 * \param[in] the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void SinterNormalSpecies::setParhamiMcKeeping
        (Mdouble alpha, Mdouble beta, Mdouble atomicVolume /*Omega*/, Mdouble surfaceEnergy /*gamma_s*/,
         Mdouble thicknessDiffusion /*deltaB*D0B*/, Mdouble activationEnergy /*QB*/, Mdouble temperature /*T*/)
{
    setSinterType(SINTERTYPE::PARHAMI_MCKEEPING);
    const Mdouble boltzmannConstant /*k_B*/ = 1.38064852e-23;
    const Mdouble gasConstant /*R_g*/ = 8.314459848;
    const Mdouble thicknessDiffusionVacancy /*DB*/ =
            thicknessDiffusion * exp(-activationEnergy / gasConstant / temperature);
    const Mdouble diffusionParameter /*DeltaB*/ =
            atomicVolume / boltzmannConstant / temperature * thicknessDiffusionVacancy;
    setInverseSinterViscosity(constants::pi / (2.0 * beta * diffusionParameter));
    setSinterAdhesion(alpha / beta * constants::pi * surfaceEnergy);
    logger(INFO,
           "Set sintering parameters: the adhesion force f_a=%*r, the rate of change of the plastic overlap is d(delta0)/dt=f_ep/(%*a^4)",
           getSinterAdhesion(),
           getInverseSinterViscosity());
}

/*!
 * \param[in] new_ a helper struct containing both the loading stiffness and the disssipation coefficient.
 */
void SinterNormalSpecies::setLoadingStiffnessAndDissipation(helpers::KAndDisp new_)
{
    setLoadingStiffness(new_.k);
    setDissipation(new_.disp);
}

/*!
 * \return the linear dissipation coefficient
 */
Mdouble SinterNormalSpecies::getDissipation() const
{
    return dissipation_;
}

/*!
 * \return value of sinterAdhesion_
 */
Mdouble SinterNormalSpecies::getSinterAdhesion() const
{
    return sinterAdhesion_;
}

/*!
 * \return value of inverseSinterViscosity_
 */
Mdouble SinterNormalSpecies::getInverseSinterViscosity() const
{
    return inverseSinterViscosity_;
}

/*!
 * \return value of sinterAdhesion_
 */
Mdouble SinterNormalSpecies::getSinterRate() const
{
    return sinterRate_;
}

Mdouble SinterNormalSpecies::getComplianceZero() const
{
    return complianceZero_;
}

Mdouble SinterNormalSpecies::getSurfTension() const
{
    return surfTension_;
}

Mdouble SinterNormalSpecies::getConstantC1() const
{
    return constantC1_;
}

Mdouble SinterNormalSpecies::getSeparationDis() const
{
    return separationDis_;
}

SINTERTYPE SinterNormalSpecies::getSinterType() const
{
    return sinterType_;
}

double SinterNormalSpecies::getTemperatureDependentSinterRate(double temperature) const
{
    return temperatureDependentSinterRate_(temperature);
}

std::function<double(double temperature)> SinterNormalSpecies::getTemperatureDependentSinterRate() const
{
    return temperatureDependentSinterRate_;
}

void SinterNormalSpecies::setTemperatureDependentSinterRate(
        std::function<double(double temperature)> temperatureDependentSinterRate)
{
    temperatureDependentSinterRate_ = temperatureDependentSinterRate;
}


/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of equal mass m
 * \param[in] tc collision time
 * \param[in] eps restitution coefficient
 * \param[in] mass effective particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 */
///\todo TW: check that the masses are described correctly here (m_eff or m_p?))
void SinterNormalSpecies::setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
{
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
 */
void SinterNormalSpecies::setStiffnessAndRestitutionCoefficient(Mdouble stiffness, Mdouble eps, Mdouble mass)
{
    loadingStiffness_ = stiffness;
    if (eps == 0.0)
    {
        dissipation_ = std::sqrt(2.0 * mass * stiffness);
    }
    else
    {
        dissipation_ =
                -std::sqrt(2.0 * mass * stiffness / (constants::sqr_pi + mathsFunc::square(log(eps)))) * log(eps);
    }
    unloadingStiffnessMax_ = loadingStiffness_;
}

///Calculates collision time for two copies of a particle of given disp, k, mass
Mdouble SinterNormalSpecies::getCollisionTime(Mdouble mass)
{
    if (mass <= 0)
    {
        std::cerr << "Error in getCollisionTime(" << mass
                  << ") mass is not set or has an unexpected value, (getCollisionTime(" << mass << "))" << std::endl;
        exit(-1);
    }
    if (loadingStiffness_ <= 0)
    {
        std::cerr << "Error in getCollisionTime(" << mass << ") stiffness=" << loadingStiffness_
                  << " is not set or has an unexpected value, (getCollisionTime(" << mass << "), with stiffness="
                  << loadingStiffness_ << ")" << std::endl;
        exit(-1);
    }
    if (dissipation_ < 0)
    {
        std::cerr << "Error in getCollisionTime(" << mass << ") dissipation=" << dissipation_
                  << " is not set or has an unexpected value, (getCollisionTime(" << mass << "), with dissipation="
                  << dissipation_ << ")" << std::endl;
        exit(-1);
    }
    Mdouble tosqrt = loadingStiffness_ / (.5 * mass) - mathsFunc::square(dissipation_ / mass);
    if (tosqrt <= -1e-8 * loadingStiffness_ / (.5 * mass))
    {
        std::cerr << "Error in getCollisionTime(" << mass
                  << ") values for mass, stiffness and dissipation would lead to an overdamped system, (getCollisionTime("
                  << mass << "), with stiffness=" << loadingStiffness_ << " and dissipation=" << dissipation_ << ")"
                  << std::endl;
        exit(-1);
    }
    return constants::pi / std::sqrt(tosqrt);
}
