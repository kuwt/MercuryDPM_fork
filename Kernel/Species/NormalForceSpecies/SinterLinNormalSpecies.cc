//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


#include "SinterLinNormalSpecies.h"
#include "Interactions/NormalForceInteractions/SinterLinInteraction.h"
#include "Species/ParticleSpecies.h"
#include "Logger.h"

class BaseParticle;

class BaseInteractable;

SinterLinNormalSpecies::SinterLinNormalSpecies()
        : BaseNormalForce()
{
    loadingStiffness_ = 0.0;
    unloadingStiffnessMax_ = 0.0;
    cohesionStiffness_ = 0.0;
    penetrationDepthMax_ = 0.0;
    dissipation_ = 0.0;
    sinterAdhesion_ = 0.0;
    sinterRate_ = 0.0;
    complianceZero_ = 0.0;
    surfTension_ = 0.0;
    fluidity = 0.0;
    separationDis_ = 0.0;
    sinterType_ = SINTER_APPROACH::FRENKEL;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SinterLinNormalSpecies::SinterLinNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
SinterLinNormalSpecies::SinterLinNormalSpecies(const SinterLinNormalSpecies& p)
        : BaseNormalForce(p)
{
    loadingStiffness_ = p.loadingStiffness_;
    unloadingStiffnessMax_ = p.unloadingStiffnessMax_;
    cohesionStiffness_ = p.cohesionStiffness_;
    penetrationDepthMax_ = p.penetrationDepthMax_;
    dissipation_ = p.dissipation_;
    sinterAdhesion_ = p.sinterAdhesion_;
    sinterRate_ = p.sinterRate_;
    complianceZero_ = p.complianceZero_;
    surfTension_ = p.surfTension_;
    fluidity = p.fluidity;
    separationDis_ = p.separationDis_;
    sinterType_ = p.sinterType_;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SinterLinNormalSpecies::SinterLinNormalSpecies(const SinterLinNormalSpecies &p) finished"<<std::endl;
#endif
}

SinterLinNormalSpecies::~SinterLinNormalSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SinterNormalSpecies::~SinterNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void SinterLinNormalSpecies::write(std::ostream& os) const
{
    os << " loadingStiffness " << loadingStiffness_;
    os << " maxUnloadingStiffness " << unloadingStiffnessMax_;
    os << " cohesionStiffness " << cohesionStiffness_;
    os << " maxPenetration " << penetrationDepthMax_;
    os << " dissipation " << dissipation_;
    os << " sinterAdhesion " << sinterAdhesion_;
    os << " sinterRate " << sinterRate_;
    os << " complianceZero " << complianceZero_;
    os << " surfTension  " << surfTension_;
    os << " complianceOne " << fluidity;
    os << " separationDis " << separationDis_;
    os << " sinterType " << (unsigned) sinterType_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void SinterLinNormalSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> loadingStiffness_;
    is >> dummy >> unloadingStiffnessMax_;
    is >> dummy >> cohesionStiffness_;
    is >> dummy >> penetrationDepthMax_;
    is >> dummy >> dissipation_;
    is >> dummy >> sinterAdhesion_;
    is >> dummy >> complianceZero_;
    is >> dummy >> surfTension_;
    is >> dummy >> sinterRate_;
    is >> dummy >> fluidity;
    is >> dummy >> separationDis_;
    unsigned type;
    is >> dummy >> type;

    sinterType_ = (SINTER_APPROACH) type;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string SinterLinNormalSpecies::getBaseName() const
{
    return "SinterLin";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void SinterLinNormalSpecies::mix(SinterLinNormalSpecies* const S,
                                      SinterLinNormalSpecies* const T)
{
    loadingStiffness_ = BaseSpecies::average(S->getLoadingStiffness(), T->getLoadingStiffness());
    unloadingStiffnessMax_ = BaseSpecies::average(S->getUnloadingStiffnessMax(), T->getUnloadingStiffnessMax());
    cohesionStiffness_ = BaseSpecies::average(S->getCohesionStiffness(), T->getCohesionStiffness());
    penetrationDepthMax_ = BaseSpecies::average(S->getPenetrationDepthMax(), T->getPenetrationDepthMax());
    dissipation_ = BaseSpecies::average(S->getDissipation(), T->getDissipation());
    sinterAdhesion_ = BaseSpecies::average(S->getSinterAdhesion(), T->getSinterAdhesion());
    complianceZero_ = BaseSpecies::average(S->getComplianceZero(), T->getComplianceZero());
    surfTension_ = BaseSpecies::average(S->getSurfTension(), T->getSurfTension());
    fluidity= BaseSpecies::average(S->getFluidity(), T->getFluidity());
    separationDis_ = BaseSpecies::average(S->getSeparationDis(), T->getSeparationDis());

}

/*!
 * \param[in] loadingStiffness      the loading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] unloadingStiffnessMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] cohesionStiffness     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void SinterLinNormalSpecies::setPlasticParameters(Mdouble loadingStiffness, Mdouble unloadingStiffnessMax,
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
Mdouble SinterLinNormalSpecies::getLoadingStiffness() const
{
    return loadingStiffness_;
}

/*!
 * \return the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble SinterLinNormalSpecies::getUnloadingStiffnessMax() const
{
    return unloadingStiffnessMax_;
}

/*!
 * \return the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
Mdouble SinterLinNormalSpecies::getCohesionStiffness() const
{
    return cohesionStiffness_;
}

/*!
 * \return the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
Mdouble SinterLinNormalSpecies::getPenetrationDepthMax() const
{
    return penetrationDepthMax_;
}

/*!
 * \param[in] loadingStiffness      the loading stiffness of the linear plastic-viscoelastic normal force.
 */
void SinterLinNormalSpecies::setLoadingStiffness(Mdouble loadingStiffness)
{
    loadingStiffness_ = loadingStiffness;
}

/*!
 * \param[in] unloadingStiffnessMax the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
 */
void SinterLinNormalSpecies::setUnloadingStiffnessMax(Mdouble unloadingStiffnessMax)
{
    unloadingStiffnessMax_ = unloadingStiffnessMax;
}

/*!
 * \param[in] cohesionStiffness     the cohesive stiffness of the linear plastic-viscoelastic normal force.
 */
void SinterLinNormalSpecies::setCohesionStiffness(Mdouble cohesionStiffness)
{
    cohesionStiffness_ = cohesionStiffness;
}

/*!
 * \param[in] penetrationDepthMax   the maximum penetration depth of the linear plastic-viscoelastic normal force.
 */
void SinterLinNormalSpecies::setPenetrationDepthMax(Mdouble penetrationDepthMax)
{
    penetrationDepthMax_ = penetrationDepthMax;
}

/*!
 * \details Calculates collision time for stiffest spring constant, divides by 50
 * \param[in] mass the optimal time step is computed to resolve a collision of two particles of this mass.
 */

Mdouble SinterLinNormalSpecies::computeTimeStep(Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    return 0.02 * constants::pi /
           std::sqrt(unloadingStiffnessMax_ / (.5 * mass) - mathsFunc::square(dissipation_ / mass));
}

/*!
 * \details should be non-negative
 * \param[in] dissipation the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
 */
void SinterLinNormalSpecies::setDissipation(Mdouble dissipation)
{
    if (dissipation < 0)
    {
        logger(ERROR, "setDissipation(%)", dissipation);
        exit(-1);
    }
    dissipation_ = dissipation;
}


void SinterLinNormalSpecies::setLoadingStiffnessAndDissipation(helpers::KAndDisp new_)
{
    setLoadingStiffness(new_.k);
    setDissipation(new_.disp);
}

void SinterLinNormalSpecies::setSinterAdhesion(Mdouble sinterAdhesion)
{
    if (sinterAdhesion < 0)
    {
        logger(ERROR, "setSinterAdhesion(%)", sinterAdhesion);
        exit(-1);
    }
    sinterAdhesion_ = sinterAdhesion;
}

void SinterLinNormalSpecies::setSinterRate(Mdouble sinterRate)
{
    if (sinterRate < 0)
    {
        logger(ERROR, "setSinterRate(%)", sinterRate);
    }
    sinterRate_ = sinterRate;
}

void SinterLinNormalSpecies::setComplianceZero(Mdouble complianceZero)
{
    if (complianceZero < 0)
    {
        logger(ERROR, "complianceZero(%)", complianceZero);
    }
    complianceZero_ = complianceZero;
}

void SinterLinNormalSpecies::setSurfTension(Mdouble surfTension)
{
    if (surfTension < 0)
    {
        logger(ERROR, "SurfTension(%)", surfTension);
    }
    surfTension_ = surfTension;
}

void SinterLinNormalSpecies::setFluidity(Mdouble complianceOne)
{
    if (complianceOne < 0)
    {
        logger(ERROR, "ComplianceOne(%)", complianceOne);
    }
    fluidity = complianceOne;
}

void SinterLinNormalSpecies::setSeparationDis(Mdouble separationDis)
{
    if (separationDis < 0)
    {
        logger(ERROR, "SeparationDis(%)", separationDis);
    }
    separationDis_ = separationDis;
}

void SinterLinNormalSpecies::setSinterType(SINTER_APPROACH sinterType)
{
    sinterType_ = sinterType;
    if (sinterType_ == SINTER_APPROACH::FRENKEL)
        logger(INFO, "Sintertype set to Frenkel");
    else if (sinterType_ == SINTER_APPROACH::VISCOELASTIC_CONTACT)
        logger(INFO, "Sintertype set to Three different mechanisms");
    else
        logger(ERROR, "Sintertype not understood");
}





/*!
 * \return the linear dissipation coefficient
 */
Mdouble SinterLinNormalSpecies::getDissipation() const
{
    return dissipation_;
}

/*!
 * \return value of sinterAdhesion_
 */
Mdouble SinterLinNormalSpecies::getSinterAdhesion() const
{
    return sinterAdhesion_;
}

/*!
 * \return value of sinterAdhesion_
 */
Mdouble SinterLinNormalSpecies::getSinterRate() const
{
    return sinterRate_;
}

Mdouble SinterLinNormalSpecies::getComplianceZero() const
{
    return complianceZero_;
}

Mdouble SinterLinNormalSpecies::getSurfTension() const
{
    return surfTension_;
}

Mdouble SinterLinNormalSpecies::getFluidity() const
{
    return fluidity;
}

Mdouble SinterLinNormalSpecies::getSeparationDis() const
{
    return separationDis_;
}

// Sinter functions:
SINTER_APPROACH SinterLinNormalSpecies::getSinterType() const
{
    return sinterType_;
}


/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of equal mass m
 * \param[in] tc collision time
 * \param[in] eps restitution coefficient
 * \param[in] mass effective particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 */
///\todo TW: check that the masses are described correctly here (m_eff or m_p?))

void SinterLinNormalSpecies::setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
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
void SinterLinNormalSpecies::setStiffnessAndRestitutionCoefficient(Mdouble stiffness, Mdouble eps, Mdouble mass)
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

void SinterLinNormalSpecies::setRestitutionCoefficient(Mdouble eps, Mdouble mass)
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

///Calculates collision time for two copies of a particle of given disp, k, mass
Mdouble SinterLinNormalSpecies::getCollisionTime(Mdouble mass) const
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

Mdouble SinterLinNormalSpecies::getRestitutionCoefficient(Mdouble mass) const
{
    if (getConstantRestitution()) mass = 1;

    return std::exp(-dissipation_ / mass * getCollisionTime(mass));
}