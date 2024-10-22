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

#include "LinearViscoelasticNormalSpecies.h"
#include "Interactions/BaseInteraction.h"
#include<cmath>
#include <Interactions/NormalForceInteractions/LinearViscoelasticInteraction.h>
#include "Logger.h"
#include "Particles/BaseParticle.h"

class BaseInteractable;

LinearViscoelasticNormalSpecies::LinearViscoelasticNormalSpecies()
        : BaseNormalForce()
{
    stiffness_ = 0.0;
    dissipation_ = 0.0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearViscoelasticNormalSpecies::LinearViscoelasticNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
LinearViscoelasticNormalSpecies::LinearViscoelasticNormalSpecies(const LinearViscoelasticNormalSpecies& p)
        : BaseNormalForce(p)
{
    stiffness_ = p.stiffness_;
    dissipation_ = p.dissipation_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearViscoelasticNormalSpecies::LinearViscoelasticNormalSpecies(const LinearViscoelasticNormalSpecies &p) finished"<<std::endl;
#endif
}

LinearViscoelasticNormalSpecies::~LinearViscoelasticNormalSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LinearViscoelasticNormalSpecies::~LinearViscoelasticNormalSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void LinearViscoelasticNormalSpecies::write(std::ostream& os) const
{
    os << " stiffness " << stiffness_
       << " dissipation " << dissipation_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void LinearViscoelasticNormalSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> stiffness_
       >> dummy >> dissipation_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string LinearViscoelasticNormalSpecies::getBaseName() const
{
    return "LinearViscoelastic";
}

///Allows the spring constant to be changed
void LinearViscoelasticNormalSpecies::setStiffness(Mdouble new_k)
{
    if (new_k >= 0)
        stiffness_ = new_k;
    else
    {
        logger(ERROR, "setStiffness(%) argument has to be non-negative!", new_k);
    }
}

///Allows the spring constant to be accessed
Mdouble LinearViscoelasticNormalSpecies::getStiffness() const
{
    return stiffness_;
}

///Allows the spring and dissipation constants to be changed simultaneously
void LinearViscoelasticNormalSpecies::setStiffnessAndDissipation(helpers::KAndDisp new_)
{
    setStiffness(new_.k);
    setDissipation(new_.disp);
}

///Allows the normal dissipation to be changed
void LinearViscoelasticNormalSpecies::setDissipation(Mdouble dissipation)
{
    if (dissipation >= 0)
    {
        dissipation_ = dissipation;
    }
    else
    {
        logger(ERROR, "Error in setDissipation(%)", dissipation);
    }
}

///Allows the normal dissipation to be accessed
Mdouble LinearViscoelasticNormalSpecies::getDissipation() const
{
    return dissipation_;
}

///Calculates collision time for two copies of a particle of given disp, k, mass
///\param[in] mass mass of a typical particle
Mdouble LinearViscoelasticNormalSpecies::getCollisionTime(Mdouble mass) const
{
    if (getConstantRestitution()) mass = 1;
    if (mass <= 0)
    {
        logger(ERROR, "Warning in getCollisionTime(%) mass is not set or has an unexpected value, "
                      "(getCollisionTime(%)", mass, mass);
    }
    if (stiffness_ <= 0)
    {
        logger(ERROR, "Error in getCollisionTime(%) stiffness=% is not set or has an unexpected value "
                      "(getCollisionTime(%), with stiffness=%)",
               mass, stiffness_, mass, stiffness_);
    }
    if (dissipation_ < 0)
    {
        logger(ERROR, "Error in getCollisionTime(%) dissipation=% is not set or has an unexpected value "
                      "(getCollisionTime(%), with dissipation=%)",
               mass, dissipation_, mass, dissipation_);
    }
    Mdouble tosqrt = stiffness_ / (.5 * mass) - mathsFunc::square(dissipation_ / mass);
    if (tosqrt <= 0)
    {
        logger(ERROR, "Warning in getCollisionTime(%) values for mass, stiffness and dissipation would lead to an "
                      "overdamped system, (getCollisionTime(%), with stiffness=% and dissipation=%",
               mass, mass, stiffness_, dissipation_);
    }
    return constants::pi / std::sqrt(tosqrt);
}

///Calculates restitution coefficient for two copies of given disp, k, mass
Mdouble LinearViscoelasticNormalSpecies::getRestitutionCoefficient(Mdouble mass) const
{
    if (getConstantRestitution()) mass = 1;
    return std::exp(-dissipation_ / mass * getCollisionTime(mass));
}

///Calculates the maximum velocity allowed for a collision of two copies of P (for higher velocities particles could pass through each other)
Mdouble LinearViscoelasticNormalSpecies::getMaximumVelocity(Mdouble radius, Mdouble mass) const
{
    return radius * std::sqrt(stiffness_ / (.5 * mass));
}

/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of P
 * \param[in] stiffness stiffness
 * \param[in] eps restitution coefficient
 * \param[in] harmonic mean of particle masses, \f$\frac{2}{1/m1+1/m2}\f$
 */
void LinearViscoelasticNormalSpecies::setStiffnessAndRestitutionCoefficient(Mdouble stiffness, Mdouble eps, Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    stiffness_ = stiffness;
    setRestitutionCoefficient(eps, mass);
}

/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of P
 * \param[in] stiffness stiffness
 * \param[in] eps restitution coefficient
 * \param[in] mass effective particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 */
void LinearViscoelasticNormalSpecies::setRestitutionCoefficient(Mdouble eps, Mdouble mass)
{
    if (getConstantRestitution()) mass = 1;
    if (eps == 0.0) {
        dissipation_ = std::sqrt(2.0 * mass * getStiffness());
    } else {
        const Mdouble logEps = log(eps);
        dissipation_ = -std::sqrt(2.0 * mass * getStiffness()
                                  / (constants::sqr_pi + mathsFunc::square(logEps))) * logEps;
    }
}

void
LinearViscoelasticNormalSpecies::setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, BaseParticle* p)
{
    Mdouble mass = p->getSpecies()->getDensity() * p->getVolume();
    setCollisionTimeAndRestitutionCoefficient(tc, eps, mass);
}


/*!
 * \details Sets k, disp such that it matches a given tc and eps for a collision of two copies of equal mass m
 * \param[in] tc collision time
 * \param[in] eps restitution coefficient
 * \param[in] mass harmonic average particle mass, \f$\frac{2}{1/m1+1/m2}\f$
 */
void LinearViscoelasticNormalSpecies::setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
{
    logger.assert_always(tc>0, "collision time has to be positive");
    logger.assert_always(eps>=0 && eps<=1, "restitution has to be in [0,1]");
    logger.assert_always(mass>0, "mass has to be positive");
    if (getConstantRestitution()) mass = 1;
    if (eps == 0.0)
    {
        stiffness_ = .5 * mass * mathsFunc::square(constants::pi / tc);
        dissipation_ = std::sqrt(2.0 * mass * stiffness_);
    }
    else
    {
        dissipation_ = -mass / tc * mathsFunc::log(eps);
        stiffness_ = .5 * mass * (mathsFunc::square(constants::pi / tc)
                                  + mathsFunc::square(dissipation_ / mass));
    }
    logger(INFO,
           "setCollisionTimeAndRestitutionCoefficient: set stiffness to % and dissipation to % to obtain a collision time of % and a restitution coefficient of % for two particles of mass %",
           getStiffness(), getDissipation(), tc, eps, mass);
}

/*!
 * \details Set k, disp such that is matches a given tc and eps for a collision of two different masses.
 * Recall the resitution constant is a function of k, disp and the mass of each particle in the collision
 * See also setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
 * \param[in] collision time
 */
void LinearViscoelasticNormalSpecies::setCollisionTimeAndRestitutionCoefficient(Mdouble collisionTime,
                                                                                Mdouble restitutionCoefficient,
                                                                                Mdouble mass1, Mdouble mass2)
{
    Mdouble reduced_mass = mass1 * mass2 / (mass1 + mass2);
    setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, 2.0 * reduced_mass);
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void
LinearViscoelasticNormalSpecies::mix(LinearViscoelasticNormalSpecies* const S, LinearViscoelasticNormalSpecies* const T)
{
    stiffness_ = BaseSpecies::average(S->getStiffness(), T->getStiffness());
    dissipation_ = BaseSpecies::average(S->getDissipation(), T->getDissipation());
}
