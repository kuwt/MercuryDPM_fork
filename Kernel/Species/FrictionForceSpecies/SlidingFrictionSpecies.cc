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

#include <cmath>
#include "SlidingFrictionSpecies.h"
#include "Species/NormalForceSpecies/LinearViscoelasticNormalSpecies.h"
#include "Species/BaseSpecies.h"
#include "Species/NormalForceSpecies/LinearPlasticViscoelasticNormalSpecies.h"

class BaseParticle;

class BaseInteractable;

SlidingFrictionSpecies::SlidingFrictionSpecies()
{
    slidingStiffness_ = 0;
    slidingDissipation_ = 0;
    slidingFrictionCoefficient_ = 0;
    slidingFrictionCoefficientStatic_ = 0;
    isSuperquadricSpecies_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SlidingFrictionSpecies::SlidingFrictionSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] s the species that is copied
 */
SlidingFrictionSpecies::SlidingFrictionSpecies(const SlidingFrictionSpecies& s)
{
    slidingStiffness_ = s.slidingStiffness_;
    slidingDissipation_ = s.slidingDissipation_;
    slidingFrictionCoefficient_ = s.slidingFrictionCoefficient_;
    slidingFrictionCoefficientStatic_ = s.slidingFrictionCoefficientStatic_;
    isSuperquadricSpecies_ = s.isSuperquadricSpecies_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SlidingFrictionSpecies::SlidingFrictionSpecies(const SlidingFrictionSpecies &p) finished"<<std::endl;
#endif
}

SlidingFrictionSpecies::~SlidingFrictionSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SlidingFrictionSpecies::~SlidingFrictionSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void SlidingFrictionSpecies::write(std::ostream& os) const
{
    //BaseSpecies::write(os);
    os << " slidingStiffness " << slidingStiffness_;
    os << " slidingDissipation " << slidingDissipation_;
    os << " slidingFrictionCoefficient " << slidingFrictionCoefficient_;
    os << " slidingFrictionCoefficientStatic " << slidingFrictionCoefficientStatic_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void SlidingFrictionSpecies::read(std::istream& is)
{
    //BaseSpecies::read(is);
    std::string dummy;
    is >> dummy >> slidingStiffness_;
    is >> dummy >> slidingDissipation_;
    is >> dummy >> slidingFrictionCoefficient_;
    is >> dummy >> slidingFrictionCoefficientStatic_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string SlidingFrictionSpecies::getBaseName() const
{
    return "SlidingFriction";
}

///Allows the spring constant to be changed
void SlidingFrictionSpecies::setSlidingStiffness(Mdouble new_kt)
{
    if (new_kt >= 0)
    {
        slidingStiffness_ = new_kt;
    }
    else
    {
        logger(ERROR, "Error in setSlidingStiffness");
    }
}

///Allows the spring constant to be accessed
Mdouble SlidingFrictionSpecies::getSlidingStiffness() const
{
    return slidingStiffness_;
}

///Allows the tangential viscosity to be changed
void SlidingFrictionSpecies::setSlidingDissipation(Mdouble new_dispt)
{
    if (new_dispt >= 0)
        slidingDissipation_ = new_dispt;
    else
    {
        logger(ERROR, "Error in setSlidingDissipation");
    }
}

///Allows the tangential viscosity to be accessed
Mdouble SlidingFrictionSpecies::getSlidingDissipation() const
{
    return slidingDissipation_;
}

///Allows the (dynamic) Coulomb friction coefficient to be changed; also sets mu_s by default
//mu has to be set to allow tangential forces (sets dispt=disp as default)
void SlidingFrictionSpecies::setSlidingFrictionCoefficient(Mdouble new_mu)
{
    if (new_mu >= 0)
    {
        slidingFrictionCoefficient_ = new_mu;
        slidingFrictionCoefficientStatic_ = slidingFrictionCoefficient_;
    }
    else
    {
        logger(ERROR, "Error in setSlidingFrictionCoefficient");
    }
}

///Allows the (dynamic) Coulomb friction coefficient to be accessed
Mdouble SlidingFrictionSpecies::getSlidingFrictionCoefficient() const
{
    return slidingFrictionCoefficient_;
}

///Allows the static Coulomb friction coefficient to be changed; also sets mu_s by default
void SlidingFrictionSpecies::setSlidingFrictionCoefficientStatic(Mdouble new_mu)
{
    if (new_mu >= 0)
    {
        slidingFrictionCoefficientStatic_ = new_mu;
    }
    else
    {
        logger(ERROR, "Error in setSlidingFrictionCoefficientStatic");
    }
}

///Allows the static Coulomb friction coefficient to be accessed
Mdouble SlidingFrictionSpecies::getSlidingFrictionCoefficientStatic() const
{
    return slidingFrictionCoefficientStatic_;
}

/*!
 * \details Returns true for any FrictionForceSpecies except EmptyFrictionSpecies, 
 * because for spherical particles, torques are only caused by tangential forces. 
 * See SpeciesHandler::useAngularDOFs for more details
 * \return true 
 */
bool SlidingFrictionSpecies::getUseAngularDOFs() const
{
    return true;
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] SFrictional the first species whose properties are mixed to create the new species
 * \param[in] TFrictional the second species whose properties are mixed to create the new species
 */
void SlidingFrictionSpecies::mix(SlidingFrictionSpecies* const S, SlidingFrictionSpecies* const T)
{
    slidingStiffness_ = BaseSpecies::average(S->getSlidingStiffness(), T->getSlidingStiffness());
    slidingDissipation_ = BaseSpecies::average(S->getSlidingDissipation(), T->getSlidingDissipation());
    slidingFrictionCoefficient_ = BaseSpecies::average(S->getSlidingFrictionCoefficient(),
                                          T->getSlidingFrictionCoefficient());
    slidingFrictionCoefficientStatic_ = BaseSpecies::average(S->getSlidingFrictionCoefficientStatic(),
                                                T->getSlidingFrictionCoefficientStatic());
}

///Sets k, disp, kt, dispt such that it matches a given tc and eps for a collision of two particles of masses m0,m1
void SlidingFrictionSpecies::setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(Mdouble tc, Mdouble eps,
                                                                                          Mdouble beta, Mdouble mass)
{
    Mdouble stiffness;
    //the dynamic cast is needed to check if the normal force species is LinearViscoelasticSpecies; otherwise this function cannot be used
    LinearViscoelasticNormalSpecies* species = dynamic_cast<LinearViscoelasticNormalSpecies*>(this);
    if (species != nullptr)
    {
        species->setCollisionTimeAndRestitutionCoefficient(tc, eps, mass);
        stiffness = species->getStiffness();
    }
    else
    {
        LinearPlasticViscoelasticNormalSpecies* species2 = dynamic_cast<LinearPlasticViscoelasticNormalSpecies*>(this);
        if (species2 != nullptr)
        {
            species2->setDissipation(-mass / tc * std::log(eps));
            species2->setLoadingStiffness(.5 * mass * (mathsFunc::square(constants::pi / tc) +
                                                       mathsFunc::square(species2->getDissipation()) / mass));
            stiffness = species2->getLoadingStiffness();
        }
        else
        {
            logger(ERROR,
                   "SlidingFrictionSpecies::setCollisionTimeAndNormalAndTangentialRestitutionCoefficient only works for "
                   "LinearViscoelasticSlidingFrictionSpecies or LinearPlasticViscoelasticSlidingFrictionSpecies");
        }
    }
    
    // from: N. G. Deen et. al. https://doi.org/10.1016/j.ces.2006.08.014
    // eq. 43 and 30
    setSlidingStiffness(2.0 / 7.0 * stiffness * (mathsFunc::square(constants::pi) + mathsFunc::square(log(beta))) /
                        (mathsFunc::square(constants::pi) + mathsFunc::square(log(eps))));
    if (beta != 0.0)
        setSlidingDissipation(-2 * log(beta) * sqrt(1.0 / 7.0 * mass * getSlidingStiffness() /
                                                    (mathsFunc::square(constants::pi) + mathsFunc::square(log(beta)))));
    else
        setSlidingDissipation(2. * sqrt(1.0 / 7.0 * mass * getSlidingStiffness()));
}

///Sets k, disp, kt, dispt such that it matches a given tc and eps for a collision of two particles of masses m0,m1
void
SlidingFrictionSpecies::setCollisionTimeAndNormalAndTangentialRestitutionCoefficientNoDispt(Mdouble tc, Mdouble eps,
                                                                                            Mdouble beta, Mdouble mass)
{
    //the dynamic cast is needed to check if the normal force species is LinearViscoelasticSpecies; otherwise this function cannot be used
    LinearViscoelasticNormalSpecies* species = dynamic_cast<LinearViscoelasticNormalSpecies*>(this);
    if (species == nullptr)
    {
        logger(ERROR, "SlidingFrictionSpecies::setCollisionTimeAndNormalAndTangentialRestitutionCoefficient only "
                      "works for LinearViscoelasticSlidingFrictionSpecies");
    }
    species->setCollisionTimeAndRestitutionCoefficient(tc, eps, mass);
    // from: V. Becker et. al. https://doi.org/10.1103/PhysRevE.77.011304
    // eq. 56
    setSlidingStiffness(2.0 / 7.0 * species->getStiffness() * mathsFunc::square(acos(-beta) / constants::pi));
    setSlidingDissipation(0);
}

void SlidingFrictionSpecies::setIsSuperquadricSpecies(bool isSuperquadricSpecies)
{
    SlidingFrictionSpecies::isSuperquadricSpecies_ = isSuperquadricSpecies;
}

bool SlidingFrictionSpecies::getIsSuperquadricSpecies() const
{
    return isSuperquadricSpecies_;
}
