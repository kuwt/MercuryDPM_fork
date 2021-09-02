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

#include <cmath>
#include <Species/NormalForceSpecies/HertzianViscoelasticNormalSpecies.h>
#include "MindlinSpecies.h"
#include "Species/NormalForceSpecies/LinearViscoelasticNormalSpecies.h"
#include "Species/BaseSpecies.h"
#include "Species/NormalForceSpecies/LinearPlasticViscoelasticNormalSpecies.h"

class BaseParticle;

class BaseInteractable;

MindlinSpecies::MindlinSpecies()
{
    slidingDissipation_ = 0;
    slidingFrictionCoefficient_ = 0;
    slidingFrictionCoefficientStatic_ = 0;
    shearModulus_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinSpecies::MindlinSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] s the species that is copied
 */
MindlinSpecies::MindlinSpecies(const MindlinSpecies& s)
{
    slidingDissipation_ = s.slidingDissipation_;
    slidingFrictionCoefficient_ = s.slidingFrictionCoefficient_;
    slidingFrictionCoefficientStatic_ = s.slidingFrictionCoefficientStatic_;
    shearModulus_ = s.shearModulus_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinSpecies::MindlinSpecies(const MindlinSpecies &p) finished"<<std::endl;
#endif
}

MindlinSpecies::~MindlinSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"MindlinSpecies::~MindlinSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] os output stream (typically the restart file)
 *
 */
void MindlinSpecies::write(std::ostream& os) const
{
    //BaseSpecies::write(os);
    os << " slidingDissipation " << slidingDissipation_;
    os << " slidingFrictionCoefficient " << slidingFrictionCoefficient_;
    os << " slidingFrictionCoefficientStatic " << slidingFrictionCoefficientStatic_;
    os << " shearModulus " << shearModulus_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void MindlinSpecies::read(std::istream& is)
{
    //BaseSpecies::read(is);
    std::string dummy;
    is >> dummy >> slidingDissipation_;
    is >> dummy >> slidingFrictionCoefficient_;
    is >> dummy >> slidingFrictionCoefficientStatic_;
    is >> dummy >> shearModulus_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string MindlinSpecies::getBaseName() const
{
    return "Mindlin";
}


///Allows the tangential viscosity to be changed
void MindlinSpecies::setSlidingDissipation(Mdouble new_dispt)
{
    if (new_dispt >= 0)
        slidingDissipation_ = new_dispt;
    else
    {
        std::cerr << "Error in setSlidingDissipation" << std::endl;
        exit(-1);
    }
}

///Allows the tangential viscosity to be accessed
Mdouble MindlinSpecies::getSlidingDissipation() const
{
    return slidingDissipation_;
}

///Allows the (dynamic) Coulomb friction coefficient to be changed; also sets mu_s by default
//mu has to be set to allow tangential forces (sets dispt=disp as default)
void MindlinSpecies::setSlidingFrictionCoefficient(Mdouble new_mu)
{
    if (new_mu >= 0)
    {
        slidingFrictionCoefficient_ = new_mu;
        slidingFrictionCoefficientStatic_ = slidingFrictionCoefficient_;
    }
    else
    {
        std::cerr << "Error in setSlidingFrictionCoefficient" << std::endl;
        exit(-1);
    }
}

///Allows the (dynamic) Coulomb friction coefficient to be accessed
Mdouble MindlinSpecies::getSlidingFrictionCoefficient() const
{
    return slidingFrictionCoefficient_;
}

///Allows the static Coulomb friction coefficient to be changed
void MindlinSpecies::setSlidingFrictionCoefficientStatic(Mdouble new_mu)
{
    if (new_mu >= 0)
    {
        slidingFrictionCoefficientStatic_ = new_mu;
    }
    else
    {
        std::cerr << "Error in setSlidingFrictionCoefficientStatic" << std::endl;
        exit(-1);
    }
}

///Allows the poisson ratio to be changed
void MindlinSpecies::setPoissonRatio(Mdouble poissonRatio)
{
    auto hertz = dynamic_cast<HertzianViscoelasticNormalSpecies*>(getBaseSpecies());
    logger.assert(hertz, "setPoissonRatio only works with HertzianViscoelastic*Species.");
    shearModulus_ = hertz->getEffectiveElasticModulus() / 2 * (1 + poissonRatio);
    logger(INFO, "Effective shear modulus set to %, based on effective elastic modulus % and Poisson's ratio %",
           getEffectiveShearModulus(),
           hertz->getEffectiveElasticModulus(), poissonRatio);
}

/////Allows the poisson ratio to be accessed
//Mdouble MindlinSpecies::getPoissonRatio() const
//{
//    auto hertz = dynamic_cast<HertzianViscoelasticNormalSpecies*>(getBaseSpecies());
//    logger.assert(hertz, "getPoissonRatio only works with HertzianViscoelastic*Species.");
//    return 2. * shearModulus_ / hertz->getEffectiveElasticModulus() - 1;
//}

///Allows the shear modulus to be changed
void MindlinSpecies::setEffectiveShearModulus(Mdouble shearModulus)
{
    logger.assert(shearModulus > 0,
                  "setEffectiveShearModulus(%): value needs to be positive",shearModulus);
    shearModulus_ = shearModulus;
}

///Allows the shear modulus to be accessed
Mdouble MindlinSpecies::getEffectiveShearModulus() const
{
    return shearModulus_;
}

///Allows the static Coulomb friction coefficient to be accessed
Mdouble MindlinSpecies::getSlidingFrictionCoefficientStatic() const
{
    return slidingFrictionCoefficientStatic_;
}

/*!
 * \details Returns true for any FrictionForceSpecies except EmptyFrictionSpecies, 
 * because for spherical particles, torques are only caused by tangential forces. 
 * See SpeciesHandler::useAngularDOFs for more details
 * \return true 
 */
bool MindlinSpecies::getUseAngularDOFs() const
{
    return true;
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] SFrictional the first species whose properties are mixed to create the new species
 * \param[in] TFrictional the second species whose properties are mixed to create the new species
 */
void MindlinSpecies::mix(MindlinSpecies* const SFrictional, MindlinSpecies* const TFrictional)
{
    slidingDissipation_ = BaseSpecies::average(SFrictional->getSlidingDissipation(),TFrictional->getSlidingDissipation());
    slidingFrictionCoefficient_ = BaseSpecies::average(SFrictional->getSlidingFrictionCoefficient(),TFrictional->getSlidingFrictionCoefficient());
    slidingFrictionCoefficientStatic_ = BaseSpecies::average(SFrictional->getSlidingFrictionCoefficientStatic(),TFrictional->getSlidingFrictionCoefficientStatic());
    shearModulus_ = BaseSpecies::average(SFrictional->getEffectiveShearModulus(), TFrictional->getEffectiveShearModulus());
}

