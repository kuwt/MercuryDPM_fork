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

#include "StandardFields.h"
#include <Particles/BaseParticle.h>

namespace CGFields
{

StandardFields::StandardFields()
{
    setZero();
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "StandardFields::StandardFields() finished" << std::endl;
#endif
}

/*!
 * \param[in,out] os the ostream into which the data is written.
 * \param[in] countVariables The number of variables in the field (including time), e.g. 1 for O, 4 for XYZ
 */
void StandardFields::writeNames(std::ostream& os, const unsigned countVariables)
{
    os << countVariables + 1 << ":VolumeFraction "; //volume
    os << countVariables + 2 << ":Density "; //mass
    os << countVariables + 3 << '-' << countVariables + 5 << ":Momentum ";
    os << countVariables + 6 << '-' << countVariables + 11 << ":MomentumFlux ";
    os << countVariables + 12 << '-' << countVariables + 20 << ":ContactStress ";
    os << countVariables + 21 << '-' << countVariables + 23 << ":InteractionForce ";
    os << countVariables + 24 << '-' << countVariables + 29 << ":ParticleSize ";
}

/*!
 * \param[in,out] os the ostream into which the data is written.
 */
void StandardFields::write(std::ostream& os) const
{
    os << volumeFraction_;
    os << ' ' << density_;
    os << ' ' << momentum_;
    os << ' ' << momentumFlux_;
    os << ' ' << contactStress_;
    os << ' ' << interactionForceDensity_;
    for (const auto p : particleSizeDensity_) os << ' ' << p;
}

/*!
 * \param[in,out] os the ostream into which the data is written.
 */
void StandardFields::output(std::ostream& os) const
{
    os << "VolumeFraction " << volumeFraction_;
    os << " Density " << density_;
    os << " Momentum " << momentum_;
    os << " MomentumFlux " << momentumFlux_;
    os << " ContactStress " << contactStress_;
    os << " InteractionForce " << interactionForceDensity_;
    os << " ParticleSize";
    for (const auto p : particleSizeDensity_) os << ' ' << p;
}

void StandardFields::setZero()
{
    volumeFraction_ = 0.0;
    density_ = 0.0;
    momentum_.setZero();
    momentumFlux_.setZero();
    contactStress_.setZero();
    interactionForceDensity_.setZero();
    for (auto& p : particleSizeDensity_) p = 0;
}

/*!
 * \return a CGField containing the square of the values in the current object
 */
StandardFields StandardFields::getSquared() const
{
    StandardFields P;
    P.volumeFraction_ = mathsFunc::square(volumeFraction_);
    P.density_ = mathsFunc::square(density_);
    P.momentum_ = Vec3D::square(momentum_);
    P.momentumFlux_ = MatrixSymmetric3D::square(momentumFlux_);
    P.contactStress_ = Matrix3D::square(contactStress_);
    P.interactionForceDensity_ = Vec3D::square(interactionForceDensity_);
    P.particleSizeDensity_ = particleSizeDensity_;
    for (auto& p : P.particleSizeDensity_) p *= p;
    return P;
}

/*!
 * \param[in] P the CGField that has to be copied
 * \return the CGField into which the values are copied
 */
StandardFields& StandardFields::operator=(const StandardFields& P)
= default;

/*!
* \param[in] P the CGField that has to be added
* \return the CGField to which the values are added
*/
StandardFields& StandardFields::operator+=(const StandardFields& P)
{
    volumeFraction_ += P.volumeFraction_;
    density_ += P.density_;
    momentum_ += P.momentum_;
    momentumFlux_ += P.momentumFlux_;
    contactStress_ += P.contactStress_;
    interactionForceDensity_ += P.interactionForceDensity_;
    for (size_t i = 0; i < particleSizeDensity_.size(); ++i) particleSizeDensity_[i] += P.particleSizeDensity_[i];
    return *this;
}

/*!
 * \param[in] P the CGField that has to be subtracted
 * \return the CGField from which the values are subtracted
 */
StandardFields& StandardFields::operator-=(const StandardFields& P)
{
    volumeFraction_ -= P.volumeFraction_;
    density_ -= P.density_;
    momentum_ -= P.momentum_;
    momentumFlux_ -= P.momentumFlux_;
    contactStress_ -= P.contactStress_;
    interactionForceDensity_ -= P.interactionForceDensity_;
    for (size_t i = 0; i < particleSizeDensity_.size(); ++i) particleSizeDensity_[i] -= P.particleSizeDensity_[i];
    return *this;
}

/*!
 * \param[in] a the scalar that we multiply with
 * \return the CGField  to which the multiplied values are written
 */
StandardFields StandardFields::operator*(const Mdouble a) const
{
    StandardFields p;
    p.volumeFraction_ = volumeFraction_ * a;
    p.density_ = density_ * a;
    p.momentum_ = momentum_ * a;
    p.momentumFlux_ = momentumFlux_ * a;
    p.contactStress_ = contactStress_ * a;
    p.interactionForceDensity_ = interactionForceDensity_ * a;
    for (size_t i = 0; i < particleSizeDensity_.size(); ++i)
        p.particleSizeDensity_[i] = particleSizeDensity_[i] * a;
    return p;
}

/*!
 * \param[in] a the scalar that we divide by
 * \return the CGField to which the divided values are written
 */
StandardFields& StandardFields::operator/=(const Mdouble a)
{
    volumeFraction_ /= a;
    density_ /= a;
    momentum_ /= a;
    momentumFlux_ /= a;
    contactStress_ /= a;
    interactionForceDensity_ /= a;
    for (auto& p : particleSizeDensity_) p /= a;
    return *this;
}

/*!
 * \param[in] phi the value of the cg function at the current CGPoint
 * \param[in] currentInteraction the fields which are produced due to the particle at its centre
 */
void StandardFields::addParticleStatistics(Mdouble phi, const StandardFields& currentInteraction)
{
    volumeFraction_ += currentInteraction.getVolumeFraction() * phi;
    density_ += currentInteraction.getDensity() * phi;
    momentum_ += currentInteraction.getMomentum() * phi;
    momentumFlux_ += currentInteraction.getMomentumFlux() * phi;
    for (size_t i = 0; i < particleSizeDensity_.size(); ++i)
        particleSizeDensity_[i] += currentInteraction.getParticleSizeDensity(i) * phi;
}

void StandardFields::addParticleDifferentialStatistics(Vec3D& dphi, const StandardFields& currentInteraction)
{
}

/*!
 * \param[in] psi the value of the line integral from C to P at the current CGPoint
 * \param[in] c the contact which is used in the line integral
 */
void StandardFields::addInteractionStatistics(Mdouble psi, const StandardFields& currentInteraction)
{
    contactStress_ += currentInteraction.getContactStress() * psi;
}

/*!
 * \param[in] phi the value of the cg function for the contact point of c and 
 * the current CGPoint
 * \param[in] c the interaction which is used in the cg function
 */
void StandardFields::addContactPointStatistics(Mdouble phi, const StandardFields& currentInteraction)
{
    interactionForceDensity_ += currentInteraction.getInteractionForceDensity() * phi;
}

/*!
 * \details If the functions returns false, addInteractionStatistics and 
 * addContactPointStatistics are never called.
 * \return True if the class contains fields that are defined as a sum over all 
 * Interactions, else false.
 */
bool StandardFields::doInteractionStatistics()
{
    return true;
}

void StandardFields::setFields(const BaseInteraction& c, IntegralType type)
{
    //P is always real, but I not
    if (type == IntegralType::I_TO_P)
    {
        contactStress_ = Matrix3D::dyadic(c.getForce(), c.getIP());
    }
    else if (type == IntegralType::I_TO_CONTACT)
    {
        contactStress_ = Matrix3D::dyadic(c.getForce(), c.getIC());
    }
    else
    {
        contactStress_ = Matrix3D::dyadic(c.getForce(), c.getCP());
    }
    interactionForceDensity_ = c.getForce();
}

void StandardFields::setFields(const BaseParticle& p)
{
    volumeFraction_ = p.getVolume();
    density_ = p.getMass();
    momentum_ = p.getVelocity() * p.getMass();
    momentumFlux_ = MatrixSymmetric3D::selfDyadic(p.getVelocity()) * p.getMass();
    Mdouble particleSize = 1.0;
    for (auto& ps : particleSizeDensity_)
    {
        ps = particleSize;
        particleSize *= p.getRadius();
    }
}

void StandardFields::setCylindricalFields(const BaseInteraction& c, IntegralType type)
{
    setFields(c, type);
    contactStress_ = contactStress_.getCylindricalTensorField(c.getContactPoint());
    interactionForceDensity_ = interactionForceDensity_.getCylindricalTensorField(c.getContactPoint());
}

void StandardFields::setCylindricalFields(const BaseParticle& p)
{
    setFields(p);
    momentum_ = momentum_.getCylindricalTensorField(p.getPosition());
    momentumFlux_ = momentumFlux_.getCylindricalTensorField(p.getPosition());
}

void StandardFields::outputStandardisedParticleSizeMomenta(std::ostream& os) const
{
    //https://en.wikipedia.org/wiki/Moment_%28mathematics%29#Higher_moments
    auto momenta = getStandardisedParticleSizeMomenta();
    os << "Particle number: " << momenta[0] << '\n';
    os << "Mean radius: " << momenta[1] << '\n';
    os << "Standard deviation: " << std::sqrt(momenta[2]) << '\n';
    os << "Skewness: " << momenta[3] << '\n'; //left- or right-sided
    os << "Kurtosis: " << momenta[4] << '\n'; //heavy or light tailed (3 for normal distribution)
    os << "5-th moment: " << momenta[5] << '\n'; //hyper skewness (0 for normal dist)
}

std::array<Mdouble, 6> StandardFields::getParticleSizeMomenta() const
{
    //https://en.wikipedia.org/wiki/Moment_(mathematics)
    auto momenta = particleSizeDensity_;
    for (size_t i = 1; i < particleSizeDensity_.size(); ++i)
    {
        momenta[i] /= particleSizeDensity_[0];
    }
    return momenta;
}

std::array<Mdouble, 6> StandardFields::getCentralParticleSizeMomenta() const
{
    //https://en.wikipedia.org/wiki/Central_moment
    auto momenta = getParticleSizeMomenta();
    Mdouble mean = momenta[1];
    momenta[5] += -5 * mean * momenta[4] + 10 * mean * mean * momenta[3]
                  - 10 * mean * mean * mean * momenta[2] + 4 * mean * mean * mean * mean * mean;
    momenta[4] += -4 * mean * momenta[3] + 6 * mean * mean * momenta[2] - 3 * mean * mean * mean * mean;
    momenta[3] += -3 * mean * momenta[2] + 2 * mean * mean * mean;
    momenta[2] += -mean * mean;
    return momenta;
}

std::array<Mdouble, 6> StandardFields::getStandardisedParticleSizeMomenta() const
{
    auto momenta = getCentralParticleSizeMomenta();
    Mdouble std = std::sqrt(momenta[2]);
    momenta[3] /= std * std * std;
    momenta[4] /= std * std * std * std;
    momenta[5] /= std * std * std * std * std;
    return momenta;
}
    
}
