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

#include "GradVelocityField.h"
#include <Particles/BaseParticle.h>

namespace CGFields
{

GradVelocityField::GradVelocityField()
{
    setZero();
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "GradVelocityField::GradVelocityField() finished" << std::endl;
#endif
}

/*!
 * \param[out] os the ostream into which the data is written.
 */
void GradVelocityField::writeNames(std::ostream& os, const unsigned countVariables)
{
    os << countVariables + 1 << ":Density ";
    os << countVariables + 2 << "-" << countVariables + 4 << ":momentum ";
    os << countVariables + 5 << "-" << countVariables + 7 << ":ddensity ";
    os << countVariables + 8 << "-" << countVariables + 16 << ":dmomentum ";
}

/*!
 * \param[out] os the ostream into which the data is written.
 */
void GradVelocityField::write(std::ostream& os) const
{
    os << density_;
    os << " " << momentum_;
    os << " " << ddensity_;
    os << " " << dmomentum_;
}

/*!
 * \param[out] os the ostream into which the data is written.
 */
void GradVelocityField::output(std::ostream& os) const
{
    os << "density " << density_;
    os << " momentum " << momentum_;
    os << " ddensity " << ddensity_;
    os << " dmomentum " << dmomentum_;
}

void GradVelocityField::setZero()
{
    density_ = 0.0;
    momentum_.setZero();
    ddensity_.setZero();
    dmomentum_.setZero();
}

/*!
 * \return a CGField containing the square of the values in the current object
 */
GradVelocityField GradVelocityField::getSquared() const
{
    GradVelocityField P;
    P.density_ = mathsFunc::square(density_);
    P.momentum_ = Vec3D::square(momentum_);
    P.ddensity_ = Vec3D::square(ddensity_);
    P.dmomentum_ = Matrix3D::square(dmomentum_);
    return P;
}

/*!
 * \param[in] P the CGField that has to be copied
 * \return the CGField into which the values are copied
 */
GradVelocityField& GradVelocityField::operator=(const GradVelocityField& P)
= default;

/*!
* \param[in] P the CGField that has to be added
* \return the CGField to which the values are added
*/
GradVelocityField& GradVelocityField::operator+=(const GradVelocityField& P)
{
    density_ += P.density_;
    momentum_ += P.momentum_;
    ddensity_ += P.ddensity_;
    dmomentum_ += P.dmomentum_;
    return *this;
}

/*!
 * \param[in] P the CGField that has to be subtracted
 * \return the CGField from which the values are subtracted
 */
GradVelocityField& GradVelocityField::operator-=(const GradVelocityField& P)
{
    density_ -= P.density_;
    momentum_ -= P.momentum_;
    ddensity_ -= P.ddensity_;
    dmomentum_ -= P.dmomentum_;
    return *this;
}

/*!
 * \param[in] a the scalar that we multiply with
 * \return the CGField  to which the multiplied values are written
 */
GradVelocityField GradVelocityField::operator*(const Mdouble a) const
{
    GradVelocityField p;
    p.density_ = density_ * a;
    p.momentum_ = momentum_ * a;
    p.ddensity_ = ddensity_ * a;
    p.dmomentum_ = dmomentum_ * a;
    return p;
}

/*!
 * \param[in] a the scalar that we divide by
 * \return the CGField to which the divided values are written
 */
GradVelocityField& GradVelocityField::operator/=(const Mdouble a)
{
    density_ /= a;
    momentum_ /= a;
    ddensity_ /= a;
    dmomentum_ /= a;
    return *this;
}

/*!
 * \param[in] phi the value of the cg function at the current CGPoint
 * \param[in] p the particle which is used in the cg function
 */
void GradVelocityField::addParticleStatistics(Mdouble phi, const GradVelocityField& currentInteraction)
{
    density_ += currentInteraction.getDensity() * phi;
    momentum_ += currentInteraction.getMomentum() * phi;
}

/*!
 * \param[in] psi the value of the line integral from C to P at the current CGPoint
 * \param[in] c the contact which is used in the line integral
 */
void GradVelocityField::addInteractionStatistics(Mdouble psi, const GradVelocityField& currentInteraction)
{
}

/*!
 * \param[in] phi the value of the cg function for the contact point of c and 
 * the current CGPoint
 * \param[in] c the interaction which is used in the cg function
 */
void GradVelocityField::addContactPointStatistics(Mdouble phi, const GradVelocityField& currentInteraction)
{
}

void GradVelocityField::addParticleDifferentialStatistics(Vec3D& dphi, const GradVelocityField& currentInteraction)
{
    ddensity_ += currentInteraction.getDensity() * dphi;
    dmomentum_ += Matrix3D::dyadic(currentInteraction.getMomentum(), dphi);
}

/*!
 * \details If the functions returns false, addInteractionStatistics and 
 * addContactPointStatistics are never called.
 * \return True if the class contains fields that are defined as a sum over all 
 * Interactions, else false.
 */
bool GradVelocityField::doInteractionStatistics()
{
    return false;
}

void GradVelocityField::setFields(const BaseInteraction& c, IntegralType type)
{
}

void GradVelocityField::setFields(const BaseParticle& p)
{
    density_ = p.getMass();
    momentum_ = p.getVelocity() * p.getMass();
}

void GradVelocityField::setCylindricalFields(const BaseInteraction& c, IntegralType type)
{
}

void GradVelocityField::setCylindricalFields(const BaseParticle& p)
{
    setFields(p);
    momentum_ = momentum_.getCylindricalTensorField(p.getPosition());
}
    
}
