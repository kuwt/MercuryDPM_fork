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

#include "XYZ.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

void XYZ::writeNames(std::ostream& os)
{
    os << "x y z ";
}

/*!
 * \param[out] os the ostream file to which the position is written
 */
void XYZ::write(std::ostream& os) const
{
    os << p_ << ' ';
}

/*!
 * \details If averaged dimensions are present (i.e. for all Coordinates except 
 * XYZ), the CGFunction has to be divided by a factor (here called volume)
 * due to integrating the variables over the averaged dimensions.
 * \param[in] min the lower limits of the mesh domain (xMin, yMin, zMin) 
 * \param[in] max the upper limits of the mesh domain (xMax, yMax, zMax) 
 * \return the volume factor
 */
Mdouble XYZ::getVolumeOfAveragedDimensions(const Vec3D& min UNUSED, const Vec3D& max UNUSED)
{
    return 1.0;
}

/*!
 * \param[in] p the position that is to be set.
 */
void XYZ::setXYZ(Vec3D p)
{
    p_ = p;
}

/*!
 * \details This function is needed to evaluate the CGFunction, as this function 
 * has the distance between the CGPoint and the Particle as an argument. To 
 * properly account for the averaging, the distance is only computed in the 
 * non-averaged directions.
 * \param[in] p the position of a particle for which the distance is computed.
 */
Mdouble XYZ::getDistanceSquared(const Vec3D& p) const
{
    return Vec3D::getLengthSquared(p_ - p);
}

/*!
 * \param[in] p vector whose length should be determined
 * \return length of the vector in the non-averaged directions
 * \todo
 */
Mdouble XYZ::getLength(const Vec3D& p)
{
    return sqrt(p.X * p.X + p.Y * p.Y + p.Z * p.Z);
}

/*!
 * \param[in] c the Interaction object from which iNormal is computed
 * \return iNormal, one of the three distances needed to calculate the line 
 * integral \f$\psi\f$ which defines the stress distribution (see image).
 * \image html LineIntegral.jpeg Illustration of the line integral
 */
Mdouble XYZ::getINormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return Vec3D::dot(c.getI()->getPosition() - p_, c.getNormal());
}

/*!
* \param[in] c the Interaction object from which iNormal is computed
* \return pNormal, one of the three distances needed to calculate the line
* integral \f$\psi\f$ which defines the stress distribution (see image).
* \image html LineIntegral.jpeg Illustration of the line integral
*/Mdouble XYZ::getPNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return Vec3D::dot(c.getP()->getPosition() - p_, c.getNormal());
}

/*!
* \param[in] c the Interaction object from which iNormal is computed
* \return cNormal, one of the three distances needed to calculate the line
* integral \f$\psi\f$ which defines the stress distribution (see image).
* \image html LineIntegral.jpeg Illustration of the line integral
*/Mdouble XYZ::getCNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return Vec3D::dot(c.getContactPoint() - p_, c.getNormal());
}

/*!
 * \param[in] c the Interaction object from which iNormal is computed
 * \param[in] pNormal the output of getPNormal needed for the computation.
 * \return iNormal, one of the three distances needed to calculate the line 
 * integral \f$\psi\f$ which defines the stress distribution (see image).
 * \image html LineIntegral.jpeg Illustration of the line integral
 */
Mdouble XYZ::getTangentialSquared(const BaseInteraction& c, Mdouble pNormal) const
{
    return Vec3D::getLengthSquared(c.getP()->getPosition() - p_) - mathsFunc::square(pNormal);
}

/*!
 * \details The prefactor of the Gauss CGFunction is set such that the integral 
 * over the non-averaged dimensions is unity.
 * \param[in] width width (equals the standard deviation in 1D) of the Gauss CGFunction.
 * \param[in] cutoff cutoff of the Gauss CGFunction
 * \return the prefactor of the Gauss CGFunction.
 */
Mdouble XYZ::getGaussPrefactor(Mdouble width, Mdouble cutoff)
{
    //Wolfram alpha: erf(c/(sqrt(2) w))-(sqrt(2/pi) c e^(-c^2/(2 w^2)))/w
    Mdouble prefactor = 1.0 / (constants::sqrt_2 * constants::sqrt_pi * width);
    Mdouble cw = cutoff / width;
    return mathsFunc::cubic(prefactor) / (
            erf(cw / constants::sqrt_2)
            - constants::sqrt_2 / constants::sqrt_pi * cw * exp(-0.5 * mathsFunc::square(cw))
    );
}

/*!
 * \details The prefactor of the Gauss line integral is set such that the integral 
 * over the non-averaged dimensions is unity.
 * \param[in] distance length of the branch vector along which the line integral is evaluated.
 * \param[in] width width (equals the standard deviation in 1D) of the Gauss CGFunction.
 * \param[in] cutoff cutoff of the Gauss CGFunction
 * \return the prefactor of the Gauss CGFunction.
 */
Mdouble XYZ::getGaussIntegralPrefactor(Mdouble distance, Mdouble width, Mdouble cutoff)
{
    Mdouble widthSqrt2 = width * constants::sqrt_2;
    Mdouble a = -cutoff;
    Mdouble b = cutoff + distance;
    //full 2D prefactor
    Mdouble prefactor = 1.0 / (constants::sqrt_2 * constants::sqrt_pi * width);
    prefactor = mathsFunc::square(prefactor) / (1.0 - exp(-0.5 * mathsFunc::square(cutoff / width)));
    return prefactor * 0.5 / (
            +erf(b / widthSqrt2) * b
            + widthSqrt2 / constants::sqrt_pi * exp(-mathsFunc::square(b / widthSqrt2))
            - erf(a / widthSqrt2) * a
            - widthSqrt2 / constants::sqrt_pi * exp(-mathsFunc::square(a / widthSqrt2))
    );
}

/*!
 * \details The volume is computed as
 * \f[volume= \int_0^1\sum_{i=1}^n c_i r/c^i 4 pi r^2 dr = 4 pi \sum_{i=1}^n c_i/(i+3) \f]
 * with 4 pi r^2 the surface area of a sphere.
 * \param[in,out] coefficients the coefficients of Polynomial CGFunctions.
 * \param[in] cutoff cutoff of the Gauss CGFunction
 */
void XYZ::normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff)
{
    Mdouble volume = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); i++)
        volume += coefficients[i] / static_cast<Mdouble>(i + 3);
    volume *= 4.0 * constants::pi * mathsFunc::cubic(cutoff);
    for (double& coefficient : coefficients)
        coefficient /= volume;
}

const unsigned XYZ::countVariables()
{
    return 3;
}

std::string XYZ::getName()
{
    return "XYZ";
}

