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

#include "R.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

void R::writeNames(std::ostream& os)
{
    os << "r ";
}

void R::write(std::ostream& os) const
{
    os << r_ << ' ';
}

Mdouble R::getVolumeOfAveragedDimensions(const Vec3D& min, const Vec3D& max)
{
    return (max.Z - min.Z);//2*constants::pi
}

void R::setR(Mdouble r)
{
    r_ = r;
}

/*!
 * \param[in] p vector whose length should be determined
 * \return length of the vector in the non-averaged directions
 * \todo
 */
Mdouble R::getLength(const Vec3D& p)
{
    return sqrt(p.X * p.X + p.Y * p.Y);
}

Mdouble R::getDistanceSquared(const Vec3D& p) const
{
    return mathsFunc::square(getLength(p) - r_);
}

Mdouble R::getINormal(const BaseInteraction& c, const Vec3D& normal) const
{
    if (Vec3D::dot(c.getContactPoint(), normal) > 0)
    {
        return (getLength(c.getI()->getPosition()) - r_);
    }
    else
    {
        return -(getLength(c.getI()->getPosition()) - r_);
    }
}

Mdouble R::getPNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    if (Vec3D::dot(c.getContactPoint(), normal) > 0)
    {
        return (getLength(c.getP()->getPosition()) - r_);
    }
    else
    {
        return -(getLength(c.getP()->getPosition()) - r_);
    }
}

Mdouble R::getCNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    if (Vec3D::dot(c.getContactPoint(), normal) > 0)
    {
        return (getLength(c.getContactPoint()) - r_);
    }
    else
    {
        return -(getLength(c.getContactPoint()) - r_);
    }
}

/*!
 * \details The volume is computed as
 * \f[volume=\int_0^1\sum_{i=1}^n c_i r^i 2 dr = 2 \sum_{i=1}^n c_i/(i+1) \f]
 * \todo
 */
void R::normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff)
{
    Mdouble volume = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); i++)
        volume += coefficients[i] / static_cast<Mdouble>(i + 1);
    volume *= 2.0 * cutoff;
    for (double& coefficient : coefficients)
        coefficient /= volume;
}

Mdouble R::getWeight()
{
    return r_;
}

Mdouble R::getDomainVolume(const Vec3D& min, const Vec3D& max)
{
    //note, the x-coordinate represents the r-coordinate here
    return constants::pi * (max.X * max.X - min.X * min.X) * (max.Z - min.Z);
}

std::string R::getName()
{
    return "R";
}
