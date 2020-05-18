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

#include "XZ.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

void XZ::writeNames(std::ostream& os)
{
    os << "x z ";
}

void XZ::write(std::ostream& os) const
{
    os << x_ << ' ' << z_ << ' ';
}

///\todo Generalise to 2D
Mdouble XZ::getVolumeOfAveragedDimensions(const Vec3D& min, const Vec3D& max)
{
    return max.Y - min.Y;
}

void XZ::setXZ(Mdouble x, Mdouble z)
{
    x_ = x;
    z_ = z;
}

Mdouble XZ::getDistanceSquared(const Vec3D& p) const
{
    return mathsFunc::square(p.X - x_) + mathsFunc::square(p.Z - z_);
}

/*!
 * \param[in] p vector whose length should be determined
 * \return length of the vector in the non-averaged directions
 * \todo
 */
Mdouble XZ::getLength(const Vec3D& p)
{
    return sqrt(p.X * p.X + p.Z * p.Z);
}

Mdouble XZ::getINormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getI()->getPosition().X - x_) * normal.X
           + (c.getI()->getPosition().Z - z_) * normal.Z;
}

Mdouble XZ::getPNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getP()->getPosition().X - x_) * normal.X
           + (c.getP()->getPosition().Z - z_) * normal.Z;
}

Mdouble XZ::getCNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getContactPoint().X - x_) * normal.X
           + (c.getContactPoint().Z - z_) * normal.Z;
}

Mdouble XZ::getTangentialSquared(const BaseInteraction& c, Mdouble pNormal) const
{
    return mathsFunc::square(c.getP()->getPosition().X - x_)
           + mathsFunc::square(c.getP()->getPosition().Z - z_)
           - mathsFunc::square(pNormal);
}

std::string XZ::getName()
{
    return "XZ";
}
