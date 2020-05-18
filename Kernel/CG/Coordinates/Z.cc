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

#include "Z.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

void Z::writeNames(std::ostream& os)
{
    os << "z ";
}

void Z::write(std::ostream& os) const
{
    os << z_ << ' ';
}

///\todo Generalise to 2D

Mdouble Z::getVolumeOfAveragedDimensions(const Vec3D& min, const Vec3D& max)
{
    return (max.X - min.X) * (max.Y - min.Y);
}

void Z::setZ(Mdouble z)
{
    z_ = z;
}

Mdouble Z::getDistanceSquared(const Vec3D& p) const
{
    return mathsFunc::square(p.Z - z_);
}

/*!
 * \param[in] p vector whose length should be determined
 * \return length of the vector in the non-averaged directions
 * \todo
 */
Mdouble Z::getLength(const Vec3D& p)
{
    return fabs(p.Z);
}

/*!
 * For all points S on the contact line from I to P, this is the minimum value for (S-R)*normal.
 * This is the lower limit of integration along the contact line
 */
Mdouble Z::getINormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getI()->getPosition().Z - z_) * normal.Z;
}

/*!
 * For all points S on the contact line from I to P, this is the maximum value for (S-R)*normal.
 * This is the upper limit of integration along the contact line
 */
Mdouble Z::getPNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getP()->getPosition().Z - z_) * normal.Z;
}

/*!
 * For all points S on the contact line from I to P, this is the maximum value for (S-R)*normal.
 * This is the upper limit of integration along the contact line
 */
Mdouble Z::getCNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getContactPoint().Z - z_) * normal.Z;
}

std::string Z::getName()
{
    return "Z";
}

