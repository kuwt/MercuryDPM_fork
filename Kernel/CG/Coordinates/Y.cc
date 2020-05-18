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

#include "Y.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

void Y::writeNames(std::ostream& os)
{
    os << "y ";
}

void Y::write(std::ostream& os) const
{
    os << y_ << ' ';
}

///\todo Generalise to 2D
Mdouble Y::getVolumeOfAveragedDimensions(const Vec3D& min, const Vec3D& max)
{
    return (max.X - min.X) * (max.Z - min.Z);
}

void Y::setY(Mdouble y)
{
    y_ = y;
}

Mdouble Y::getDistanceSquared(const Vec3D& p) const
{
    return mathsFunc::square(p.Y - y_);
}

/*!
 * \param[in] p vector whose length should be determined
 * \return length of the vector in the non-averaged directions
 * \todo
 */
Mdouble Y::getLength(const Vec3D& p)
{
    return fabs(p.Y);
}

Mdouble Y::getINormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getI()->getPosition().Y - y_) * normal.Y;
}

Mdouble Y::getPNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getP()->getPosition().Y - y_) * normal.Y;
}

Mdouble Y::getCNormal(const BaseInteraction& c, const Vec3D& normal) const
{
    return (c.getContactPoint().Y - y_) * normal.Y;
}

std::string Y::getName()
{
    return "Y";
}
