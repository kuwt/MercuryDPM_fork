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

#include "Combtooth.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include "Math/ExtendedMath.h"

Combtooth::Combtooth()
{
    axis_ = Vec3D(0, 0, 1);
    position_ = Vec3D(0, 0, 0);
    radius_ = 1;
}

Combtooth::Combtooth(const Combtooth& other) : BaseWall(other)
{
    axis_ = other.axis_;
    position_ = other.position_;
    radius_ = other.radius_;
    
    /* Normalise axis_ */
    axis_ /= axis_.getLength();
}

Combtooth::~Combtooth() = default;

void Combtooth::set(Vec3D axis, Vec3D position, Mdouble radius)
{
    axis_ = axis / axis.getLength();
    position_ = position;
    radius_ = radius;
}

Combtooth* Combtooth::copy() const
{
    return new Combtooth(*this);
}

bool Combtooth::getDistanceAndNormal(const BaseParticle& p,
                                     Mdouble& distance, Vec3D& normal_return) const
{
    /* define shortcuts */
    const Mdouble x0 = p.getPosition().X;
    const Mdouble y0 = p.getPosition().Y;
    const Mdouble z0 = p.getPosition().Z;
    const Mdouble ra = p.getWallInteractionRadius(this); // note, not getRadius()
    
    // distance between x0 and the *surface* (not the axis)
    distance = sqrt(
            pow((p.getPosition() - position_).getLength(), 2)
            - pow(Vec3D::dot(p.getPosition() - position_, axis_), 2)
    ) - radius_;
    if (distance >= p.getWallInteractionRadius(this))
        return false;
    else
    {
        Vec3D axisContactPoint; // the point on the axis closest to the particle
        axisContactPoint = position_ + Vec3D::dot(p.getPosition() - position_, axis_) * axis_;
        normal_return = (axisContactPoint - p.getPosition()); // inward-pointing normal
        normal_return /= normal_return.getLength();
        return true;
    }
}

BaseInteraction* Combtooth::getInteractionWith(BaseParticle* p,
                                                            unsigned timeStamp, InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    if (getDistanceAndNormal(*p, distance, normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal); // outward-pointing normal to cylinder
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        /// \todo Quick hack JMF2 please clean up with teh new way
        c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        return c;
    }
    else
        return nullptr;
}

void Combtooth::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> axis_
       >> dummy >> position_
       >> dummy >> radius_;
}

void Combtooth::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " axis " << axis_
       << " position " << position_
       << " radius " << radius_;
}

std::string Combtooth::getName() const
{
    return "Combtooth";
}
