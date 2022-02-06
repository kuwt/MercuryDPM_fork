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

#include "ArcWall.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include "Math/ExtendedMath.h"

ArcWall::ArcWall()
{
    axis_ = Vec3D(0, 0, 1);
    pos_ = Vec3D(0, 1, 0);
    radius_ = 1;
    centreline_ = Vec3D(0, -1, 0);
    semiangle_ = 20. * constants::degree;
}

ArcWall::ArcWall(const ArcWall& other) : BaseWall(other)
{
    axis_ = other.axis_;
    pos_ = other.pos_;
    radius_ = other.radius_;
    centreline_ = other.centreline_;
    semiangle_ = other.semiangle_;
}

void ArcWall::set(Vec3D axis, Vec3D pos, Mdouble radius, Vec3D centreline, Mdouble semiangle)
{
    axis_ = axis;
    pos_ = pos;
    radius_ = radius;
    centreline_ = centreline;
    semiangle_ = semiangle;
    
    /* Normalise axis_ and centreline_, and make centreline_ perpendicular to axis_ */
    axis_.normalise();
    centreline_ = centreline_ - Vec3D::dot(centreline_, axis_) * axis_;
    centreline_.normalise();
}

ArcWall* ArcWall::copy() const
{
    return new ArcWall(*this);
}

bool ArcWall::getDistanceAndNormal(const BaseParticle& p,
                                   Mdouble& distance,
                                   Vec3D& normal) const
{
    Vec3D nearestPointOnAxis = pos_ + Vec3D::dot(p.getPosition() - pos_, axis_) * axis_;
    Mdouble distanceFromAxis = (p.getPosition() - nearestPointOnAxis).getLength();
    // distance between x0 and the surface
    distance = radius_ - distanceFromAxis;
    if (distance >= p.getWallInteractionRadius(this))
        return false;
    
    // the normal points into the wall, i.e. out from the axis
    normal = p.getPosition() - nearestPointOnAxis;
    normal.normalise();
    
    // If the arc is actually a full circle then there is definitely an
    // interaction, regardless of direction.
    if (semiangle_ >= 180. * constants::degree)
        return true;
    
    // Otherwise it depends on whether the interaction direction is
    // sufficiently aligned with the centreline.
    return acos(Vec3D::dot(normal, centreline_)) < semiangle_;
}

BaseInteraction* ArcWall::getInteractionWith(BaseParticle* p,
                                             unsigned timeStamp,
                                             InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    if (getDistanceAndNormal(*p, distance, normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        /// \todo Hacked please fix @Thomas
        return c;
    }
    else
        return nullptr;
}

void ArcWall::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> axis_
       >> dummy >> pos_
       >> dummy >> radius_
       >> dummy >> centreline_
       >> dummy >> semiangle_;
}

void ArcWall::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " axis " << axis_
       << " pos " << pos_
       << " radius " << radius_
       << " centreline " << centreline_
       << " semiangle " << semiangle_;
}

std::string ArcWall::getName() const
{
    return "ArcWall";
}
