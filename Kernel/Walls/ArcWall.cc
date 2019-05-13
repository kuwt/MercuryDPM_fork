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
    semiangle_ = 20 * constants::pi / 180.;
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
    axis_ = axis_ / axis_.getLength();
    centreline_ = centreline_ - Vec3D::dot(centreline_, axis_) * axis_;
    centreline_ /= centreline_.getLength();
}

ArcWall* ArcWall::copy() const
{
    return new ArcWall(*this);
}

bool ArcWall::getDistanceAndNormal(const BaseParticle& p,
                                   Mdouble& distance, Vec3D& normal_return) const
{
    // distance between x0 and the *surface* (not the axis)
    distance = radius_ - sqrt(
            pow((p.getPosition() - pos_).getLength(), 2)
            - pow(Vec3D::dot(p.getPosition() - pos_, axis_), 2)
    );
    if (distance >= p.getWallInteractionRadius(this))
        return false;
    else
    {
        Vec3D axisContactPoint; // the point on the axis closest to the particle
        axisContactPoint = pos_ + Vec3D::dot(p.getPosition() - pos_, axis_) * axis_;
        // outward-pointing normal, for an inward-pointing force
        normal_return = p.getPosition() - axisContactPoint;
        normal_return /= normal_return.getLength();
        
        /* There is an interaction iff the interaction normal is within
         * semiangle of the centreline. */
        double angle_between = acos(Vec3D::dot(normal_return, centreline_));
        if (angle_between < semiangle_ || semiangle_ >= constants::pi)
            return true;
        else
            return false;
    }
}

BaseInteraction* ArcWall::getInteractionWith(BaseParticle* p,
                                                          unsigned timeStamp, InteractionHandler* interactionHandler)
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
