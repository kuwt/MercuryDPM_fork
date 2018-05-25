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
    const Mdouble ra = p.getWallInteractionRadius(); // note, not getRadius()
    
    // distance between x0 and the *surface* (not the axis)
    distance = sqrt(
            pow((p.getPosition() - position_).getLength(), 2)
            - pow(Vec3D::dot(p.getPosition() - position_, axis_), 2)
    ) - radius_;
    if (distance >= p.getWallInteractionRadius())
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

std::vector<BaseInteraction*> Combtooth::getInteractionWith(BaseParticle* p,
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
        return {c};
    }
    else
        return {};
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
