#include "VChute.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include "Math/ExtendedMath.h"

/* A finite-length V-chute. */

VChute::VChute()
{
    l_ = 1.0;
    w_ = 1.0;
    alpha_ = 20.0 / 180.0 * constants::pi;
}

VChute::VChute(const VChute& other) : BaseWall(other)
{
    l_ = other.l_;
    w_ = other.w_;
    alpha_ = other.alpha_;
}

VChute::VChute(Mdouble length, Mdouble width, Mdouble alpha)
{
    l_ = length;
    w_ = width;
    alpha_ = alpha;
}

VChute::~VChute() = default;

void VChute::set(Mdouble length, Mdouble width, Mdouble alpha)
{
    l_ = length;
    w_ = width;
    alpha_ = alpha;
}

VChute* VChute::copy() const
{
    return new VChute(*this);
}

bool VChute::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    /* define shortcuts */
    const Mdouble x0 = p.getPosition().X;
    const Mdouble y0 = p.getPosition().Y;
    const Mdouble z0 = p.getPosition().Z;
    const Mdouble ra = p.getWallInteractionRadius(this); // note, not getRadius()
    
    /* Has the particle flown off the ends of the chute? */
    if (x0 < -ra || x0 > l_ + ra)
        return false;
    
    /* Don't bother bounding in the y, z directions. The wall is linear, so
     * calculating collisions isn't that expensive anyway. */
    /*
    if (fabs(y0) > w_ + ra
     || z0 < -ra || z0 > tan(alpha_)*(w_ + ra) + ra ) {
        return false;
    }
    */
    
    /* Get the distance and normal. */
    // Newton is unnecessary, because linear problem.
    // The special case y0==0 is problematic because it means double-contact, but
    // let's just hope it doesn't come up. (TODO)
    Mdouble q;
    Mdouble distanceSquared;
    if (y0 > 0)
    {
        q = y0 * pow(cos(alpha_), 2) + z0 * tan(alpha_);
        distanceSquared = pow(
                y0 * (pow(cos(alpha_), 2) - 1) + z0 * sin(alpha_) * cos(alpha_), 2)
                          + pow(
                y0 * sin(alpha_) * cos(alpha_) + z0 * (pow(sin(alpha_), 2) - 1), 2);
    }
    else if (y0 < 0)
    {
        q = y0 * pow(cos(alpha_), 2) - z0 * tan(alpha_);
        distanceSquared = pow(
                y0 * (pow(cos(alpha_), 2) - 1) - z0 * sin(alpha_) * cos(alpha_), 2)
                          + pow(
                -y0 * sin(alpha_) * cos(alpha_) + z0 * (pow(sin(alpha_), 2) - 1), 2);
    }
    else
    {
        /* The case y=0 requires special handling... */
        
        /* In the past, this case simply triggered an error and died. */
        // return false; 
        
        /* If the particle has y=0 then it is in contact with both sides of the
         * V-chute; both sides will exert a force on the particle
         * simultaneously. This is difficult to handle. I shall assume that the
         * resultant force is exactly upwards, as though the particle was in
         * contact with a flat plane (z=0).
         * TODO try to make this more realistic
         */
        std::cerr << "y0 == 0, dangerous" << std::endl;
        q = 0;
        distanceSquared = pow(z0, 2);
    }
    
    //If distance is too large there is no contact
    if (distanceSquared >= pow(p.getWallInteractionRadius(this), 2))
    {
        return false;
    }
    
    Vec3D ContactPoint;
    distance = sqrt(distanceSquared);
    ContactPoint.X = x0;
    ContactPoint.Y = q;
    ContactPoint.Z = fabs(q) * tan(alpha_);
    normal_return = ContactPoint - p.getPosition();
    normal_return /= normal_return.getLength();
    return true;
    
}

BaseInteraction*
VChute::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
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
        /// \todo Hacked please fix
        return c;
    }
    else
        return nullptr;
}

/*!
 * \param[in,out] is The input stream from which the Coil is read.
 */
void VChute::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> l_
       >> dummy >> w_
       >> dummy >> alpha_;
}

/*!
 * \param[in,out] os The outpus stream to which the Coil is written.
 */
void VChute::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " Length " << l_
       << " Width " << w_
       << " alpha " << alpha_;
}

/*!
 * \return The string "VChute".
 */
std::string VChute::getName() const
{
    return "VChute";
}
