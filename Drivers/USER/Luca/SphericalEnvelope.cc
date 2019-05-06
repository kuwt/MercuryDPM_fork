
#include "SphericalEnvelope.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

// Last update 22.03.19
// - creation
// - eventually move this into an elliproidal shell

// Default constructor
SphericalEnvelope::SphericalEnvelope()
{
   start_.setZero();
   r_ = 1.0;

   logger(DEBUG, "SphericalEnvelope() constructor finished.");
}

// Copy constructor
SphericalEnvelope::SphericalEnvelope(const SphericalEnvelope& other)
: BaseWall(other)
{
   start_ = other.start_;
   r_ = other.r_;

   logger(DEBUG, "SphericalEnvelope(const SphericalEnvelope&) copy constructor finished.");
}

// Initialized constructor
// Parameters in: centre of the sphere and its radius
SphericalEnvelope::SphericalEnvelope(Vec3D start, Mdouble radius)
{
   start_ = start;
   r_ = radius;

   logger(DEBUG, "SphericalEnvelope(Vec3D, Mdouble) constructor finished.");
}

// Destructor
SphericalEnvelope::~SphericalEnvelope()
{
   logger(DEBUG, "~SphericalEnvelope() finished, destroyed the screw.");
}

// A pointer to the copy of the envelope
SphericalEnvelope* SphericalEnvelope::copy() const
{
   return new SphericalEnvelope(*this);
}

// Function to set the envelope parameters
void SphericalEnvelope::set(Vec3D start, Mdouble radius)
{
   start_ = start;
   r_ = radius;

   std::cout << "\n\nSpherical envelope parameters set\n";
   std::cout << "Start: " << start_ << "\n";
   std::cout << "Radius: " << r_ << "\n";
}

// Function to set the radius
void SphericalEnvelope::setRadius(Mdouble radius)
{
   r_ = radius;
}

// Function to get the radius
double SphericalEnvelope::getRadius(Mdouble radius)
{
   return r_;
}

// Input stream
void SphericalEnvelope::read(std::istream& is)
{
   BaseWall::read(is);
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> r_;
}

// Old input stream (no idea what this is)
void SphericalEnvelope::oldRead(std::istream& is)
{
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> r_;
}

// Writer
void SphericalEnvelope::write(std::ostream& os) const
{
   BaseWall::write(os);
   os << " Start " << start_
   << " Radius " << r_;
}

// Returns the string "SphericalEnvelope"
std::string SphericalEnvelope::getName() const
{
   return "SphericalEnvelope";
}

// The contact detection algorithm
bool SphericalEnvelope::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
   // radial position of the particle with respect to the envelope centre in spherical coordinates
   Mdouble r_P = (p.getPosition() - start_).getLength();

   // collision check
   if (r_P > r_ - p.getWallInteractionRadius()) // collision
   {
      distance = r_ - r_P;
      // this error check allows for unphysical compressions up to twice the particle radius due to the heavy compression needed
      // if (distance < - p.getWallInteractionRadius()) std::cout << std::endl << "COLLISION ERROR WITH THE SPHERICAL SURFACE" << std::endl;

      normal_return = p.getPosition()/r_P;

      return true;
   }
   else // no collision
   {
      return false;
   }
}

// Checks for the interaction between a particle p at a time timeStamp
// In case of interaction returns a pointer to the BaseInteraction happened between the Envelope and the BaseParticle at time timeStamp
std::vector<BaseInteraction *>
SphericalEnvelope::getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler)
{
   Mdouble distance;
   Vec3D normal;
   std::vector<BaseInteraction*> interactions;
   if (getDistanceAndNormal(*p,distance,normal))
   {
      BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
      c->setNormal(-normal);
      c->setDistance(distance);
      c->setOverlap(p->getRadius() - distance);
      c->setContactPoint(p->getPosition()-(p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
      interactions.push_back(c);
   }
   return interactions;
}
