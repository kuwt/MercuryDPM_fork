#include "CompressionPistonSurface.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// This is the texture of the surface of the compression piston

// Default constructor
CompressionPistonSurface::CompressionPistonSurface()
{
   position_ = 1.0;
   height_ = 1.0;
   thickness_ = 1.0;
   angle_ = 0.0;

   logger(DEBUG, "CompressionPistonSurface() constructor finished.");
}

// Copy constructor
CompressionPistonSurface::CompressionPistonSurface(const CompressionPistonSurface& other)
: BaseWall(other)
{
   position_ = other.position_;
   height_ = other.height_;
   thickness_ = other.thickness_;
   angle_ = other.angle_;

   logger(DEBUG, "CompressionPistonSurface(const CompressionPistonSurface&) copy constructor finished.");
}

// The good constructor
CompressionPistonSurface::CompressionPistonSurface(Mdouble position, Mdouble height, Mdouble thickness, Mdouble angle)
{
   position_ = position;
   height_ = height;
   thickness_ = thickness;
   angle_ = angle;

   logger(DEBUG, "CompressionPistonSurface(Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
CompressionPistonSurface::~CompressionPistonSurface()
{
   logger(DEBUG, "~CompressionPistonSurface() finished, destroyed the Compression Piston.");
}

// A pointer to the copy of the piston
CompressionPistonSurface* CompressionPistonSurface::copy() const
{
   return new CompressionPistonSurface(*this);
}

// Function to set the piston parameters
void CompressionPistonSurface::set(Mdouble position, Mdouble height, Mdouble thickness, Mdouble angle)
{
   position_ = position;
   height_ = height;
   thickness_ = thickness;
   angle_ = angle;

   std::cout << "\nCompression piston parameters set.\n";
   std::cout << "Initial position: " << position_ << "\n";
   std::cout << "Height: " << height_ << "\n";
   std::cout << "Thickness: " << thickness_ << "\n";
   std::cout << "Angle: " << angle_ << "\n";
}

// Moves the piston
void CompressionPistonSurface::move(Mdouble dh)
{
   position_ += dh;
}

// Rotates the piston
void CompressionPistonSurface::rotate(Mdouble dTheta)
{
   angle_ += dTheta;
}

// Returns the actual piston position
double CompressionPistonSurface::getPosition() const
{
   return position_;
}

// Returns the actual piston angle
double CompressionPistonSurface::getAngle() const
{
   return angle_;
}

// Reads a compressing piston surface from an input stream, for example a restart file.
void CompressionPistonSurface::read(std::istream& is)
{
   BaseWall::read(is);
   std::string dummy;
   is >> dummy >> position_
   >> dummy >> height_
   >> dummy >> thickness_
   >> dummy >> angle_;
}

// Writes this compressing piston surface to an output stream, for example a restart file.
void CompressionPistonSurface::write(std::ostream& os) const
{
   BaseWall::write(os);
   os << " Position " << position_
   << " Height " << height_
   << " Thickness " << thickness_
   << " Angle " << angle_;
}

// Returns the name of the driver
std::string CompressionPistonSurface::getName() const
{
   return "CompressionPistonSurface";
}

// The contact detection algorithm
// Takes as input a pointer to a particle and evaluates if there is a collision.
// If so updates the distance between the particle centre and the collision point, the collision normal and returns true.
// Returns false if there is no collision.
bool CompressionPistonSurface::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
   // Checks for a possible collision along the Z direction
   if (p.getPosition().Z + p.getRadius() < position_ - height_) return false;

   // Checks for a possible collision along the XY plane
   if (fabs(p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_)) > p.getRadius() + thickness_) return false;

   // Collision cases
   if (p.getPosition().Z > position_ - height_) // Pure radial collision
   {
      distance = fabs(p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_)) - thickness_;
      normal_return.X = -sin(angle_);
      normal_return.Y = cos(angle_);
      normal_return.Z = 0.0;

      // reverse the normal direction if the particle is "on the other (radial) side" of the wall
      if (p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_) > 0.0)
      {
         normal_return.X *= -1.0;
         normal_return.Y *= -1.0;
      }

      return true;
   }
   else if (fabs(p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_)) < thickness_) // Pure axial collision
   {
      distance = position_ - height_ - p.getPosition().Z;
      normal_return.X = 0.0;
      normal_return.Y = 0.0;
      normal_return.Z = -1.0;

      return true;
   }
   else // Mixed collision
   {
      distance = sqrt(pow(p.getPosition().Z - position_ + height_, 2.0) + pow(fabs(p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_)) - thickness_, 2.0));

      // check if there is collision with the blade's edge
      if (distance > p.getRadius()) return false;

      normal_return.X = -sin(angle_)*(fabs(p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_)) - thickness_);
      normal_return.Y = cos(angle_)*(fabs(p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_)) - thickness_);
      normal_return.Z = -1.0*(p.getPosition().Z - position_ + height_);
      normal_return /= normal_return.getLength();

      // reverse the normal direction if the particle is "on the other (radial) side" of the wall
      if (p.getPosition().Y*cos(angle_) - p.getPosition().X*sin(angle_) > 0.0)
      {
         normal_return.X *= -1.0;
         normal_return.Y *= -1.0;
      }

      return true;
   }

   std::cout << "Collison failure: the function did not return properly" << std::endl;
   return false;

}

// Checks for the interaction between a particle p at a time timeStamp.
// In case of interaction returns a pointer to the BaseInteraction happened between the inter-well wall and the
// BaseParticle at time timeStamp
std::vector<BaseInteraction *> CompressionPistonSurface::getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler)
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
      c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
      interactions.push_back(c);
   }
   return interactions;
}
