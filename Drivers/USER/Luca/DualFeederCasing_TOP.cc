
#include "DualFeederCasing_TOP.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

// Last update 19.03.19
// - creation

// Default constructor
DualFeederCasing_TOP::DualFeederCasing_TOP()
{
   start_.setZero();
   l_ = 1.0;
   r_ = 1.0;
   delta_ = 0.0;

   logger(DEBUG, "DualFeederCasing_TOP() constructor finished.");
}

// Copy constructor
DualFeederCasing_TOP::DualFeederCasing_TOP(const DualFeederCasing_TOP& other)
: BaseWall(other)
{
   start_ = other.start_;
   l_ = other.l_;
   r_ = other.r_;
   delta_ = other.delta_;

   logger(DEBUG, "DualFeederCasing_TOP(const DualFeederCasing_TOP&) copy constructor finished.");
}

// Initialized constructor
// Parameters in: initial starting point of the axis, length, radius, screw axes distance.
DualFeederCasing_TOP::DualFeederCasing_TOP(Vec3D start, Mdouble length, Mdouble radius, Mdouble distance)
{
   start_ = start;
   l_ = length;
   r_ = radius;
   delta_ = distance;

   logger(DEBUG, "DualFeederCasing_TOP(Vec3D, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
DualFeederCasing_TOP::~DualFeederCasing_TOP()
{
   logger(DEBUG, "~DualFeederCasing_TOP() finished, destroyed the screw.");
}

// A pointer to the copy of the casing
DualFeederCasing_TOP* DualFeederCasing_TOP::copy() const
{
   return new DualFeederCasing_TOP(*this);
}

// Function to set the casing parameters
void DualFeederCasing_TOP::set(Vec3D start, Mdouble length, Mdouble radius, Mdouble distance)
{
   start_ = start;
   l_ = length;
   r_ = radius;
   delta_ = distance;

   std::cout << "\n\nTOP casing crew parameters set\n";
   std::cout << "Start: " << start_ << "\n";
   std::cout << "Length: " << l_ << "\n";
   std::cout << "Radius: " << r_ << "\n";
   std::cout << "Distance: " << delta_ << "\n";
}

// Function to set the radius
void DualFeederCasing_TOP::setRadius(Mdouble radius)
{
   r_ = radius;
}

// Function to set the distance
void DualFeederCasing_TOP::setDistance(Mdouble distance)
{
   delta_ = distance;
}

// Input stream
void DualFeederCasing_TOP::read(std::istream& is)
{
   BaseWall::read(is);
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> l_
   >> dummy >> r_
   >> dummy >> delta_;
}

// Old input stream (no idea what this is)
void DualFeederCasing_TOP::oldRead(std::istream& is)
{
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> l_
   >> dummy >> r_
   >> dummy >> delta_;
}

// Writer
void DualFeederCasing_TOP::write(std::ostream& os) const
{
   BaseWall::write(os);
   os << " Start " << start_
   << " Length " << l_
   << " Radius " << r_
   << " Distance " << delta_;
}

// Returns the string "DualFeederCasing_TOP"
std::string DualFeederCasing_TOP::getName() const
{
   return "DualFeederCasing_TOP";
}

// The contact detection algorithm
bool DualFeederCasing_TOP::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
   // the axes of the screws are parallel to the z axis, and lay in the plane y = 0.0
   // the screw axes are located at points c = (+-delta_, 0.0, z), so the only important distance is the radial from +-c
   // if the particle lays below the screw axes there is no collision
   if (p.getPosition().Y < start_.Y) return false;

   // if the particle is outside the length containing the screw there is no collision
   if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius()) return false;
   if (p.getPosition().Z < start_.Z - p.getWallInteractionRadius()) return false;

   // radial position of the particle with respect to the closest screw axis
   Mdouble r_P = sqrt(pow(fabs(p.getPosition().X - start_.X) - delta_, 2.0) + pow(p.getPosition().Y - start_.Y, 2.0));

   // trigonometric functions relative to the particle angle (evaluated with respect to the closest sxcrew axis)
   Mdouble cosTheta_P = (fabs(p.getPosition().X - start_.X) - delta_)/r_P;
   Mdouble sinTheta_P = (p.getPosition().Y - start_.Y)/r_P;

   // getting the sign of the X component of the normal right
   if (p.getPosition().X - start_.X < 0.0) cosTheta_P *= -1.0;

   // radial vector at the particle position (with respect to the closest screw axis)
   Vec3D radialVector;
   radialVector.X = cosTheta_P;
   radialVector.Y = sinTheta_P;
   radialVector.Z = 0.0;

   // distinguishing between the cases
   if (fabs(p.getPosition().X - start_.X) < p.getWallInteractionRadius()) // collision with the edgy part where the cylinders merge together
   {
      // Y position of the merging point (the X is zero by definition)
      Mdouble m_Y = sqrt(pow(r_, 2.0) - pow(delta_, 2.0));

      // distance from the merging point
      distance = sqrt(pow(p.getPosition().X - start_.X, 2.0) + pow(p.getPosition().Y - start_.Y - m_Y, 2.0));

      // if the distance is negative prints an error message
      if (distance < 0.0) std::cout << std::endl << "COLLISION ERROR WITH THE CYLINDER MERGING POINT" << std::endl;
      if (distance > p.getWallInteractionRadius()) return false;

      // Y vector
      Vec3D yVector;
      yVector.X = 0.0;
      yVector.Y = 1.0;
      yVector.Z = 0.0;

      // the normal to the collision is a linear combination of radial and Y normals, weighted via the relative position of the particle with respect to the merging point
      normal_return = fabs(p.getPosition().Y - start_.Y - m_Y)*yVector + fabs(p.getPosition().X - start_.X)*radialVector;
      normal_return /= normal_return.getLength();

      return true;
   }
   // else if (fabs(p.getPosition().X - start_.X) > r_ + delta_ - p.getWallInteractionRadius()) // collision with the flat parts where the cylinders merge with the vertical walls
   // {
   //    // X position of the merging point (the Y is zero by definition)
   //    Mdouble n_X = r_ + delta_;
   //
   //    // distance from the merging point
   //    distance = sqrt(pow(fabs(p.getPosition().X - start_.X) - n_X, 2.0) + pow(p.getPosition().Y - start_.Y, 2.0));
   // 
   //    // if the distance is negative prints an error message
   //    if (distance < 0.0) std::cout << std::endl << "COLLISION ERROR WITH THE WALLS MERGING POINT" << std::endl;
   //    if (distance > p.getWallInteractionRadius()) return false;
   //
   //    // X vector
   //    Vec3D xVector;
   //    xVector.X = 1.0;
   //    xVector.Y = 0.0;
   //    xVector.Z = 0.0;
   //    if (p.getPosition().X - start_.X < 0.0) xVector.X *= -1.0;
   //
   //    // the normal to the collision is a linear combination of radial and X normals, weighted via the relative position of the particle with respect to the merging point
   //    normal_return = fabs(p.getPosition().Y - start_.Y)*xVector + fabs(fabs(p.getPosition().X - start_.X) - n_X)*radialVector;
   //    normal_return /= normal_return.getLength();
   //
   //    return true;
   // }
   else // collision with the internal surface of the cylinders
   {
      // distance from the closest cylindrical surface based on closest screw axis
      distance = r_ - sqrt(pow(fabs(p.getPosition().X - start_.X) - delta_, 2.0) + pow(p.getPosition().Y - start_.Y, 2.0));

      // if the distance is negative prints an error message
      if (distance < 0.0) std::cout << std::endl << "COLLISION ERROR WITH THE CYLINDRICAL SURFACE" << std::endl;
      if (distance > p.getWallInteractionRadius()) return false;

      // the normal is purely radial
      normal_return = radialVector;
      normal_return /= normal_return.getLength();

      return true;
   }

   // error message if some case has not been accounted for
   std::cout << "\nERROR: THE SCREW COLLISION FUNCTION DID NOT RETURN PROPERLY\n";

   return false;
}

// Checks for the interaction between a particle p at a time timeStamp
// In case of interaction returns a pointer to the BaseInteraction happened between the Screw and the BaseParticle at time timeStamp
std::vector<BaseInteraction *>
DualFeederCasing_TOP::getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler)
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
