
#include "ScrewCC.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

// Last update --.--.--

/*

*/

// Default constructor
ScrewCC::ScrewCC()
{
   start_.setZero();
   l_ = 1.0;
   rMax_ = 1.0;
   rMin_ = 1.0;
   n_ = 1.0;
   omega_ = 1.0;
   thickness_ = 0.0;
   delta_ = 0.0;
   pitch_ = 1.0;
   h_ = 0.5/constants::pi;
   offset_ = 0.0;

   logger(DEBUG, "ScrewCC() constructor finished.");
}

// Copy constructor
ScrewCC::ScrewCC(const ScrewCC& other)
: BaseWall(other)
{
   start_ = other.start_;
   l_ = other.l_;
   rMax_ = other.rMax_;
   rMin_ = other.rMax_;
   n_ = other.n_;
   omega_ = other.omega_;
   thickness_ = other.thickness_;
   delta_ = 0.5*other.thickness_;
   pitch_ = other.l_/other.n_;
   h_ = 0.5*other.l_/(constants::pi * other.n_);
   offset_ = other.offset_;

   logger(DEBUG, "ScrewCC(const ScrewCC&) copy constructor finished.");
}

// Initialized constructor
// Parameters in: initial starting point of the screw, length, blade radius, number of turns, angular velocity, thickness of the blade.
ScrewCC::ScrewCC(Vec3D start, Mdouble l, Mdouble r1, Mdouble r0, Mdouble n, Mdouble omega, Mdouble thickness)
{
   start_ = start;
   l_ = l;
   rMax_ = r1;
   rMin_ = r0;
   n_ = n;
   omega_ = omega;
   thickness_ = thickness;
   delta_ = 0.5*thickness_;
   pitch_ = l_/n_;
   h_ = 0.5*l_/(constants::pi * n_);
   offset_ = 0.0;

   logger(DEBUG, "ScrewCC(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
ScrewCC::~ScrewCC()
{
   logger(DEBUG, "~ScrewCC() finished, destroyed the screw.");
}

// A pointer to the copy of the screw
ScrewCC* ScrewCC::copy() const
{
   return new ScrewCC(*this);
}

// Function to set the screw parameters
void ScrewCC::set(Vec3D start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness)
{
   start_ = start;
   l_ = length;
   rMax_ = bladeRadius;
   rMin_ = shaftRadius;
   n_ = numberOfTurns;
   omega_ = omega;
   thickness_ = thickness;
   delta_ = 0.5 * thickness_;
   pitch_ = l_/n_;
   h_ = 0.5 * l_/(constants::pi * n_);
   offset_ = 0.0;

   std::cout << "\n\nScrew parameters set.\n";
   std::cout << "Start : " << start_ << "\n";
   std::cout << "Length : " << l_ << "\n";
   std::cout << "Blade radius : " << rMax_ << "\n";
   std::cout << "Shaft radius : " << rMin_ << "\n";
   std::cout << "Number of turns : " << n_ << "\n";
   std::cout << "Angular velocity : " << omega_ << "\n";
   std::cout << "Thickness : " << thickness_ << "\n";
   std::cout << "Half thickness : " << delta_ << "\n";
   std::cout << "Pitch length : " << pitch_ << "\n";
   std::cout << "Rescaled length : " << h_ << "\n";
   std::cout << "Initial angular offset : " << offset_ << "\n";
}

// Function to change the screw radius
void ScrewCC::setRadius(Mdouble bladeRadius)
{
   rMax_ = bladeRadius;
}

// Function to change the thickness (and the half-thickness as well)
void ScrewCC::setThickness(Mdouble thickness)
{
   thickness_ = thickness;
   delta_ = 0.5 * thickness_;
}

// Function to change the screw rotation velocity
void ScrewCC::setOmega(Mdouble omega)
{
   omega_ = omega;
}

// Function to make the screw rotate by adding \omega*t to the angular offset
void ScrewCC::rotate(Mdouble dt)
{
   offset_ += omega_ * dt;
}

// Function that returns the current angular offset
Mdouble ScrewCC::getOffset()
{
   return offset_;
}

// Function to change the screw angular offset
void ScrewCC::setOffset(Mdouble off)
{
   offset_ = off;
}

// Function to make the screw rotate by incrementing the angular offset
void ScrewCC::incrementOffset(Mdouble off)
{
   offset_ += off;
}

// Input stream
void ScrewCC::read(std::istream& is)
{
   BaseWall::read(is);
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> l_
   >> dummy >> rMax_
   >> dummy >> rMin_
   >> dummy >> n_
   >> dummy >> omega_
   >> dummy >> thickness_
   >> dummy >> offset_;
}

// Old input stream (no idea what this is)
void ScrewCC::oldRead(std::istream& is)
{
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> l_
   >> dummy >> rMax_
   >> dummy >> rMin_
   >> dummy >> n_
   >> dummy >> omega_
   >> dummy >> thickness_
   >> dummy >> offset_;
}

// Writer
void ScrewCC::write(std::ostream& os) const
{
   BaseWall::write(os);
   os << " start " << start_
   << " Length " << l_
   << " Max radius " << rMax_
   << " Min radius " << rMin_
   << " Revolutions " << n_
   << " Omega " << omega_
   << " Thickness " << thickness_
   << " Offset " << offset_;
}

// Returns the string "ScrewCC"
std::string ScrewCC::getName() const
{
   return "ScrewCC";
}

// The contact detection algorithm
// All the quantities in cylindrical coordinates are referred to the screw axis
bool ScrewCC::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
   // computing particle position in cylindrical coordinates
   // IMPORTANT: the angle must be defined in the interval [0, +2*pi[ radians!
   Mdouble r_P = sqrt(pow(p.getPosition().X - start_.X, 2.0) + pow(p.getPosition().Y - start_.Y, 2.0));
   Mdouble theta_P = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
   if (theta_P < 0.0) theta_P += 2.0*constants::pi;
   Mdouble z_P = p.getPosition().Z - start_.Z;

   // if the particle is outside the cylinder containing the screw there is no collision
   if (r_P > rMax_ + p.getWallInteractionRadius()) return false;
   if (z_P > l_ + p.getWallInteractionRadius()) return false;
   if (z_P < - p.getWallInteractionRadius()) return false;

   // trigonometric functions of the helix angle at the particle position
   Mdouble cosEta = r_P/sqrt(pow(r_P, 2.0) + pow(h_, 2.0));

   // half-pitch of the screw
   Mdouble halfPitch_ = 0.5*pitch_;

   // squared rescaled inverse half-thickness lambda = (1 - thickness/2/(p/2))^2 = (1 - delta/(h_*pi))^2
   Mdouble lambda_ = pow(1.0 - delta_/halfPitch_, 2.0);



   // stuff needed for normal computation
   // the radial variable R(z,xi) as a function of axial position and angle (from the screw parametrization)
   Mdouble r_blade = rMin_ + (rMax_ - rMin_)*(fmod(z_P - h_*(theta_P + offset_), halfPitch_)/halfPitch_ - 1.0)/lambda_;

   // derivative of R(z,xi) as a function of z
   Mdouble dRdZ = 2.0*(rMax_ - rMin_)*(fmod(z_P - h_*(theta_P + offset_), halfPitch_)/halfPitch_ - 1.0)/lambda_/halfPitch_;



   // computing radial and axial distance of particle from blade surface
   // --- AXIAL DISTANCE ---
   Mdouble deltaZ = fmod(z_P - h_*(theta_P + offset_ + 2.0*constants::pi*((int)(z_P/pitch_))),pitch_);
   if (deltaZ > 0.5*pitch_) deltaZ -= pitch_;
   if (deltaZ < -0.5*pitch_) deltaZ += pitch_;

   // determining correct sign of the normal
   Mdouble normalSign = -1.0;
   if (deltaZ < 0.0) normalSign = 1.0;

   // --- DISTANCE ALONG NORMAL TO HELICOID ---
   Mdouble deltaN = fabs(deltaZ)*cosEta - delta_;

   // trigonometric functions relative to the particle angle
   Mdouble cosTheta_P = (p.getPosition().X - start_.X)/r_P;
   Mdouble sinTheta_P = (p.getPosition().Y - start_.Y)/r_P;



   // radial vector at the particle position
   Vec3D radialVector;
   radialVector.X = cosTheta_P;
   radialVector.Y = sinTheta_P;
   radialVector.Z = 0.0;

   // normal to the blade at the particle position
   Vec3D normalVector;
   normalVector.X = normalSign*dRdZ*h_*sinTheta_P;
   normalVector.Y = -normalSign*dRdZ*h_*cosTheta_P;
   normalVector.Z = normalSign*r_blade*dRdZ;

   // adds a radial component
   normalVector -= r_blade*radialVector;

   // normalizes and takes the right direction of the vector
   normalVector = normalVector/normalVector.getLength();



   // check which type of collision it is
   if (r_P >= rMax_) // collision with the blade edge
   {
      if (deltaN <= 0.0) // radial collision
      {
         // the distance between the collision point and the particle centre
         distance = r_P - rMax_;

         // if the distance is negative prints an error message
         if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SIDE: OVERLAP > PARTICLE RADIUS\n";

         // the collision has no other component other than the normal to the blade
         normal_return = -radialVector;

         return true;
      }
      else if (deltaN < p.getWallInteractionRadius()) // complex collision
      {
         // radial distance from screw edge
         Mdouble deltaR = r_P - rMax_;

         // distance between the particle centre and the and the edge
         distance = sqrt(pow(deltaN, 2) + pow(deltaR, 2));

         // if the particle-blade_edge distance is higher than the particle radius there is no collision
         if (distance > p.getWallInteractionRadius()) return false;

         // if the distance is negative prints an error message
         if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW EDGE: OVERLAP > PARTICLE RADIUS\n";

         // the normal to the edge through the particle centre
         // IMPORTANT: it needs to be normalized!
         normal_return = deltaN*normalVector - deltaR*radialVector;
         normal_return /= normal_return.getLength();

         return true;
      }
      else // no collision
      {
         return false;
      }
   }
   else // collision with the blade surface
   {
      // distance between particle and radial projection on blade (computed via helicoid(z_P,theta_P))
      Mdouble deltaR = fabs(r_P - r_blade);

      // distance between particle and axial projection on blade (computed via helicoid(r_P,theta_P))
      deltaN = h_*(theta_P + offset_ + constants::pi*(1.0 + normalSign*sqrt((r_P - rMin_)/(rMax_ - rMin_))))*cosEta;

      // if some of the distances are negative prints an error message
      if (deltaR < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SURFACE: r_P < r_0\n";
      if (deltaN < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SURFACE: DeltaN < 0\n";

      // computes the approximated distance from the blade surface
      // the distance between the collision point and the particle centre
      distance = sqrt(pow(deltaR, 2.0) + pow(deltaN, 2.0));

      // checkes for collision
      if (distance > p.getWallInteractionRadius()) return false;

      normal_return = normalVector;

      return true;
   }

   // error message if some case has not been accounted for
   std::cout << "\nERROR: THE SCREW COLLISION FUNCTION DID NOT RETURN PROPERLY\n";

   return false;
}

// Checks for the interaction between a particle p at a time timeStamp
// In case of interaction returns a pointer to the BaseInteraction happened between the Screw and the BaseParticle at time timeStamp
std::vector<BaseInteraction *>
ScrewCC::getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler)
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
