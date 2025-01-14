#include "HopperAndCasing.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

/*
   IMPORTANT: this implementation still needs the back wall at z=0.0 to work.
   The mixed collision between casing and back wall is somehow not working yet.
*/

HopperAndCasing::HopperAndCasing()
{
   start_.setZero(); // the origin of the geometry
   radius_ = 1.0; // the radius of the screw casing
   totalLength_ = 1.0; // the total length
   hopperLength_ = 0.5; // the length of teh hopper in the screw axial direction
   hopperHeight_ = 1.0;
   logger(DEBUG, "HopperAndCasing() constructor finished.");
}

HopperAndCasing::HopperAndCasing(const HopperAndCasing& other)
: BaseWall(other)
{
   start_ = other.start_;
   radius_ = other.radius_;
   totalLength_ = other.totalLength_;
   hopperLength_ = other.hopperLength_;
   hopperHeight_ = other.hopperHeight_;
   logger(DEBUG, "HopperAndCasing(const Shaft&) copy constructor finished.");
}

HopperAndCasing::HopperAndCasing(Vec3D start, Mdouble r, Mdouble lTot, Mdouble l, Mdouble h)
{
   start_ = start;
   radius_ = r;
   totalLength_ = lTot;
   hopperLength_ = l;
   hopperHeight_ = h;
   logger(DEBUG, "HopperAndCasing(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

HopperAndCasing::~HopperAndCasing()
{
   logger(DEBUG, "~HopperAndCasing() finished, destroyed the HopperAndCasing.");
}

HopperAndCasing* HopperAndCasing::copy() const
{
   return new HopperAndCasing(*this);
}

void HopperAndCasing::set(Vec3D start, Mdouble r, Mdouble lTot, Mdouble l, Mdouble h)
{
   start_ = start;
   radius_ = r;
   totalLength_ = lTot;
   hopperLength_ = l;
   hopperHeight_ = h;

   std::cout << std::endl << std::endl << "HopperAndCasing parameters set." << std::endl;
   std::cout << "Start : " << start_ << std::endl;
   std::cout << "Radius : " << radius_ << std::endl;
   std::cout << "Total length : " << totalLength_ << std::endl;
   std::cout << "Hopper length : " << hopperHeight_ << std::endl;
}

// the contact detection algorithm
bool HopperAndCasing::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
   /* divide the geometry in 3 parts:
      1) the pure hopper
      2) the pure casing
      3) the transition part underneath the hopper and before the casing
   */

   // rescale the particle position with respect of the starting point
   Vec3D pPosition = p.getPosition() - start_;

   // PART 1
   if (pPosition.Z < hopperLength_ && pPosition.Y > radius_)
   {
      // collision with back or front wall of the hopper
      if (fabs(pPosition.Z - 0.5*hopperLength_) > 0.5*hopperLength_ - p.getWallInteractionRadius(this))
      {
         // pure collision
         if (fabs(pPosition.X) < radius_ - p.getWallInteractionRadius(this))
         {
            distance = 0.5*hopperLength_ - fabs(pPosition.Z - 0.5*hopperLength_);
            normal_return.X = 0.0;
            normal_return.Y = 0.0;
            normal_return.Z = 1.0;

            // taking the correct sign of the normal
            if (pPosition.Z < 0.5*hopperLength_) normal_return.Z *= -1.0;
            return true;
         }
         else // mixed collision
         {
            distance = sqrt(pow(radius_ - fabs(pPosition.X),2.0) + pow(0.5*hopperLength_ - fabs(pPosition.Z - 0.5*hopperLength_),2.0));
            normal_return.X = (p.getWallInteractionRadius(this) - radius_ + fabs(pPosition.X));
            normal_return.Y = 0.0;
            normal_return.Z = (p.getWallInteractionRadius(this) - 0.5*hopperLength_ + fabs(pPosition.Z - 0.5*hopperLength_));
            normal_return /= normal_return.getLength();

            // taking the correct sign of the normal
            if (pPosition.X < 0.0) normal_return.X *= -1.0;
            if (pPosition.Z < 0.5*hopperLength_) normal_return.Z *= -1.0;
            return true;
         }
      }

      // collison with side walls of the hopper
      if (fabs(pPosition.X) > radius_ - p.getWallInteractionRadius(this))
      {
         // only the pure collision has to be evaluated, since the mixed one is already accounted for
         distance = radius_ - fabs(pPosition.X);
         normal_return.X = 1.0;
         normal_return.Y = 0.0;
         normal_return.Z = 0.0;

         // taking the correct sign of the normal
         if (pPosition.X < 0.0) normal_return.X *= -1.0;
         return true;
      }

      // no collision detected
      return false;
   }

   // PART 2
   if (pPosition.Z > hopperLength_)
   {
      // squared distance from the screw axis
      Mdouble rho = sqrt(pow(pPosition.X,2.0) + pow(pPosition.Y,2.0));

      // collision from the casing inside
      if (rho > radius_ - p.getWallInteractionRadius(this))
      {
         distance = radius_ - rho;
         normal_return.X = pPosition.X/rho;
         normal_return.Y = pPosition.Y/rho;
         normal_return.Z = 0.0;

         return true;
      }

      // no collision detected
      return false;
   }

   // PART 3
   if (pPosition.Z <= hopperLength_ && pPosition.Y <= radius_)
   {
      // bottom part
      if (pPosition.Y < 0.0)
      {
         // squared distance from the screw axis
         Mdouble rho = sqrt(pow(pPosition.X,2.0) + pow(pPosition.Y,2.0));

         // collision with the bottom of the casing
         if (rho > radius_ - p.getWallInteractionRadius(this))
         {
            // pure collision with the half-cilinder on the bottom
            if (pPosition.Z > p.getWallInteractionRadius(this))
            {
               distance = radius_ - rho;
               normal_return.X = pPosition.X/rho;
               normal_return.Y = pPosition.Y/rho;
               normal_return.Z = 0.0;

               return true;
            }
            else // mixed collision
            {
               distance = sqrt(pow(pPosition.Z,2.0) + pow(radius_ - rho,2.0));
               normal_return.X = pPosition.X/rho*(p.getWallInteractionRadius(this) - radius_ + rho);
               normal_return.Y = pPosition.Y/rho*(p.getWallInteractionRadius(this) - radius_ + rho);
               normal_return.Z = -(p.getWallInteractionRadius(this) - pPosition.Z);
               normal_return /= normal_return.getLength();

               return true;
            }
         }

         // collision with the back wall of the hopper
         if (pPosition.Z < p.getWallInteractionRadius(this))
         {
            // only the pure collision has to be evaluated, since the mixed one is already accounted for
            distance = pPosition.Z;
            normal_return.X = 0.0;
            normal_return.Y = 0.0;
            normal_return.Z = -1.0;

            return true;
         }

         // no collision detected
         return false;
      }
      else // top part
      {
         // collision with the back wall
         if (pPosition.Z < p.getWallInteractionRadius(this))
         {
            // pure collision
            if (fabs(pPosition.X) < radius_ - p.getWallInteractionRadius(this))
            {
               distance = pPosition.Z;
               normal_return.X = 0.0;
               normal_return.Y = 0.0;
               normal_return.Z = -1.0;

               return true;
            }
            else // mixed collision
            {
               distance = sqrt(pow(radius_ - pPosition.X,2.0) + pow(pPosition.Z,2.0));
               normal_return.X = (p.getWallInteractionRadius(this) - radius_ + fabs(pPosition.X));
               normal_return.Y = 0.0;
               normal_return.Z = -(p.getWallInteractionRadius(this) - pPosition.Z);
               normal_return /= normal_return.getLength();

               // taking the correct sign of the normal
               if (pPosition.X < 0.0) normal_return.X *= -1.0;
               return true;
            }
         }

         // collision with the front wall
         if (pPosition.Z > hopperLength_ - p.getWallInteractionRadius(this))
         {
            // distance from the screw axis
            Mdouble rho = sqrt(pow(pPosition.X,2.0) + pow(pPosition.Y,2.0));

            // collision with top of cylindrical casing
            if (rho < radius_ && rho > radius_ - p.getWallInteractionRadius(this))
            {
               distance = sqrt(pow(hopperLength_ - pPosition.Z,2.0) + pow(radius_ - rho,2.0));
               normal_return.X = pPosition.X/rho*(p.getWallInteractionRadius(this) - radius_ + rho);
               normal_return.Y = pPosition.Y/rho*(p.getWallInteractionRadius(this) - radius_ + rho);
               normal_return.Z = (p.getWallInteractionRadius(this) - hopperLength_ + pPosition.Z);
               normal_return /= normal_return.getLength();

               return true;
            }

            // collision with lateral walls
            if (fabs(pPosition.X) > radius_ - p.getWallInteractionRadius(this))
            {
               distance = sqrt(pow(radius_ - fabs(pPosition.X),2.0) + pow(hopperLength_ - pPosition.Z,2.0));
               normal_return.X = (p.getWallInteractionRadius(this) - radius_ + fabs(pPosition.X));
               normal_return.Y = 0.0;
               normal_return.Z = (p.getWallInteractionRadius(this) - hopperLength_ + pPosition.Z);
               normal_return /= normal_return.getLength();

               // taking the correct sign of the normal
               if (pPosition.X < 0.0) normal_return.X *= -1.0;
               return true;
            }

            // pure collision
            if (rho > radius_)
            {
               distance = hopperLength_ - pPosition.Z;
               normal_return.X = 0.0;
               normal_return.Y = 0.0;
               normal_return.Z = 1.0;

               return true;
            }
         }

         // pure collision with side walls
         if (fabs(pPosition.X) > radius_ - p.getWallInteractionRadius(this))
         {
            // only the pure collision has to be evaluated, since the mixed ones are already accounted for
            distance = radius_ - fabs(pPosition.X);
            normal_return.X = 1.0;
            normal_return.Y = 0.0;
            normal_return.Z = 0.0;

            // taking the correct sign of the normal
            if (pPosition.X < 0.0) normal_return.X *= -1.0;
            return true;
         }

         // no collision detected
         return false;
      }
   }

   // // the particle was no t in the specified domain (likely to be an error)
   // std::cout << std::endl << "FATAL: COLLISION DIDN'T RETURN PROPERLY" << std::endl;
   // std::cout << pPosition.X/radius_ << " " << pPosition.Y/hopperHeight_ << " " << pPosition.Z/hopperLength_ << "   " << distance/p.getWallInteractionRadius(this) << "   " << distance << "   " << normal_return << std::endl;
   // exit(0);
   return false;
}

// std::cout << pPosition << "   " << distance/p.getWallInteractionRadius(this) << "   " << distance << "   " << normal_return << std::endl;

void HopperAndCasing::read(std::istream& is)
{
   BaseWall::read(is);
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> radius_
   >> dummy >> totalLength_
   >> dummy >> hopperLength_
   >> dummy >> hopperHeight_;
}

void HopperAndCasing::write(std::ostream& os) const
{
   BaseWall::write(os);
   os << " Start " << start_
   << " Radius " << radius_
   << " TotalLength " << totalLength_
   << " HopperLength " << hopperLength_
   << " HopperHeight " << hopperHeight_;
}

std::string HopperAndCasing::getName() const
{
   return "HopperAndCasing";
}

std::vector<BaseInteraction *>
HopperAndCasing::getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler)
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
