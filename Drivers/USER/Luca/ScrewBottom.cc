#include "ScrewBottom.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

/*
   IMPORTANT: this implementation still needs the back wall at z=0.0 to work.
   The mixed collision between casing and back wall is somehow not working yet.
*/

ScrewBottom::ScrewBottom()
{
   start_.setZero(); // the origin of the geometry
   radius_ = 1.0; // the radius of the screw casing
   length_ = 1.0; // the total length
   logger(DEBUG, "HopperAndCasing() constructor finished.");
}

ScrewBottom::ScrewBottom(const ScrewBottom& other)
: BaseWall(other)
{
   start_ = other.start_;
   radius_ = other.radius_;
   length_ = other.length_;
   logger(DEBUG, "ScrewBottom(const ScrewBottom&) copy constructor finished.");
}

ScrewBottom::ScrewBottom(Vec3D start, Mdouble r, Mdouble l)
{
   start_ = start;
   radius_ = r;
   length_ = l;
   logger(DEBUG, "ScrewBottom(Vec3D, Mdouble, Mdouble) constructor finished.");
}

ScrewBottom::~ScrewBottom()
{
   logger(DEBUG, "~ScrewBottom() finished, destroyed the ScrewBottom.");
}

ScrewBottom* ScrewBottom::copy() const
{
   return new ScrewBottom(*this);
}

void ScrewBottom::set(Vec3D start, Mdouble r, Mdouble l)
{
   start_ = start;
   radius_ = r;
   length_ = l;

   std::cout << std::endl << std::endl << "ScrewBottom parameters set." << std::endl;
   std::cout << "Start : " << start_ << std::endl;
   std::cout << "Radius : " << radius_ << std::endl;
   std::cout << "Length : " << length_ << std::endl;
}

// the contact detection algorithm
bool ScrewBottom::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
   // the collision happens when the particle touches the bottom semi-cylinder

   // rescale the particle position with respect of the starting point
   Vec3D pPosition = p.getPosition() - start_;

   if (pPosition.Y < 0.0)
   {
      // squared distance from the screw axis
      Mdouble rho = sqrt(pow(pPosition.X,2.0) + pow(pPosition.Y,2.0));

      // collision from the casing inside
      if (rho > radius_ - p.getWallInteractionRadius())
      {
         distance = radius_ - rho;
         normal_return.X = pPosition.X/rho;
         normal_return.Y = pPosition.Y/rho;
         normal_return.Z = 0.0;

         return true;
      }

      return false;
   }
   else
   {
      if (fabs(pPosition.X) > radius_ - p.getWallInteractionRadius())
      {
         distance = radius_ - fabs(pPosition.X);
         normal_return.X = pPosition.X;
         normal_return.Y = pPosition.Y;
         normal_return.Z = 0.0;
         normal_return /= normal_return.getLength();

         return true;
      }

      return false;
   }

   std::cout << "COLLISION ERROR!" << std::endl;
}

void ScrewBottom::read(std::istream& is)
{
   BaseWall::read(is);
   std::string dummy;
   is >> dummy >> start_
   >> dummy >> radius_
   >> dummy >> length_;
}

void ScrewBottom::write(std::ostream& os) const
{
   BaseWall::write(os);
   os << " Start " << start_
   << " Radius " << radius_
   << " Length " << length_;
}

std::string ScrewBottom::getName() const
{
   return "ScrewBottom";
}

std::vector<BaseInteraction *>
ScrewBottom::getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler)
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
