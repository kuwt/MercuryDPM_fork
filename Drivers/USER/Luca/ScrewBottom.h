#ifndef ScrewBottom_H
#define ScrewBottom_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class ScrewBottom : public BaseWall
{
public:

   ScrewBottom();

   ScrewBottom(const ScrewBottom& other);

   ScrewBottom(Vec3D start, Mdouble r, Mdouble l);

   ~ScrewBottom();

   ScrewBottom* copy() const final;

   void set(Vec3D start, Mdouble r, Mdouble l);

   bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

   void read(std::istream& is) override;

   void write(std::ostream& os) const override;

   std::string getName() const final;

   std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler);

private:
   Vec3D start_; // the origin of the geometry
   Mdouble radius_; // the radius of the screw casing
   Mdouble length_; // the total length
};

#endif
