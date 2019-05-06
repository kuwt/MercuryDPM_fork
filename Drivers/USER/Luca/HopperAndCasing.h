#ifndef HopperAndCasing_H
#define HopperAndCasing_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class HopperAndCasing : public BaseWall
{
public:

   HopperAndCasing();

   HopperAndCasing(const HopperAndCasing& other);

   HopperAndCasing(Vec3D start, Mdouble r, Mdouble lTot, Mdouble l, Mdouble h);

   ~HopperAndCasing();

   HopperAndCasing* copy() const final;

   void set(Vec3D start, Mdouble r, Mdouble lTot, Mdouble l, Mdouble h);

   bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

   void read(std::istream& is) override;

   void write(std::ostream& os) const override;

   std::string getName() const final;

   std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler);

private:
   Vec3D start_; // the origin of the geometry
   Mdouble radius_; // the radius of the screw casing
   Mdouble totalLength_; // the total length
   Mdouble hopperLength_; // the length of teh hopper in the screw axial direction
   Mdouble hopperHeight_; // the height of the hopper over the middle plane containing the screw axis
};

#endif
