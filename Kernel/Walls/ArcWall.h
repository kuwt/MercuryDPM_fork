#ifndef ARCWALL_H
#define ARCWALL_H

#include "Walls/BaseWall.h"
#include <math.h>

class WallHandler;
class BaseParticle;
class BaseWall;

/*!
 * \brief A wall that is the inside of an arc of a cylinder.
 * \details The ArcWall is specified by the cylinder's axis (in turn by a
 * position vector and a direction vector), its radius, a 'centreline direction'
 * (a unit vector normal to the direction vector) for the centreline of the arc,
 * and a semiangle. 
 * For a wall that is the outside of an arc, use AxisymmetricIntersectionOfWalls
 * or Combtooth instead.
 */
class ArcWall : public BaseWall
{
    public:
        /*!
         * \brief Default constructor
         */
        ArcWall();

        /*!
         * \brief Copy constructor
         */
        ArcWall(const ArcWall& aw);

        /*!
         * \brief Default destructor
         */
        virtual ~ArcWall();

        /*!
         * \brief Set 
         */
        void set(Vec3D axis, Vec3D pos, Mdouble radius, Vec3D centreline, Mdouble semiangle);

        ArcWall* copy() const override;

        bool getDistanceAndNormal(const BaseParticle& p,
                Mdouble& distance, Vec3D& normal_return) const override;

        std::vector<BaseInteraction*> getInteractionWith(BaseParticle* p,
                unsigned timeStamp, InteractionHandler* interactionHandler) override;

        void read(std::istream& is) override;
        void write(std::ostream& os) const override;
        std::string getName() const override;

    private:
        Vec3D axis_;
        Vec3D pos_;
        Mdouble radius_;
        Vec3D centreline_;
        Mdouble semiangle_;

};
#endif
