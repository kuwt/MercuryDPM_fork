#ifndef COMBTOOTH_H
#define COMBTOOTH_H

#include "BaseWall.h"
#include "Math/Vector.h"
#include <cmath>

class Combtooth : public BaseWall
{
public:
    /*!
     * \brief Default constructor
     */
    Combtooth();
    
    /*!
     * \brief Copy constructor
     */
    Combtooth(const Combtooth& ct);
    
    /*!
     * \brief Default destructor
     */
    ~Combtooth() override;
    
    /*!
     * \brief Set
     */
    void set(Vec3D axis, Vec3D position, Mdouble radius);
    
    /*!
     * \brief Copy
     */
    Combtooth* copy() const override;
    
    bool getDistanceAndNormal(const BaseParticle& p,
                              Mdouble& distance, Vec3D& normal_return) const override;
    
    std::vector<BaseInteraction*> getInteractionWith(BaseParticle* p,
                                                     unsigned timeStamp,
                                                     InteractionHandler* interactionHandler) override;
    
    void read(std::istream& is) override;
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Vec3D axis_; // unit vector pointing in direction of axis
    Vec3D position_; // position vector of a point that the axis goes through
    Mdouble radius_;
    
};

#endif
