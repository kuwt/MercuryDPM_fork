#ifndef VCHUTE_H
#define VCHUTE_H

#include "BaseWall.h"
#include "Math/Vector.h"

/* A finite-length V-chute. */

class VChute : public BaseWall
{
public:
    VChute();
    
    VChute(const VChute& other);
    
    VChute(Mdouble length, Mdouble width, Mdouble alpha);
    
    ~VChute() override;
    
    void set(Mdouble length, Mdouble width, Mdouble alpha);
    
    VChute* copy() const override;
    
    // wheeee virtual mc virty-virtual
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    BaseInteraction*
    getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    void read(std::istream& is) override;
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Mdouble l_; // length
    Mdouble w_; // width
    Mdouble alpha_; // alpha, in radians (V-ness)
};

#endif
