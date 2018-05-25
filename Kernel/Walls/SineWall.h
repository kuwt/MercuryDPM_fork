/* A half-infinite surface which has a sinusoidal direction in the x direction,
 * and no variation in the y direction.
 *
 * The equation of the surface is
 *   z = sw_amp * sin( sw_wavn * x + sw_phshift )  for x < l
 * Because this is nonlinear, we must solve for the interaction location using
 * Newton's method. 
 * */

#ifndef SINEWALL_H
#define SINEWALL_H

#include "BaseWall.h"
#include "Math/Vector.h"

class SineWall : public BaseWall
{
public:
    /*!
     * \brief Default constructor, sets a chute with default parameters.
     */
    SineWall();
    
    /*!
     * \brief Copy constructor
     */
    SineWall(const SineWall& other);
    
    /*!
     * \brief Constructor in which all parameters are set.
     */
    SineWall(Mdouble length, Mdouble sw_wavn, Mdouble sw_phshift, Mdouble sw_amp);
    
    ~SineWall() override;
    
    void set(Mdouble length, Mdouble sw_wavn, Mdouble sw_phshift, Mdouble sw_amp);
    
    SineWall* copy() const override;
    
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const override;
    
    std::vector<BaseInteraction*>
    getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    void read(std::istream& is) override;
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Mdouble l_;
    Mdouble sw_amp_; // amplitude of sinusoidal variations
    Mdouble sw_wavn_; // (angular) wavenumber of sinusoidal variations
    Mdouble sw_phshift_; // phase shift of sinusoidal variations
};

#endif
