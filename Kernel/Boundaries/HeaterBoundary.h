#ifndef HEATERBOUNDARY_H
#define HEATERBOUNDARY_H

#include "DPMBase.h"
#include "BaseBoundary.h"
#include "Math/Vector.h"
#include "Math/RNG.h"

class DPMBase;

class ParticleHandler;

class BaseParticle;

class RNG;

/*!
 * \class HeaterBoundary
 * \brief Supplies a 'constant heat flux' to a cuboidal region (specified by two
 * corner points) by adding a random velocity at each time step to each particle
 * therein, increasing the granular temperature (velocity variance). 
 */

class HeaterBoundary : public BaseBoundary
{
public:
    HeaterBoundary();
    
    HeaterBoundary(const HeaterBoundary& other);
    
    ~HeaterBoundary() override;
    
    HeaterBoundary* copy() const override;
    
    void set(Vec3D posMin, Vec3D posMax, Vec3D specificHeatStrength);
    
    void set2D(Vec3D posMin, Vec3D posMax, Mdouble specificHeatStrength);
    
    void set3D(Vec3D posMin, Vec3D posMax, Mdouble specificHeatStrength);
    
    void setStrength(Vec3D specificHeatStrength);
    
    void setStrength2D(Mdouble specificHeatStrength);
    
    void setStrength3D(Mdouble specificHeatStrength);
    
    /*!
     * \brief Returns the volume of the HeaterBoundary cuboid.
     */
    inline Mdouble getVolume() const;
    
    /*!
     * \brief Returns a negative value if and only if the position is inside the
     * boundary (and therefore the particle to be heated).
     */
    Mdouble getDistance(const Vec3D& position) const;
    
    /*!
     * \brief Runs at the end of each time step.
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!
     * \brief Checks if a given particle is inside the HeaterBoundary. If
     * so, heats it.
     */
    bool checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH);
    
    /*!
     * \brief Reads some boundary properties from an std::istream.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes the boundary properties to an std::ostream.
     */
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Vec3D posMin_, posMax_;
    
    /*!
     * \brief The specific heat strengths have units of energy (mass)^{-1} (time)^{-1) = m^2 s^{-3}
     * For isotropic heating, take a vector proportional to (s/3, s/3, s/3).
     */
    Vec3D specificHeatStrength_;
    
};

#endif
