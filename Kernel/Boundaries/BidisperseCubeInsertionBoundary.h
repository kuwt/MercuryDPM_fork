#ifndef BOUNDARIES_BIDISPERSECUBEINSERTIONBOUNDARY_H
#define BOUNDARIES_BIDISPERSECUBEINSERTIONBOUNDARY_H

#include "CubeInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"
#include "Math/Vector.h"

/*!
 * \class BidisperseCubeInsertionBoundary
 * \brief Like a CubeInsertionBoundary but the particles generated are one
 * of two types.
 */

class BidisperseCubeInsertionBoundary : public CubeInsertionBoundary
{
public:
    /*!
     * \brief Constructor; sets everything to 0.
     */
    BidisperseCubeInsertionBoundary();
    
    /*!
     * \brief Copy constructor with deep copy.
     */
    BidisperseCubeInsertionBoundary(const BidisperseCubeInsertionBoundary& other);
    
    /*!
     * \brief Destructor: default destructor.
     */
    ~BidisperseCubeInsertionBoundary() override;
    
    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    BidisperseCubeInsertionBoundary* copy() const override;
    
    /*!
     * \brief Sets the properties of this bidisperse cuboidal insertion boundary
     */
    void set(BaseParticle* particleToCopyA, BaseParticle* particleToCopyB, double probA, int maxFailed, Vec3D posMin,
             Vec3D posMax, Vec3D velMin, Vec3D velMax);
    
    /*!
     * \brief Get the particles that need to be copied 
     */
    BaseParticle* getParticleToCopyA() const;
    
    BaseParticle* getParticleToCopyB() const;
    
    /*!
     * \brief Generates a particle with random position, radius and velocity 
     */
    BaseParticle* generateParticle(RNG& random) override;
    
    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

private:
    BaseParticle* particleToCopyA_;
    BaseParticle* particleToCopyB_;
    double probA_;
};

#endif
