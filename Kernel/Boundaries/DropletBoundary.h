#ifndef DROPLETBOUNDARY_H
#define DROPLETBOUNDARY_H

#include "DPMBase.h"
#include "BaseBoundary.h"
#include "Math/Vector.h"
#include "Math/RNG.h"

class DPMBase;

class ParticleHandler;

class BaseParticle;

class RNG;

/*!
 * \class DropletBoundary
 * \brief Supplies a 'constant heat flux' to a cuboidal region (specified by two
 * corner points) by adding a random velocity at each time step to each particle
 * therein, increasing the granular temperature (velocity variance).
 * \details Note, you need to create a species for the droplets that has liquidVolumeMax-0, or the contact happens at a non-zero distance
 */

class DropletBoundary : public BaseBoundary
{
public:

    struct Droplet {
        Vec3D position;
        Vec3D velocity;
        double radius;
    };

    std::vector<Droplet> droplets_;

private:

    std::function<void(DropletBoundary&)> generateDroplets_
        = [] (DropletBoundary&) {};

public:

    DropletBoundary() {}
    
    DropletBoundary(const DropletBoundary& other) {
        droplets_ = other.droplets_;
        generateDroplets_ = other.generateDroplets_;
    }
    
    ~DropletBoundary() override {
        logger(VERBOSE, "A DropletBoundary has been destroyed.");
    }
    
    DropletBoundary* copy() const override {
        return new DropletBoundary(*this);
    }

    std::string getName() const override {
        return "DropletBoundary";
    }

    /*!
     * \brief Runs at the end of each time step.
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;

    /*!
     * \brief Reads some boundary properties from an std::istream.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes the boundary properties to an std::ostream.
     */
    void write(std::ostream& os) const override;

    void setGenerateDroplets(std::function<void(DropletBoundary&)> generateDroplets) {
        generateDroplets_=generateDroplets;
    }

    void writeVTK(std::fstream& file) override;

};

#endif
