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

        // Used for when to repel droplet at walls.
        Vec3D force;
    };

    std::vector<Droplet> droplets_;

    // volume of liquid stored in droplets
    double dropletVolume = 0;
    // volume of liquid absorbed by particles
    double absorbedVolume = 0;
    // volume of liquid lost at walls
    double lostVolume = 0;

private:

    std::function<void(DropletBoundary&)> generateDroplets_
        = [] (DropletBoundary&) {};

public:

    DropletBoundary() {}
    
    DropletBoundary(const DropletBoundary& other) {
        droplets_ = other.droplets_;
        generateDroplets_ = other.generateDroplets_;
        removeDropletsAtWalls_ = other.removeDropletsAtWalls_;
        dropletSpecies_ = other.dropletSpecies_;
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

    // this is the number of timesteps between checks whether the droplet is in contact with a particle or wall.
    unsigned checkCount = 3;

    void setRemoveDropletsAtWalls(bool removeDroplets) {
        removeDropletsAtWalls_ = removeDroplets;
    }

    void setDropletSpecies(const ParticleSpecies* species) {
        dropletSpecies_ = species;
    }

    void actionsBeforeTimeLoop() override
    {
        // When no droplet species were set, use the last species from the handler.
        // This is fine for when the droplets are removed at the walls, if not, a warning is shown.
        if (!dropletSpecies_)
        {
            dropletSpecies_ = getHandler()->getDPMBase()->speciesHandler.getLastObject();
            if (!removeDropletsAtWalls_)
                logger(WARN, "DropletBoundary: Droplets should repel from wall, but no droplet species was set. Using last species from the species handler instead.");
        }

        // When droplets should repel from the wall, it requires to be checked every time step.
        if (!removeDropletsAtWalls_)
            checkCount = 1;
    }

private:
    bool removeDropletsAtWalls_ = true;
    const ParticleSpecies* dropletSpecies_ = nullptr;
};

#endif
