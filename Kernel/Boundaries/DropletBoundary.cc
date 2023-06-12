#include <Particles/LiquidFilmParticle.h>
#include <MercuryBase.h>
#include "DropletBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"

void DropletBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
    MercuryBase* dpm = dynamic_cast<MercuryBase*>(getHandler()->getDPMBase());
    logger.assert_always(dpm,"You can only run DropletBoundary with Mercury2D or Mercury3D, not DPMBase");
    // generate new droplets
    generateDroplets_(*this);
    // integrate Newtons equation of motion
    Vec3D g = dpm->getGravity();
    double dt = dpm->getTimeStep();
    for (auto& d : droplets_) {
        Mdouble m = dropletSpecies_->getMassFromRadius(d.radius);
        d.velocity += dt * (d.force + m * g) / m;
        d.force.setZero();
        d.position += dt*d.velocity;
    }
    // check for interaction with particles; this is costly, so we only do iit every 50 time steps
    if (dpm->getNumberOfTimeSteps()%checkCount==0) {
        LiquidFilmParticle p;
        p.setSpecies(dropletSpecies_);
        for (auto &d : droplets_) {
            p.setPosition(d.position);
            p.setVelocity(d.velocity);
            p.setRadius(d.radius);
            auto q = dpm->hGridFindParticleContacts(&p);
            if (q.size() > 0) {
                double liquidVolume = std::pow(2.0 * d.radius, 3) * constants::pi / 6.0;
                for (auto particle : q) {
                    // should that be a static cast?
                    auto lfp = static_cast<LiquidFilmParticle *>(particle);
                    lfp->addLiquidVolume(liquidVolume/q.size());
                }
                absorbedVolume += liquidVolume;
                dropletVolume -= liquidVolume;
                d.radius = 0;
            }
            // check for interactions with walls
            for (BaseWall* w : dpm->wallHandler) {
                //Checks if the particle is interacting with the current wall
                BaseInteraction* i = w->getInteractionWith(&p, dpm->getNumberOfTimeSteps() + 1, &dpm->interactionHandler);
                // if there is a wall interacting with the droplet
                if (i != nullptr) {
                    if (removeDropletsAtWalls_) {
                        //set droplet radius to zero
                        double liquidVolume = std::pow(2.0 * d.radius, 3) * constants::pi / 6.0;
                        lostVolume += liquidVolume;
                        dropletVolume -= liquidVolume;
                        d.radius = 0;
                    } else {
                        // For some reason the species have to be set here, otherwise the force is always 0.
                        i->setSpecies(dropletSpecies_);
                        i->computeForce();
                        d.force += i->getForce();
                    }
                }
            }
        }
        // erase all droplets of radius 0
        droplets_.erase(std::remove_if(droplets_.begin(), droplets_.end(),
                                       [](const Droplet &d) { return d.radius == 0; }), droplets_.end());
    }
}

/*!
 * \details Reads a number of boundary properties from the given std::istream.
 * \param[in,out] is   the istream
 */
void DropletBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    size_t n;
    Vec3D position, velocity;
    double radius;
    is >> dummy >> checkCount;
    is >> dummy >> n;
    droplets_.reserve(n);
    droplets_.resize(0);
    for (int i = 0; i < n; ++i) {
        is >> position >> velocity >> radius;
        droplets_.emplace_back(Droplet{position,velocity,radius});
    }
}

/*!
 * \details Writes the boundary properties to an std::ostream. 
 * \param[out] os   the ostream the properties are to be written to.
 */
void DropletBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " checkCount " << checkCount;
    os << " n " << droplets_.size();
    for (auto& d : droplets_) {
        os << " " << d.position;
        os << " " << d.velocity;
        os << " " << d.radius;
    }
}

void DropletBoundary::writeVTK(std::fstream& file) {
    file << "<Piece NumberOfPoints=\"" << droplets_.size() << "\" NumberOfCells=\"" << 0 << "\">\n";
    file << "<Points>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& d : droplets_) file << '\t' << (float)d.position.X << ' ' << (float)d.position.Y << ' ' << (float)d.position.Z << '\n';
    file << "  </DataArray>\n";
    file << "</Points>\n";
    file << "<PointData  Vectors=\"vector\">\n";
    file << "  <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& d : droplets_) file << '\t' << (float)d.velocity.X << ' ' << (float)d.velocity.Y << ' ' << (float)d.velocity.Z << '\n';
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Float32\" Name=\"Radius\" format=\"ascii\">\n";
    for (const auto& d : droplets_) file << '\t' << (float)d.radius << '\n';
    file << "  </DataArray>\n";
    file << "</PointData>\n";
}
