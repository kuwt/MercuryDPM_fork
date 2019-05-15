//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Logger.h>
#include "ParticleSpecies.h"
#include "DPMBase.h"
#include "Particles/BaseParticle.h"
#include "Particles/SuperQuadricParticle.h"

class BaseParticle;

class BaseInteractable;

ParticleSpecies::ParticleSpecies()
{
    density_ = 1.0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ParticleSpecies::ParticleSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p the species that is copied
 */
ParticleSpecies::ParticleSpecies(const ParticleSpecies& p)
{
    density_ = p.density_;
    temperatureDependentDensity_ = p.temperatureDependentDensity_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ParticleSpecies::ParticleSpecies(const ParticleSpecies &p) finished"<<std::endl;
#endif
}

ParticleSpecies::~ParticleSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"ParticleSpecies::~ParticleSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void ParticleSpecies::write(std::ostream& os) const
{
    //note we inherit from BaseObject, not BaseParticle
    BaseObject::write(os);
    os << " density " << density_;
    if (getConstantRestitution()) {
        os << " constantRestitution " << getConstantRestitution();
    }
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void ParticleSpecies::read(std::istream& is)
{
    BaseSpecies::read(is);
    std::string dummy;
    bool constantRestitution;
    is >> dummy >> density_;
    if (helpers::readOptionalVariable(is, "constantRestitution", constantRestitution)) {
        setConstantRestitution(constantRestitution);
    }
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string ParticleSpecies::getBaseName() const
{
    return "Particle";
}

/*!
 * \param[in] density the particle density
 */
void ParticleSpecies::setDensity(Mdouble density)
{
    logger.assert_always(density >= 0, "[ParticleSpecies::setDensity(%)] value cannot be negative", density);
    density_ = density;
    if (getHandler()) getHandler()->getDPMBase()->particleHandler.computeAllMasses(getIndex());
}

/*!
 * \return the particle density
 */
Mdouble ParticleSpecies::getDensity() const
{
    return density_;
}

Mdouble ParticleSpecies::getMassFromRadius(const Mdouble radius) const
{
    return getDensity() * getVolumeFromRadius(radius);
}

///\todo this should depend on the particle shape; thus, it should be a static function of BaseParticle
Mdouble ParticleSpecies::getVolumeFromRadius(const Mdouble radius) const
{
    if (getHandler() == nullptr)
    {
        logger(ERROR,
               "[Species::VolumeFromRadius()] No handler has been set, therefore, I can't figure out the dimensions.");
        return 0;
    }
    
    unsigned int particleDimensions = getHandler()->getDPMBase()->getParticleDimensions();
    if (particleDimensions == 3)
    {
        return 4.0 / 3.0 * constants::pi * radius * radius * radius;
    }
    else if (particleDimensions == 2)
    {
        return constants::pi * radius * radius;
    }
    else if (particleDimensions == 1)
    {
        return 2.0 * radius;
    }
    else
    {
        logger(ERROR, "[Species::VolumeFromRadius()] the dimension of the particle is wrongly set to %",
               particleDimensions);
        return 0.0;
    }
}

///Compute BaseParticle mass function, which required a reference to the Species vector. It computes the Particles mass, Inertia and the inverses.
/// this function is called, if BaseParticleHandler::addObject, SpeciesHandler::addObject, ParticleSpecies::setDensity, BaseParticle::setRadius or DPMBase::setParticleDimensions is called
void ParticleSpecies::computeMass(BaseParticle* p) const
{
    if (!p->isFixed())
    {
        switch (p->getParticleDimensions())
        {
            case 3:
            {
                p->invMass_ = 1.0 / (4.0 / 3.0 * constants::pi * p->getRadius() * p->getRadius() * p->getRadius() *
                                     getDensity());
                p->invInertia_ =
                        MatrixSymmetric3D(1, 0, 0, 1, 0, 1) / (.4 * p->getMass() * mathsFunc::square(p->getRadius()));
                //
                if (p->getName() == "SuperQuadricParticle")
                {
                    SuperQuadricParticle* SE = dynamic_cast<SuperQuadricParticle*>(p);
                    Vec3D axes = SE->getAxes();
                    Mdouble eps1 = SE->getExponentEps1();
                    Mdouble eps2 = SE->getExponentEps2();
                    Mdouble volume = SE->getVolume();
                    
                    Mdouble help1 = mathsFunc::beta(1.5 * eps2, 0.5 * eps2);
                    Mdouble help2 = mathsFunc::beta(0.5 * eps1, 2.0 * eps1 + 1.0);
                    Mdouble help3 = mathsFunc::beta(0.5 * eps2, 0.5 * eps2 + 1);
                    Mdouble help4 = mathsFunc::beta(1.5 * eps1, eps1 + 1.0);
                    
                    p->invMass_ = 1.0 / (volume * getDensity());
                    p->invInertia_.XX = 1.0 / (getDensity() * (0.5 * axes.X * axes.Y * axes.Z * eps1 * eps2) *
                                               (axes.Y * axes.Y * help1 * help2
                                                + 4.0 * axes.Z * axes.Z * help3 * help4));
                    p->invInertia_.XY = 0.0;
                    p->invInertia_.XZ = 0.0;
                    p->invInertia_.YY = 1.0 / (getDensity() * (0.5 * axes.X * axes.Y * axes.Z * eps1 * eps2) *
                                               (axes.X * axes.X * help1 * help2
                                                + 4.0 * axes.Z * axes.Z * help3 * help4));
                    p->invInertia_.YZ = 0.0;
                    p->invInertia_.ZZ = 1.0 / (getDensity() * (0.5 * axes.X * axes.Y * axes.Z * eps1 * eps2) *
                                               (axes.X * axes.X + axes.Y * axes.Y) * help1 * help2);
                }
                break;
            }
            
            case 2:
            {
                p->invMass_ = 1.0 / (constants::pi * p->getRadius() * p->getRadius() * getDensity());
                p->invInertia_ =
                        MatrixSymmetric3D(1, 0, 0, 1, 0, 1) / (.5 * p->getMass() * mathsFunc::square(p->getRadius()));
                break;
            }
            
            case 1:
            {
                p->invMass_ = 1.0 / (2.0 * p->getRadius() * getDensity());
                p->invInertia_ = MatrixSymmetric3D(1, 0, 0, 1, 0, 1) / std::numeric_limits<Mdouble>::quiet_NaN();
                break;
            }
            
            default:
            {
                logger(ERROR, "ParticleSpecies::computeMass()] the dimension of the particle is not set");
            }
        }
    }
}

const std::function<double(double)>& ParticleSpecies::getTemperatureDependentDensity() const
{
    return temperatureDependentDensity_;
}

void ParticleSpecies::setTemperatureDependentDensity(
        const std::function<double(double)>& temperatureDependentDensity)
{
    temperatureDependentDensity_ = temperatureDependentDensity;
//    density_ = temperatureDependentDensity_(0);
//    logger(INFO,"Setting initial density to %",temperature_);
}

/*!
 * \return A pointer to the to the lightest BaseParticle (by mass) in this ParticleHandler.
 */
Mdouble ParticleSpecies::getLargestInverseParticleMassLocal() const
{
    Mdouble maxInvMass = 0;
    logger.assert(getHandler() != nullptr && getHandler()->getDPMBase() != nullptr,"speciesHandler must be set");
    for (BaseParticle* const p : getHandler()->getDPMBase()->particleHandler)
    {
        if (p->getSpecies()==this && !(p->isFixed() || p->isMPIParticle() || p->isPeriodicGhostParticle()) && p->getInvMass() > maxInvMass)
        {
            maxInvMass = p->getInvMass();
        }
    }
    return maxInvMass;
}

Mdouble ParticleSpecies::getSmallestParticleMass() const
{
#ifdef MERCURY_USE_MPI
    Mdouble maxInvMass = 0;
    Mdouble invMassLocal = getLargestInverseParticleMassLocal();
    //Obtain the global value
    MPIContainer::Instance().allReduce(invMassLocal, maxInvMass, MPI_MAX);
    //return value
    return 1.0 / maxInvMass;
#else
    return 1.0 / getLargestInverseParticleMassLocal();
#endif

}

/**
 * \details Sets #maxInteractionDistance_
 * @param interactionDistance the interaction distance that has been changed
 */
void ParticleSpecies::setMaxInteractionDistance(Mdouble interactionDistance) {
    // if maxInteractionDistance_ has increased it's simple
    if (interactionDistance>=maxInteractionDistance_) {
        maxInteractionDistance_ = interactionDistance;
    } else /*else we need to recompute*/ {
        maxInteractionDistance_ = 0;
        int j = getIndex();
        for (int i=0; i<getHandler()->getSize(); ++i) {
            maxInteractionDistance_ = std::max(maxInteractionDistance_,getHandler()->getMixedObject(i,j)->getInteractionDistance());
        }
    }
}

const BaseSpecies* ParticleSpecies::getMixedSpecies(const ParticleSpecies* s) const {
    return (getIndex()==s->getIndex())?this:getHandler()->getMixedObject(getIndex(),s->getIndex());
}
