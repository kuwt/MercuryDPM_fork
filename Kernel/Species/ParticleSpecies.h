//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef ParticleSpecies_H
#define ParticleSpecies_H

#include "BaseSpecies.h"
#include "Math/ExtendedMath.h"
#include "Interactions/BaseInteraction.h"

class BaseInteractable;

//Note the getVelocity can for some Species be dependent on which point on the Species is meant.
class ParticleSpecies : public BaseSpecies
{
public:
    typedef BaseInteraction InteractionType;
    typedef BaseSpecies MixedSpeciesType;
    
    ///\brief The default constructor.
    ParticleSpecies();

    ///\brief The default copy constructor.
    ParticleSpecies(const ParticleSpecies& p);

    ParticleSpecies(BaseNormalForce* normalForce, BaseFrictionForce* frictionForce, BaseAdhesiveForce* adhesiveForce);

    ///\brief The default destructor.
    ~ParticleSpecies();
    
    /*!
     * \brief Creates a deep copy of the object from which it is called.
     * \details See BaseSpecies::copy for details
     */
    ParticleSpecies* copy() const override = 0;
    
    /*!
     * \brief Creates a new MixedSpecies with the same force properties as the 
     * Species from which it is called. See Species::copyMixed for details
     */
    virtual BaseSpecies* copyMixed() const = 0;
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is) override;
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const override;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;
    
    ///Allows density_ to be changed
    ///\todo recalculate masses when setting dim_particle or rho
    void setDensity(Mdouble density);
    
    ///\todo TW: should getMassFromRadius be removed? IFCD: it is used in at 
    /// least one driver (AxisymmetricHopper).
    Mdouble getMassFromRadius(Mdouble radius) const;

    Mdouble getMassFromRadius(const Mdouble radius,SpeciesHandler& speciesHandler);

    Mdouble getVolumeFromRadius(Mdouble radius) const;
    
    ///Allows density_ to be accessed
    Mdouble getDensity() const;
    
    /*!
     * \brief Compute Particle mass function, which required a reference to the Species vector. It computes the Particles mass, Inertia and the inverses.
     */
    void computeMass(BaseParticle* p) const;
    
    void setTemperatureDependentDensity(const std::function<double(double)>& temperatureDependentDensity);
    
    const std::function<double(double)>& getTemperatureDependentDensity() const;
    
    /*!
     * \brief Computes mass of the lightest particle (by mass) belonging to this species.
     * This computation calls getLightestInverseParticleMassLocal, such that the computation is done on each node.
     */
    Mdouble getSmallestParticleMass() const;
    
    ///\brief returns the largest separation distance at which adhesive short-range forces can occur.
    /*!
     * \details returns the largest separation distance (negative overlap) at which
     * (adhesive) short-range forces can occur (needed for contact detection).
     * Defined in each of the AdhesiveForceSpecies
     * It is defined as a virtual function here to allow the function
     * to be called from a BaseSpecies pointer (which is the kind of pointer
     * used for MixedSpecies).
     */
    Mdouble getMaxInteractionDistance() const {return maxInteractionDistance_;}
    
    /**
     * \brief Sets #maxInteractionDistance_
     * @param interactionDistance the interaction distance that has been changed
     */
    void setMaxInteractionDistance(Mdouble interactionDistance=0);
    
    const BaseSpecies* getMixedSpecies(const ParticleSpecies* s) const;
    
private:
    
    /*!
     * \brief Computes inverse mass of the lightest particle (by mass) belonging to this species.
     * If MPI is used, this computation is done locally on each node.
     */
    Mdouble getLargestInverseParticleMassLocal() const;

    /*!
     * \brief The mass density.
     */
    Mdouble density_;
    
    /*!
     * Change this function to let the particles expand due to temperature.
     * The default value (empty) stands for constant density.
     */
    std::function<double(double temperature)> temperatureDependentDensity_;
    
    /**
     * Returns the max distance between particles of this species and any other species below which adhesive forces can occur (needed for contact detection)
     */
    Mdouble maxInteractionDistance_;
    
    
};

#endif
