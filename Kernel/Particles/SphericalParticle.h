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

#ifndef SphericalParticle_H
#define SphericalParticle_H

#include "BaseParticle.h"

/*!
 * \class SphericalParticle
 * \brief A basic particle. 
 * \details Before r3648, this SphericalParticle class was known as BaseParticle. 
 */
class SphericalParticle final : public BaseParticle
{
public:
    
    /*!
     * \brief Particle constructor
     */
    SphericalParticle() = default;
    
    /*!
     * \brief Particle constructor, setting the species.
     * This is useful to create e.g. template particles for the InsertionBoundaries.
     * It is marked explicit because you dont want to use it with the assignment operator, particle = species;
     */
    explicit SphericalParticle(const ParticleSpecies* s) : BaseParticle(s) {
    
    };

    /*!
     * \brief Particle constructor
     */
    SphericalParticle(const SphericalParticle& p) = default;
    
    ~SphericalParticle() override = default;
    
    /*!
     * \brief Particle copy method. Calls copy constructor of this class (useful for polymorphism)
     */
    SphericalParticle* copy() const override {
        return new SphericalParticle(*this);
    }
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override {
        return "BaseParticle";
    }

    bool isSphericalParticle() const override {return true;}
};
#endif
