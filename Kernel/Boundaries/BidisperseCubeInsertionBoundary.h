//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
