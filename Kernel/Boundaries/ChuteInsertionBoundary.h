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

#ifndef BOUNDARIES_CHUTEINSERTIONBOUNDARY_H
#define BOUNDARIES_CHUTEINSERTIONBOUNDARY_H

#include "InsertionBoundary.h"
#include "Math/Vector.h"

class BaseParticle;

class RNG;

/*!
 * \class ChuteInsertionBoundary
 * \brief Used for modeling chute inflow. Inherits from InsertionBoundary.
 */
class ChuteInsertionBoundary : public InsertionBoundary
{
public:
    
    /*!
     * \brief Default constructor. 
     */
    ChuteInsertionBoundary();
    
    ///\brief Copy constructor.
    ChuteInsertionBoundary(const ChuteInsertionBoundary& other);
    
    /*!
     * \brief Copy method; creates a copy on the heap.
     */
    ChuteInsertionBoundary* copy() const override;
    
    /*!
     * \brief Sets all boundary properties at once.
     */
    void set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, double
    fixedParticleRadius, double inflowVelocity, double inflowVelocityVariance);
    
    /*!
     * \brief Sets all boundary properties at once.
     */
    void set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, double
    fixedParticleRadius, double inflowVelocity, double inflowVelocityVariance);
    
    /*!
     * \brief old style set function which assumes a uniform psd.
     */
    void set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, Mdouble rMin,
             Mdouble rMax, double fixedParticleRadius, double inflowVelocity, double inflowVelocityVariance);
    
    void placeParticle(BaseParticle* p, RNG& random) override;
    
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
    
    /*!
     * \brief The two extremal corners of the cuboidal insertion boundary
     */
    Vec3D posMin_, posMax_;
    
    /*!
     * \brief radius of the fixed bottom particles, mean particle velocity in 
     * X-direction, and allowed maximum randomly added/subtracted velocities in all three
     * directions while generating particles (expressed as \% of inflowVelocity_).
     * NB: see also documentation of ChuteInsertionBoundary::generateParticle().
     */
    double fixedParticleRadius_, inflowVelocity_, inflowVelocityVariance_;
    
};

#endif
