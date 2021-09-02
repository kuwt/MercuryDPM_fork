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

#ifndef BOUNDARIES_CUBEINSERTIONBOUNDARY_H
#define BOUNDARIES_CUBEINSERTIONBOUNDARY_H

#include "InsertionBoundary.h"
#include "Math/Vector.h"

class BaseParticle;

class RNG;

/*!
 * \class CubeInsertionBoundary
 * \brief It's an insertion boundary which has cuboidal shape (yes, 'CuboidalInsertionBoundary'
 * would have been the correct name).
 */

class CubeInsertionBoundary : public InsertionBoundary
{
public:
    
    /*!
     * \brief Constructor; sets everything to 0.
     */
    CubeInsertionBoundary();
    
    /*!
     * \brief Copy constructor with deep copy.
     */
    CubeInsertionBoundary(const CubeInsertionBoundary& other);
    
    /*!
     * \brief Destructor: default destructor.
     */
    ~CubeInsertionBoundary() override;
    
    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    CubeInsertionBoundary* copy() const override;
    
    /*!
     * \brief Sets the properties of the CubeInsertionBoundary
     */
    void
    set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax,
        double radMin, double radMax);

    void
    set(BaseParticle& particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax,
        double radMin, double radMax);
    /*!
     * \brief Sets the range of particle radii that may be generated.
     */
    void setRadiusRange(Mdouble radMin, Mdouble radMax);
    
    /*!
     * \brief Sets the geometry (position and velocity distribution) of the
     * CubeInsertionBoundary
     */
    void setGeometry(Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax);
    
    /*!
     * \brief Generates a random position, velocity for the particle p
     */
    void placeParticle(BaseParticle* p, RNG& random) override;
    
    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief deprecated version of CubeInsertionBoundary::read().
     */
    MERCURY_DEPRECATED
    void oldRead(std::istream& is);
    
    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;


private:
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

protected:
    
    /*!
     * \brief Minimal and maximal positions defining the boundary's boundaries,
     * and minimum and maximum velocity of the particles to be inserted.
     */
    Vec3D posMin_, posMax_, velMin_, velMax_;
};

#endif
