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

#ifndef MERCURYDPM_SPHEREINSERTIONBOUNDARY_H
#define MERCURYDPM_SPHEREINSERTIONBOUNDARY_H

#include "InsertionBoundary.h"
#include "Math/Vector.h"

/*!
 * \class SphereInsertionBoundary
 * \brief Inserts particles in a sphere with a given radius, polar and azimuthal radius at a specific origin which is 0 by default.
 */
class SphereInsertionBoundary : public InsertionBoundary
{
public:

    /*!
     * \brief Constructor; sets everything to 0.
     */
    SphereInsertionBoundary();

    /*!
     * \brief Copy constructor with deep copy.
     */
    SphereInsertionBoundary(const SphereInsertionBoundary &other);

    /*!
     * \brief Destructor: default destructor.
     */
    ~SphereInsertionBoundary() override;

    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    SphereInsertionBoundary *copy() const override;

    /*!
     * \brief Sets the properties of the InsertionBoundary for mutliple different particle types
     */
    void set(BaseParticle *particleToCopy, unsigned int maxFailed, Mdouble rMin, Mdouble rMax, Mdouble thetaMin = 0,
             Mdouble thetaMax = 2.0*constants::pi, Mdouble phiMin = 0, Mdouble phiMax = constants::pi,
             Vec3D velMin = {0, 0, 0}, Vec3D velMax = {0, 0, 0});

    void
    set(BaseParticle &particleToCopy, unsigned int maxFailed, Mdouble rMin, Mdouble rMax, Mdouble thetaMin = 0,
        Mdouble thetaMax = 2.0*constants::pi, Mdouble phiMin = 0, Mdouble phiMax = constants::pi,
        Vec3D velMin = {0, 0, 0}, Vec3D velMax = {0, 0, 0});

    /*!
     * \brief Sets the properties of the InsertionBoundary for a single particle type
     */
    void
    set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Mdouble rMin, Mdouble rMax, Mdouble thetaMin = 0,
        Mdouble thetaMax = 2.0*constants::pi, Mdouble phiMin = 0, Mdouble phiMax = constants::pi,
        Vec3D velMin = {0, 0, 0}, Vec3D velMax = {0,0,0});

    /*!
     * \brief old style set function which assumes a uniform psd.
     */
    void set(BaseParticle* particleToCopy, unsigned int maxFailed, Mdouble radMin,
             Mdouble radMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax, Vec3D velMin, Vec3D velMax,
             Mdouble particleRMin, Mdouble particleRMax);

    /*!
     * \brief old style set function which assumes a uniform psd.
     */
    void set(BaseParticle& particleToCopy, unsigned int maxFailed, Mdouble radMin, Mdouble radMax, Mdouble thetaMin,
             Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax, Vec3D velMin, Vec3D velMax, Mdouble particleRMin,
             Mdouble particleRMax);

    /*!
     * \brief Sets the geometry (position and orientation) of the
     * SphereInsertionBoundary
     */
    void setGeometry(Mdouble rMin, Mdouble rMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax);

    /*!
     * \brief Returns the origin of the spheres coordinate system
     */
    Vec3D getOrigin() const;

    /*!
     * \brief shift the origin of the sphere in the cartesian coordinate system
     */
    void shiftBoundary(Vec3D shift) override;

    /*!
     * \brief Generates a random position, velocity for the particle p
     */
    void placeParticle(BaseParticle* p, RNG& random) override;

    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;

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
     * \brief Minimal and maximal azimuth angle in range [0,2*pi] defining the circle of the sphere in which insertion should happen.
     */
    Mdouble thetaMin_, thetaMax_;

    /*!
     * \brief Minimal and maximal polar angle in range [0,pi] defining the circle of the sphere in which insertion should happen.
     */
    Mdouble phiMin_, phiMax_;
    /*!
     * \brief Minimal and maximal radius defining the inner and outer circle of the sphere in which insertion should happen.
     */
    Mdouble radMin_, radMax_;

    /*!
     * \brief origin of the sphere.
     */
    Vec3D origin_;
};


#endif //MERCURYDPM_SPHEREINSERTIONBOUNDARY_H
