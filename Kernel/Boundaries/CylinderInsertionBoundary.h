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

#ifndef BOUNDARIES_CYLINDERINSERTIONBOUNDARY_H
#define BOUNDARIES_CYLINDERINSERTIONBOUNDARY_H

#include "InsertionBoundary.h"
#include "Math/Vector.h"

class BaseParticle;

class RNG;

/*!
 * \class CylinderInsertionBoundary
 * \brief Inserts particles in a cylinder with a given radius, height and azimuthal radius at a specific origin which is 0 by default.
 */

class CylinderInsertionBoundary : public InsertionBoundary
{
public:

    /*!
     * \brief Constructor; sets everything to 0.
     */
    CylinderInsertionBoundary();

    /*!
     * \brief Copy constructor with deep copy.
     */
    CylinderInsertionBoundary(const CylinderInsertionBoundary &other);

    /*!
     * \brief Destructor: default destructor.
     */
    ~CylinderInsertionBoundary() override;

    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    CylinderInsertionBoundary *copy() const override;

    /*!
     * \brief Sets the properties of the InsertionBoundary for mutliple different particle types
     */
    void set(BaseParticle *particleToCopy, unsigned int maxFailed, Mdouble rMin, Mdouble rMax, Mdouble hMin, Mdouble hMax,
             Vec3D normal = {0, 0, 1}, Mdouble thetaMin = -constants::pi, Mdouble thetaMax = constants::pi, Vec3D velMin = {0, 0, 0}, Vec3D velMax = {0, 0, 0});

    void
    set(BaseParticle &particleToCopy, unsigned int maxFailed, Mdouble rMin, Mdouble rMax, Mdouble hMin, Mdouble hMax,
        Vec3D normal = {0, 0, 1}, Mdouble thetaMin = -constants::pi, Mdouble thetaMax = constants::pi, Vec3D velMin = {0, 0, 0}, Vec3D velMax = {0, 0, 0});

    /*!
     * \brief Sets the properties of the InsertionBoundary for a single particle type
     */
    void
    set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Mdouble rMin,
        Mdouble rMax, Mdouble hMin, Mdouble hMax, Vec3D normal = {0, 0, 1}, Mdouble thetaMin = -constants::pi, Mdouble thetaMax = constants::pi, Vec3D velMin = {0, 0, 0},
        Vec3D velMax = {0,0,0});

    /*!
     * \brief old style set function which assumes a uniform psd.
     */
    void set(BaseParticle* particleToCopy, unsigned int maxFailed, Mdouble radMin,
             Mdouble radMax, Mdouble hMin, Mdouble hMax, Vec3D normal, Mdouble thetaMin, Mdouble thetaMax, Vec3D velMin, Vec3D velMax,
             Mdouble particleRMin, Mdouble particleRMax);

    /*!
     * \brief old style set function which assumes a uniform psd.
     */
    void set(BaseParticle& particleToCopy, unsigned int maxFailed, Mdouble radMin, Mdouble radMax, Mdouble hMin,
             Mdouble hMax, Vec3D normal, Mdouble thetaMin, Mdouble thetaMax, Vec3D velMin, Vec3D velMax, Mdouble particleRMin,
             Mdouble particleRMax);
    
    /*!
     * \brief Sets the geometry (position and orientation) of the
     * CylinderInsertionBoundary
     */
    void setGeometry(Mdouble rMin, Mdouble rMax, Mdouble hMin, Mdouble hMax, Mdouble thetaMin, Mdouble thetaMax, Vec3D normal);

    /*!
     *
     * \brief sets the orientation of the Cylinder
     */
    void setOrientation(Vec3D normal);

    /*!
     * \brief gets the rotationMatrix for a rotation around the X axis by a certain angle
     */
    Matrix3D getRotationMatrixX(Mdouble rotationAngle);

    /*!
     * \brief gets the rotationMatrix for a rotation around the Y axis by a certain angle
     */
    Matrix3D getRotationMatrixY(Mdouble rotationAngle);

    /*!
     * \brief gets the rotationMatrix for a rotation around the Z axis by a certain angle
     */
    Matrix3D getRotationMatrixZ(Mdouble rotationAngle);

    /*!
     * \brief Returns the origin of the cylinders coordinate system
     */
    Vec3D getOrigin() const;

    /*!
     * \brief Returns the orientation matrix of the cylinder
     */
    Matrix3D getOrientationMatrix() const;

    /*!
     * \brief shift the origin of the cylinder in the cartesian coordinate system
     */
    void shiftBoundary(Vec3D shift) override;

    /*!
     * \brief rotate the cylinder around its origin
     */
    void rotateBoundary(Vec3D angle) override;

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
     * \brief Minimal and maximal azimuth angle in range [-pi,pi] defining the circle of the cylinder in which insertion should happen.
     */
    Mdouble thetaMin_, thetaMax_;

    /*!
     * \brief Minimal and maximal radius defining the inner and outer circle of the cylinder in which insertion should happen.
     */
    Mdouble radMin_, radMax_;

    /*!
     * \brief Minimal and maximal height defining the cylinder in which insertion should happen.
     */
    Mdouble hMin_, hMax_;

    /*!
     * \brief Orientation matrix which rotates the cylinder based on the chosen normal axis (the height axis of the cylinder).
     */
    Matrix3D orientationMatrix_;

    /*!
     * \brief origin of the cylinder.
     */
    Vec3D origin_;
};

#endif
