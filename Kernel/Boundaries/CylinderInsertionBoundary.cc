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

#include <random>
#include "CylinderInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"

/*!
 * \details Default constructor; sets every data member to 0.
 */
CylinderInsertionBoundary::CylinderInsertionBoundary() : InsertionBoundary()
{
    radMin_ = 0.0;
    radMax_ = 0.0;
    hMin_ = 0.0;
    hMax_ = 0.0;
    thetaMin_ = 0.0;
    thetaMax_ = 0.0;
    origin_.setZero();
    orientationMatrix_.setZero();
}

/*!
 * \details Copy constructor
 */
CylinderInsertionBoundary::CylinderInsertionBoundary(const CylinderInsertionBoundary& other)
        : InsertionBoundary(other)
{
    radMin_ = other.radMin_;
    radMax_ = other.radMax_;
    hMin_ = other.hMin_;
    hMax_ = other.hMax_;
    thetaMin_ = other.thetaMin_;
    thetaMax_ = other.thetaMax_;
    origin_ = other.origin_;
    orientationMatrix_ = other.orientationMatrix_;
}

/*!
 * \details Destructor. Since there are no pointers in this class, there is no 
 *          need for any actions here.
 */
CylinderInsertionBoundary::~CylinderInsertionBoundary()
= default;

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer. 
 * \return      pointer to the copy on the heap
 */
CylinderInsertionBoundary* CylinderInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "CylinderInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new CylinderInsertionBoundary(*this);
}

/*!
 * \details Sets all the properties of the cylindrical insertion boundary.
 * \param[in] particleToCopy    vector of pointers to the BaseParticles which are used as a basis
 *                              for the particle types (species) to be inserted
 * \param[in] maxFailed         The maximum number of times the insertion of a
 *                              particle may be tried and failed before the insertion
 *                              of particles is considered done
 *                              NB: this property is used in the parent's
 *                              InsertionBoundary::checkBoundaryBeforeTimeStep().
 * \param[in] radMin            minimal radius defining of the cylinder boundary
 * \param[in] radMax            maximal radius of the cylinder boundary
 * \param[in] hMin              minimal height of the cylinder boundary
 * \param[in] hMax              maximal height of the cylinder boundary
 * \param[in] thetaMin          minimal azimuth angle of the cylinder boundary
 * \param[in] thetaMax          maximal azimuth angle of the cylinder boundary
 * \param[in] normal            normal position of the Z-axis (height) of the cylinder
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 */
void CylinderInsertionBoundary::set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble hMin, Mdouble hMax, Vec3D normal, Mdouble thetaMin, Mdouble thetaMax,
                                    Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
    maxFailed_ = maxFailed;
    setParticleToCopy(particleToCopy);
    setGeometry(radMin, radMax, hMin, hMax, thetaMin, thetaMax, normal);
}

/*!
 * \details Sets all the properties of the cylindrical insertion boundary.
 * \param[in] particleToCopy    Pointer to the BaseParticle which is used as a basis
 *                              for the particles to be inserted
 * \param[in] maxFailed         The maximum number of times the insertion of a
 *                              particle may be tried and failed before the insertion
 *                              of particles is considered done
 *                              NB: this property is used in the parent's
 *                              InsertionBoundary::checkBoundaryBeforeTimeStep().
 * \param[in] radMin            minimal radius defining of the cylinder boundary
 * \param[in] radMax            maximal radius of the cylinder boundary
 * \param[in] hMin              minimal height of the cylinder boundary
 * \param[in] hMax              maximal height of the cylinder boundary
 * \param[in] thetaMin          minimal azimuth angle of the cylinder boundary
 * \param[in] thetaMax          maximal azimuth angle of the cylinder boundary
 * \param[in] normal            normal position of the Z-axis (height) of the cylinder
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 */
void CylinderInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble hMin, Mdouble hMax, Vec3D normal, Mdouble thetaMin,
                                    Mdouble thetaMax, Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
    maxFailed_ = maxFailed;
    setParticleToCopy(particleToCopy);
    setGeometry(radMin, radMax, hMin, hMax, thetaMin, thetaMax, normal);
}

void CylinderInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble hMin, Mdouble hMax, Vec3D normal, Mdouble thetaMin, Mdouble thetaMax,
                                    Vec3D velMin, Vec3D velMax)
{
    set(&particleToCopy, maxFailed, radMin, radMax, hMin, hMax, normal, thetaMin, thetaMax, velMin, velMax);
}

/*!
 * \details old style set function which also assumes a uniform psd. Note if you want a specific PSD do not use but
 * this is quicker for a uniform in size PSD
 */
void CylinderInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble hMin, Mdouble hMax, Vec3D normal, Mdouble thetaMin, Mdouble thetaMax,
                                    Vec3D velMin, Vec3D velMax, Mdouble particleRMin, Mdouble particleRMax)
{
    PSD uniformPSD;
    uniformPSD.setDistributionUniform(particleRMin, particleRMax, 2);
    setPSD(uniformPSD);
    set(particleToCopy, maxFailed, radMin, radMax, hMin, hMax, normal, thetaMin, thetaMax, velMin, velMax);
}

void CylinderInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble hMin, Mdouble hMax, Vec3D normal, Mdouble thetaMin, Mdouble thetaMax,
                                    Vec3D velMin, Vec3D velMax, Mdouble particleRMin, Mdouble particleRMax)
{
    PSD uniformPSD;
    uniformPSD.setDistributionUniform(particleRMin, particleRMax, 2);
    setPSD(uniformPSD);
    set(particleToCopy, maxFailed, radMin, radMax, hMin, hMax, normal, thetaMin, thetaMax, velMin, velMax);
}

/*!
 * \details set the geometry of the Cylindrical insertion boundary by setting the position and orientation.
 * \param[in] rMin              minimal radius defining an inner radius of the cylinder boundary
 * \param[in] rMax              maximal outer radius of the cylinder boundary
 * \param[in] hMin              minimal height of the cylinder boundary
 * \param[in] hMax              maximal height of the cylinder boundary
 * \param[in] thetaMin          minimal azimuth angle of the cylinder boundary
 * \param[in] thetaMax          maximal azimuth angle of the cylinder boundary
 * \param[in] normal            normal position of the Z-axis (height) of the cylinder
 */
void CylinderInsertionBoundary::setGeometry(Mdouble rMin, Mdouble rMax, Mdouble hMin, Mdouble hMax, Mdouble thetaMin, Mdouble thetaMax, Vec3D normal)
{
    radMin_ = rMin;
    radMax_ = rMax;
    hMin_ = hMin;
    hMax_ = hMax;
    thetaMin_ = thetaMin;
    thetaMax_ = thetaMax;
    setOrientation(normal);
}


/*!
 * \details set the orientation of the cylinder by rotating it around the X, Y and Z axis. The order of rotation is first X, then Y, then Z.
 * \param[in] normal            Defines the angle of rotation in rad
 */
void CylinderInsertionBoundary::setOrientation(Vec3D normal)
{
    logger.assert_always(!normal.isZero(),"no normal vector defined. Please define a normal vector to determine the coordinate system for the cylinder");
    logger.assert_always(normal.X == 0 || normal.X == 1,"the normal in X direction has to be either set to 0 or 1");
    logger.assert_always(normal.Y == 0 || normal.Y == 1,"the normal in Y direction has to be either set to 0 or 1");
    logger.assert_always(normal.Z == 0 || normal.Z == 1,"the normal in Z direction has to be either set to 0 or 1");
    logger.assert_always(normal.getLength() == 1, "The normal can only point to a single direction. Please set the normal to point to the X, Y or Z coordinate.");
    if (normal.Z)
    {
        // do nothing
        orientationMatrix_ = getRotationMatrixZ(0);
    }
    if (normal.Y)
    {
        // rotate around X to get the Z coordinate pointing to Y
        orientationMatrix_ = getRotationMatrixX(constants::pi/2.0);
    }
    if (normal.X)
    {
        // rotate around Y to get the Z coordinate pointing to X
        orientationMatrix_ = getRotationMatrixY(constants::pi/2.0);
    }
}


/*!
 * \details gets the rotationMatrix for a rotation around the X axis by a certain angle
 * \param[in] rotationAngle  the angle by which the rotation should happen
 * \return the rotation matrix in X direction
 */
Matrix3D CylinderInsertionBoundary::getRotationMatrixX(Mdouble rotationAngle) {
    Matrix3D rotationMatrix;
    rotationMatrix.XX = 1.0;
    rotationMatrix.YY = cos(rotationAngle);
    rotationMatrix.YZ = -sin(rotationAngle);
    rotationMatrix.ZY = sin(rotationAngle);
    rotationMatrix.ZZ = cos(rotationAngle);
    return rotationMatrix;
}

/*!
 * \details gets the rotationMatrix for a rotation around the Y axis by a certain angle
 * \param[in] rotationAngle  the angle by which the rotation should happen
 * \return the rotation matrix in Y direction
 */
Matrix3D CylinderInsertionBoundary::getRotationMatrixY(Mdouble rotationAngle) {
    Matrix3D rotationMatrix;
    rotationMatrix.XX = cos(rotationAngle);
    rotationMatrix.XZ = sin(rotationAngle);
    rotationMatrix.YY = 1.0;
    rotationMatrix.ZX = -sin(rotationAngle);
    rotationMatrix.ZZ = cos(rotationAngle);
    return rotationMatrix;
}

/*!
 * \details gets the rotationMatrix for a rotation around the Z axis by a certain angle
 * \param[in] rotationAngle  the angle by which the rotation should happen
 * \return the rotation matrix in Z direction
 */
Matrix3D CylinderInsertionBoundary::getRotationMatrixZ(Mdouble rotationAngle) {
    Matrix3D rotationMatrix;
    rotationMatrix.XX = cos(rotationAngle);
    rotationMatrix.XY = -sin(rotationAngle);
    rotationMatrix.YX = sin(rotationAngle);
    rotationMatrix.YY = cos(rotationAngle);
    rotationMatrix.ZZ = 1.0;
    return rotationMatrix;
}


/*!
 * \details rotate the cylinder around its origin by an XYZ rotation. The cylinder first gets rotated around X axis by
 * angle.X then around the Y axis by angle.Y and finally around the Z axis by angle.Z.
 * \param[in] angle        Defines the angles of rotation in rad
 */
void CylinderInsertionBoundary::rotateBoundary(Vec3D angle)
{
    orientationMatrix_ = getRotationMatrixZ(angle.Z) * getRotationMatrixY(angle.Y) * getRotationMatrixX(angle.X)*orientationMatrix_;
}


/*!
 * \details places a particle with random position (although within the boundary,
 * of course), velocity and radius and returns its pointer.
 * \param[in] random    Random number generator
 */
void CylinderInsertionBoundary::placeParticle(BaseParticle* p, RNG& random)
{
    // create vectors for cartesian position and velocity and cylindrical position
    Vec3D pos, vel, posCylindrical;

    // draw random numbers in cartesian coordinates and make sure they are inside the cylinder sphere
    do
    {
        pos.X = random.getRandomNumber(-radMax_,radMax_);
        pos.Y = random.getRandomNumber(-radMax_,radMax_);
        pos.Z = random.getRandomNumber(hMin_, hMax_);
        posCylindrical = pos.getCylindricalCoordinates();
    }
    // make sure particle stay inside the sphere of our cylinder
    while (posCylindrical.X > radMax_ ||
           posCylindrical.X < radMin_ ||
           posCylindrical.Y < thetaMin_ ||
           posCylindrical.Y >= thetaMax_);

    // rotate cylindrical coordinates by the orientationMatrix
    pos = orientationMatrix_*pos;

    // set the velocity in cartesian coordinates
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);

    p->setPosition(pos + origin_);
    p->setVelocity(vel);
}

/*!
 * \details shifts the cylinder by a translation vector in the cartesian coordinate system
 * \param[in] shift        translation vector which shifts the cylinder by certain X,Y,Z coordinates
 */
void CylinderInsertionBoundary::shiftBoundary(Vec3D shift)
{
    origin_ += shift;
}

/*!
 * \details Returns the origin of the cylinders coordinate system
 * \return origin_      the origin of the cylinders coordinate system
 */
Vec3D CylinderInsertionBoundary::getOrigin() const
{
    return origin_;
}

/*!
 * \details Returns the orientation matrix of the cylinder
 * \return orientationMatrix_  the orientation matrix which rotates the particles inserted into the cylinder to the right spot
 */
Matrix3D CylinderInsertionBoundary::getOrientationMatrix() const
{
    return orientationMatrix_;
}


/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream
 */
void CylinderInsertionBoundary::read(std::istream& is)
{
    InsertionBoundary::read(is);
    std::string dummy;
    is >> dummy >> radMin_
       >> dummy >> radMax_;
    is >> dummy >> hMin_
       >> dummy >> hMax_;
    is >> dummy >> thetaMin_
       >> dummy >> thetaMax_;
    is >> dummy >> orientationMatrix_;
    is >> dummy >> velMin_
       >> dummy >> velMax_;
    is >> dummy >> origin_;
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void CylinderInsertionBoundary::write(std::ostream& os) const
{
    InsertionBoundary::write(os);
    os << " radMin " << radMin_
       << " radMax " << radMax_
       << " hMin " << hMin_
       << " hMax " << hMax_
       << " thetaMin " << thetaMin_
       << " thetaMax " << thetaMax_
       << " orientationMatrix " << orientationMatrix_
       << " velMin " << velMin_
       << " velMax " << velMax_
       << " origin " << origin_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'CylinderInsertionBoundary'
 */
std::string CylinderInsertionBoundary::getName() const
{
    return "CylinderInsertionBoundary";
}
