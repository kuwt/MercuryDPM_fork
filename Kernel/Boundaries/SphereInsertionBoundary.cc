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
#include "SphereInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"

/*!
 * \details Default constructor; sets every data member to 0.
 */
SphereInsertionBoundary::SphereInsertionBoundary()
{
    radMin_ = 0.0;
    radMax_ = 0.0;
    phiMin_ = 0.0;
    phiMax_ = 0.0;
    thetaMin_ = 0.0;
    thetaMax_ = 0.0;
    origin_.setZero();
}

/*!
 * \details Copy constructor
 */
SphereInsertionBoundary::SphereInsertionBoundary(const SphereInsertionBoundary &other)
        : InsertionBoundary(other)
{
    radMin_ = other.radMin_;
    radMax_ = other.radMax_;
    phiMin_ = other.phiMin_;
    phiMax_ = other.phiMax_;
    thetaMin_ = other.thetaMin_;
    thetaMax_ = other.thetaMax_;
    origin_ = other.origin_;
}

/*!
 * \details Destructor. Since there are no pointers in this class, there is no
 *          need for any actions here.
 */
SphereInsertionBoundary::~SphereInsertionBoundary()
= default;

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 * \return      pointer to the copy on the heap
 */
SphereInsertionBoundary* SphereInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "SphereInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new SphereInsertionBoundary(*this);
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
 * \param[in] radMin            minimal radius defining of the sphere boundary
 * \param[in] radMax            maximal radius of the sphere boundary
 * \param[in] phiMin              Minimal polar angle of the sphere boundary
 * \param[in] phiMax              maximal polar angle of the sphere boundary
 * \param[in] thetaMin          minimal azimuth angle of the sphere boundary
 * \param[in] thetaMax          maximal azimuth angle of the sphere boundary
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 */
void SphereInsertionBoundary::set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax,
                                    Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
    maxFailed_ = maxFailed;
    setParticleToCopy(particleToCopy);
    setGeometry(radMin, radMax, thetaMin, thetaMax, phiMin, phiMax);
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
 * \param[in] radMin            minimal radius defining of the sphere boundary
 * \param[in] radMax            maximal radius of the sphere boundary
 * \param[in] phiMin              Minimal polar angle of the sphere boundary
 * \param[in] phiMax              maximal polar angle of the sphere boundary
 * \param[in] thetaMin          minimal azimuth angle of the sphere boundary
 * \param[in] thetaMax          maximal azimuth angle of the sphere boundary
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 */
void SphereInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax,
                                    Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
    maxFailed_ = maxFailed;
    setParticleToCopy(particleToCopy);
    setGeometry(radMin, radMax, thetaMin, thetaMax, phiMin, phiMax);
}

void SphereInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax,
                                    Vec3D velMin, Vec3D velMax)
{
    set(&particleToCopy, maxFailed, radMin, radMax, thetaMin, thetaMax, phiMin, phiMax, velMin, velMax);
}

/*!
 * \details old style set function which also assumes a uniform psd. Note if you want a specific PSD do not use but
 * this is quicker for a uniform in size PSD
 */
void SphereInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax,
                                    Vec3D velMin, Vec3D velMax, Mdouble particleRMin, Mdouble particleRMax)
{
    PSD uniformPSD;
    uniformPSD.setDistributionUniform(particleRMin, particleRMax, 2);
    setPSD(uniformPSD);
    set(particleToCopy, maxFailed, radMin, radMax, thetaMin, thetaMax, phiMin, phiMax, velMin, velMax);
}

void SphereInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Mdouble radMin,
                                    Mdouble radMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax,
                                    Vec3D velMin, Vec3D velMax, Mdouble particleRMin, Mdouble particleRMax)
{
    PSD uniformPSD;
    uniformPSD.setDistributionUniform(particleRMin, particleRMax, 2);
    setPSD(uniformPSD);
    set(particleToCopy, maxFailed, radMin, radMax, thetaMin, thetaMax, phiMin, phiMax, velMin, velMax);
}

/*!
 * \details set the geometry of the Cylindrical insertion boundary by setting the position and orientation.
 * \param[in] radMin            minimal radius defining of the sphere boundary
 * \param[in] radMax            maximal radius of the sphere boundary
 * \param[in] phiMin            minimal polar angle of the sphere boundary
 * \param[in] phiMax            maximal polar angle of the sphere boundary
 * \param[in] thetaMin          minimal azimuth angle of the sphere boundary
 * \param[in] thetaMax          maximal azimuth angle of the sphere boundary
 */
void SphereInsertionBoundary::setGeometry(Mdouble rMin, Mdouble rMax, Mdouble thetaMin, Mdouble thetaMax, Mdouble phiMin, Mdouble phiMax)
{
    radMin_ = rMin;
    radMax_ = rMax;
    phiMin_ = phiMin;
    phiMax_ = phiMax;
    thetaMin_ = thetaMin;
    thetaMax_ = thetaMax;
}

/*!
 * \details places a particle with random position (although within the boundary,
 * of course), velocity and radius and returns its pointer.
 * \param[in] random    Random number generator
 */
void SphereInsertionBoundary::placeParticle(BaseParticle* p, RNG& random)
{
    // create vectors for cartesian position and velocity and cylindrical position
    Vec3D pos, vel, posSpherical;

    // draw random numbers in cartesian coordinates and make sure they are inside the sphere sphere
    do
    {
        pos.X = random.getRandomNumber(-radMax_,radMax_);
        pos.Y = random.getRandomNumber(-radMax_,radMax_);
        pos.Z = random.getRandomNumber(-radMax_,radMax_);
        posSpherical = pos.getSphericalCoordinates();
    }
        // make sure particle stay inside the sphere of our sphere
    while (posSpherical.X > radMax_ ||
           posSpherical.X < radMin_ ||
           posSpherical.Y < phiMin_ ||
           posSpherical.Y >= phiMax_ ||
           posSpherical.Z < thetaMin_||
           posSpherical.Z >= thetaMax_);

    // set the velocity in cartesian coordinates
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);

    p->setPosition(pos + origin_);
    p->setVelocity(vel);
}

/*!
 * \details shifts the sphere by a translation vector in the cartesian coordinate system
 * \param[in] shift        translation vector which shifts the sphere by certain X,Y,Z coordinates
 */
void SphereInsertionBoundary::shiftBoundary(Vec3D shift)
{
    origin_ += shift;
}

/*!
 * \details Returns the origin of the spheres coordinate system
 * \return origin_      the origin of the spheres coordinate system
 */
Vec3D SphereInsertionBoundary::getOrigin() const
{
    return origin_;
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream
 */
void SphereInsertionBoundary::read(std::istream& is)
{
    InsertionBoundary::read(is);
    std::string dummy;
    is >> dummy >> radMin_
       >> dummy >> radMax_;
    is >> dummy >> phiMin_
       >> dummy >> phiMax_;
    is >> dummy >> thetaMin_
       >> dummy >> thetaMax_;
    is >> dummy >> velMin_
       >> dummy >> velMax_;
    is >> dummy >> origin_;
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void SphereInsertionBoundary::write(std::ostream& os) const
{
    InsertionBoundary::write(os);
    os << " radMin " << radMin_
       << " radMax " << radMax_
       << " phiMin " << phiMin_
       << " phiMax " << phiMax_
       << " thetaMin " << thetaMin_
       << " thetaMax " << thetaMax_
       << " velMin " << velMin_
       << " velMax " << velMax_
       << " origin " << origin_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'SphereInsertionBoundary'
 */
std::string SphereInsertionBoundary::getName() const
{
    return "SphereInsertionBoundary";
}

