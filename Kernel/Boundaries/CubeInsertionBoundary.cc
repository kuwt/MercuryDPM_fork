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
#include "CubeInsertionBoundary.h"
#include "Particles/BaseParticle.h"
#include "Math/RNG.h"

/*!
 * \details Default constructor; sets every data member to 0.
 */
CubeInsertionBoundary::CubeInsertionBoundary() : InsertionBoundary()
{
    posMin_ = Vec3D(0.0, 0.0, 0.0);
    posMax_ = Vec3D(0.0, 0.0, 0.0);
}

/*!
 * \details Copy constructor
 */
CubeInsertionBoundary::CubeInsertionBoundary(const CubeInsertionBoundary& other)
        : InsertionBoundary(other)
{
    posMin_ = other.posMin_;
    posMax_ = other.posMax_;
}

/*!
 * \details Destructor. Since there are no pointers in this class, there is no 
 *          need for any actions here.
 */
CubeInsertionBoundary::~CubeInsertionBoundary()
= default;

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer. 
 * \return      pointer to the copy on the heap
 */
CubeInsertionBoundary* CubeInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "CubeInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new CubeInsertionBoundary(*this);
}

/*!
 * \details Sets all the properties of the cuboidal insertion boundary.
 * \param[in] particleToCopy    vector of pointers to the BaseParticles which are used as a basis
 *                              for the particle types (species) to be inserted
 * \param[in] maxFailed         The maximum number of times the insertion of a
 *                              particle may be tried and failed before the insertion
 *                              of particles is considered done
 *                              NB: this property is used in the parent's
 *                              InsertionBoundary::checkBoundaryBeforeTimeStep().
 * \param[in] posMin            First defining corner of cuboidal insertion boundary
 * \param[in] posMax            Second defining corner of cuboidal insertion boundary
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 */
void CubeInsertionBoundary::set(std::vector<BaseParticle*> particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D
posMax, Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
    maxFailed_ = maxFailed;
    setParticleToCopy(particleToCopy);
    setGeometry(posMin, posMax);
}

/*!
 * \details Sets all the properties of the cuboidal insertion boundary.
 * \param[in] particleToCopy    Pointer to the BaseParticle which is used as a basis
 *                              for the particles to be inserted
 * \param[in] maxFailed         The maximum number of times the insertion of a
 *                              particle may be tried and failed before the insertion
 *                              of particles is considered done
 *                              NB: this property is used in the parent's
 *                              InsertionBoundary::checkBoundaryBeforeTimeStep().
 * \param[in] posMin            First defining corner of cuboidal insertion boundary
 * \param[in] posMax            Second defining corner of cuboidal insertion boundary
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 */
void CubeInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D
posMax, Vec3D velMin, Vec3D velMax)
{
    velMin_ = velMin;
    velMax_ = velMax;
    maxFailed_ = maxFailed;
    setParticleToCopy(particleToCopy);
    setGeometry(posMin, posMax);
}

void CubeInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                Vec3D velMin, Vec3D velMax)
{
    set(&particleToCopy, maxFailed, posMin, posMax, velMin, velMax);
}

/*!
 * \details old style set function which also assumes a uniform psd. Note if you want a specific PSD do not use but
 * this is quicker for a uniform in size PSD
 */
void CubeInsertionBoundary::set(BaseParticle* particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                Vec3D velMin, Vec3D velMax, Mdouble rMin, Mdouble rMax)
{
    PSD uniformPSD;
    uniformPSD.setDistributionUniform(rMin, rMax, 2);
    setPSD(uniformPSD);
    set(particleToCopy, maxFailed, posMin, posMax, velMin, velMax);
}

void CubeInsertionBoundary::set(BaseParticle& particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                Vec3D velMin, Vec3D velMax, Mdouble rMin, Mdouble rMax)
{
    PSD uniformPSD;
    uniformPSD.setDistributionUniform(rMin, rMax, 2);
    setPSD(uniformPSD);
    set(particleToCopy, maxFailed, posMin, posMax, velMin, velMax);
}

/*!
 * \details set the geometry of the Cuboidal insertion boundary.
 * \param[in] posMin            First defining corner of cuboidal insertion boundary
 * \param[in] posMax            Second defining corner of cuboidal insertion boundary
 */
void CubeInsertionBoundary::setGeometry(Vec3D posMin, Vec3D posMax)
{
    posMin_ = posMin;
    posMax_ = posMax;
}

/*!
 * \details places a particle with random position (although within the boundary,
 * of course), velocity and radius and returns its pointer.
 * \param[in] random    Random number generator
 */
// Irana: where is the delete of P? There is a new in copy()
// Bram: @Irana: The Particle is only CREATED in the InsertionBoundary. Deletion  
// is done either by a DeletionBoundary, or at the end of a program (in that case, 
// ~ParticleHandler -> ~BaseHandler -> BaseHandler::clear(), where first the individual
// objects are deleted, followed by the clearance of the std::vector with object pointers).
//
// JMFT: @Irana, @Bram: As of 2549, the delete is done in InsertionBoundary.
// The InsertionBoundary calls generateParticle and placeParticle, which are
// geometry-specific and implemented by the child classes, such as
// CubeInsertionBoundary. If the particle is wanted then it is copied into the
// particleHandler.  In any case, the InsertionBoundary's version of the
// particle (created here) gets deleted by the InsertionBoundary.
void CubeInsertionBoundary::placeParticle(BaseParticle* p, RNG& random)
{
    Vec3D pos, vel;
    pos.X = random.getRandomNumber(posMin_.X, posMax_.X);
    pos.Y = random.getRandomNumber(posMin_.Y, posMax_.Y);
    pos.Z = random.getRandomNumber(posMin_.Z, posMax_.Z);
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);
    p->setPosition(pos);
    p->setVelocity(vel);
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream
 */
void CubeInsertionBoundary::read(std::istream& is)
{
    InsertionBoundary::read(is);
    std::string dummy;
    is >> dummy >> posMin_
       >> dummy >> posMax_;
    is >> dummy >> velMin_
       >> dummy >> velMax_;
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void CubeInsertionBoundary::write(std::ostream& os) const
{
    InsertionBoundary::write(os);
    os << " posMin " << posMin_
       << " posMax " << posMax_
       << " velMin " << velMin_
       << " velMax " << velMax_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'CubeInsertionBoundary'
 */
std::string CubeInsertionBoundary::getName() const
{
    return "CubeInsertionBoundary";
}