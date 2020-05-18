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

#include "BidisperseCubeInsertionBoundary.h"

BidisperseCubeInsertionBoundary::BidisperseCubeInsertionBoundary()
{
    particleToCopyA_ = nullptr;
    particleToCopyB_ = nullptr;
    probA_ = 0;
}

BidisperseCubeInsertionBoundary::BidisperseCubeInsertionBoundary(const BidisperseCubeInsertionBoundary& other)
{
    particleToCopyA_ = other.particleToCopyA_;
    particleToCopyB_ = other.particleToCopyB_;
    probA_ = other.probA_;
}

/*!
 * \details Destructor. This class introduces two new pointers (beyond those
 * introduced by its ancestors) so we need to free them. 
 */
BidisperseCubeInsertionBoundary::~BidisperseCubeInsertionBoundary()
{
    delete particleToCopyA_;
    delete particleToCopyB_;
}

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer. 
 * \return      pointer to the copy on the heap
 */
BidisperseCubeInsertionBoundary* BidisperseCubeInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "BidisperseCubeInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new BidisperseCubeInsertionBoundary(*this);
}


void BidisperseCubeInsertionBoundary::set(BaseParticle* particleToCopyA, BaseParticle* particleToCopyB, double probA,
                                          int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax)
{
    particleToCopyA_ = particleToCopyA;
    particleToCopyB_ = particleToCopyB;
    probA_ = probA;
    CubeInsertionBoundary::set(particleToCopyA, maxFailed,
                               posMin, posMax, velMin, velMax, 0, 0);
}

BaseParticle* BidisperseCubeInsertionBoundary::getParticleToCopyA() const
{
    return particleToCopyA_;
}

BaseParticle* BidisperseCubeInsertionBoundary::getParticleToCopyB() const
{
    return particleToCopyB_;
}

/*!
 * \brief Generates a particle with random position and velocity 
 */
BaseParticle* BidisperseCubeInsertionBoundary::generateParticle(RNG& random)
{
    double check = random.getRandomNumber(0, 1);
    BaseParticle* P = (check < probA_) ? particleToCopyA_->copy()
                                       : particleToCopyB_->copy();
    Vec3D pos, vel;
    pos.X = random.getRandomNumber(posMin_.X, posMax_.X);
    pos.Y = random.getRandomNumber(posMin_.Y, posMax_.Y);
    pos.Z = random.getRandomNumber(posMin_.Z, posMax_.Z);
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);
    P->setPosition(pos);
    P->setVelocity(vel);
    return P;
}

/*!
 * \details Reads the boundary properties from an istream
 * \param[in,out] is        the istream
 */
void BidisperseCubeInsertionBoundary::read(std::istream& is)
{
    CubeInsertionBoundary::read(is);
    std::string dummy;
    is >> dummy >> probA_;
}

/*!
 * \details Writes boundary's properties to an ostream
 * \param[in] os    the ostream
 */
void BidisperseCubeInsertionBoundary::write(std::ostream& os) const
{
    CubeInsertionBoundary::write(os);
    os << " probA " << probA_;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'CubeInsertionBoundary'
 */
std::string BidisperseCubeInsertionBoundary::getName() const
{
    return "BidisperseCubeInsertionBoundary";
}

