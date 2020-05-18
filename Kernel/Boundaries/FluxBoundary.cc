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

#include <limits>

#include "FluxBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"

/*!
 * \details Default constructor (calls the parent-constructor of BaseBoundary as well)
 */
FluxBoundary::FluxBoundary()
        : BaseBoundary()
{
    distance_ = std::numeric_limits<double>::quiet_NaN();
    scaleFactor_ = std::numeric_limits<double>::quiet_NaN();
    numberOfParticlesCrossedForw_ = 0;
    numberOfParticlesCrossedBack_ = 0;
    massCrossedForw_ = 0;
    massCrossedBack_ = 0;
    volumeCrossedForw_ = 0;
    volumeCrossedBack_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"FluxBoundary::FluxBoundary() finished"<<std::endl;
#endif
}

/*!
 * \details Destructor
 */
FluxBoundary::~FluxBoundary()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"FluxBoundary::~FluxBoundary() finished"<<std::endl;
#endif
}

/*!
 * \details Copy function, which creates a copy and returns a pointer to that copy
 * (on the heap)
 * \return pointer to the copy
 */
FluxBoundary* FluxBoundary::copy() const
{
    return new FluxBoundary(*this);
}

/*!
 * \details Defines the placing of the (2D) boundary based on the given normal 
 * and distance.
 * \param[in] normal boundary normal vector
 * \param[in] distance 'distance' between the origin and the boundary,
 * such that the following relation is satisfied:
 * \f[     
 * \mathbf{r} \cdot \mathbf{\hat{n}} = d
 * \f]
 * in which \f$ \mathbf{\hat{n}} \f$ and \f$ d \f$ are the given normal vector and
 * distance, respectively. 
 * NB: If the distance is the ACTUAL distance from the origin, the normal vector 
 * must be of UNIT LENGTH for the placing of the boundary to be done correctly.
 */
void FluxBoundary::set(const Vec3D& normal, Mdouble distance)
{
    scaleFactor_ = 1. / std::sqrt(Vec3D::dot(normal, normal));
    normal_ = normal * scaleFactor_;
    distance_ = distance * scaleFactor_;
}

/*!
 * \details Resets the various counts to zero.
 */
void FluxBoundary::reset()
{
    numberOfParticlesCrossedForw_ = 0;
    numberOfParticlesCrossedBack_ = 0;
    massCrossedForw_ = 0;
    massCrossedBack_ = 0;
    volumeCrossedForw_ = 0;
    volumeCrossedBack_ = 0;
}

/*!
 * \details Resets the boundary's 'distance' from the origin to be the one given.
 * \param[in] distance the new 'distance' between boundary and origin.
 * see also comments of FluxBoundary::set().
 */
void FluxBoundary::move(Mdouble distance)
{
    distance_ = distance * scaleFactor_;
}

/*!
 * \details Calculates the shortest distance between the wall and given position.
 * \param[in] position the position of which the distance should be calculated.
 */
Mdouble FluxBoundary::getDistance(const Vec3D& position) const
{
    return distance_ - Vec3D::dot(position, normal_);
}

void FluxBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
    for (auto p = pH.begin(); p != pH.end(); ++p)
        checkBoundaryAfterParticleMoved(*p, pH);
}


/*!
 * \details Checks if particle has passed the boundary, and if so, counts it.
 * \param[in] p pointer to the particle which is to be checked
 * \param[out] pH the particle's ParticleHandler
 * \return FALSE, since the particle is not deleted by this boundary
 */
bool FluxBoundary::checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH)
{
    if (getDistance(p->getPosition() - p->getDisplacement()) > 0 && getDistance(p->getPosition()) < 0)
    {
        numberOfParticlesCrossedForw_++;
        massCrossedForw_ += p->getMass();
        volumeCrossedForw_ += p->getVolume();
    }
    else if (getDistance(p->getPosition() - p->getDisplacement()) < 0 && getDistance(p->getPosition()) > 0)
    {
        numberOfParticlesCrossedBack_++;
        massCrossedBack_ += p->getMass();
        volumeCrossedBack_ += p->getVolume();
    }
    
    return false;
}

/*!
 * \details Returns the number of particles that have crossed in either
 * direction. All of these quantities have to be divided by two, because
 * DPMBase.cc calls DPMBase::checkInteractionWithBoundaries() TWICE in each time
 * step. (The intention is that the first round resolves InsertionBoundaries and
 * the second resolves DeletionBoundaries, but this means FluxBoundaries are
 * updated twice.) */
unsigned int FluxBoundary::getNumberOfParticlesCrossedForw() const
{
    return numberOfParticlesCrossedForw_ / 2;
}

double FluxBoundary::getMassOfParticlesCrossedForw() const
{
    return massCrossedForw_ * 0.5;
}

double FluxBoundary::getVolumeOfParticlesCrossedForw() const
{
    return volumeCrossedForw_ * 0.5;
}

unsigned int FluxBoundary::getNumberOfParticlesCrossedBack() const
{
    return numberOfParticlesCrossedBack_ / 2;
}

double FluxBoundary::getMassOfParticlesCrossedBack() const
{
    return massCrossedBack_ * 0.5;
}

double FluxBoundary::getVolumeOfParticlesCrossedBack() const
{
    return volumeCrossedBack_ * 0.5;
}

unsigned int FluxBoundary::getNumberOfParticlesCrossedNet() const
{
    return (numberOfParticlesCrossedForw_ - numberOfParticlesCrossedBack_) / 2;
}

double FluxBoundary::getMassOfParticlesCrossedNet() const
{
    return (massCrossedForw_ - massCrossedBack_) * 0.5;
}

double FluxBoundary::getVolumeOfParticlesCrossedNet() const
{
    return (volumeCrossedForw_ - volumeCrossedBack_) * 0.5;
}

/*!
 * \details Reads a number of boundary properties from the given std::istream.
 * \param[in,out] is   the istream
 */
void FluxBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> normal_
       >> dummy >> scaleFactor_
       >> dummy >> distance_;
}

/*!
 * \details the deprecated version of the read-method. Should not be used by new 
 * users!
 * \deprecated Should be gone by Mercury 2.0. Use FluxBoundary::read() instead.
 */
void FluxBoundary::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> normal_ >> dummy >> scaleFactor_ >> dummy >> distance_;
}

/*!
 * \details Writes the boundary properties to an std::ostream. 
 * \param[out] os   the ostream the properties are to be written to.
 */
void FluxBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " normal " << normal_
       << " scaleFactor " << scaleFactor_
       << " distance " << distance_;
}

/*!
 * \details Returns the object's class name (i.e. 'FluxBoundary').
 * \return the object's class name
 */
std::string FluxBoundary::getName() const
{
    return "FluxBoundary";
}

