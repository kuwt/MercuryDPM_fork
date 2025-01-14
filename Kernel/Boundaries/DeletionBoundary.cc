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

#include <limits>

#include "DeletionBoundary.h"
#include "ParticleHandler.h"
#include "Particles/BaseParticle.h"
#include <algorithm>
#include "Domain.h"
#include "BoundaryHandler.h"
#include "DPMBase.h"

/*!
 * \details Default constructor (calls the parent-constructor of BaseBoundary as well)
 */
DeletionBoundary::DeletionBoundary()
        : BaseBoundary()
{
    distance_ = std::numeric_limits<double>::quiet_NaN();
    scaleFactor_ = std::numeric_limits<double>::quiet_NaN();
    isActivated_ = true;
    trackOutflow_ = false;
    numberOfParticlesDeleted_ = 0;
    massDeleted_ = 0;
    volumeDeleted_ = 0;
    
    logger(DEBUG, "DeletionBoundary::DeletionBoundary() finished", true);
}

/*!
 * \details Destructor
 */
DeletionBoundary::~DeletionBoundary()
{
    tracker.close();
    logger(DEBUG, "DeletionBoundary::~DeletionBoundary() finished", true);
}

/*!
 * \details Copy function, which creates a copy and returns a pointer to that copy
 * (on the heap)
 * \return pointer to the copy
 */
DeletionBoundary* DeletionBoundary::copy() const
{
    return new DeletionBoundary(*this);
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
void DeletionBoundary::set(const Vec3D& normal, Mdouble distance)
{
    scaleFactor_ = 1. / std::sqrt(Vec3D::dot(normal, normal));
    normal_ = normal * scaleFactor_;
    distance_ = distance * scaleFactor_;
}

/*!
 * \details Resets the boundary's 'distance' from the origin to be the one given.
 * \param[in] distance the new 'distance' between boundary and origin.
 * see also comments of DeletionBoundary::set().
 */
void DeletionBoundary::move(Mdouble distance)
{
    distance_ = distance * scaleFactor_;
}

/*!
 * \details The default DeletionBoundary is a half-space (specified by a
 * plane) but the geometry can be overridden (e.g. by CubeDeletionBoundary).
 * Calculates the shortest distance between the wall and given position.
 * \param[in] position the position of which the distance should be calculated.
 */
Mdouble DeletionBoundary::getDistance(const Vec3D& position) const
{
    return distance_ - Vec3D::dot(position, normal_);
}

/*!
 * \details Checks if particle has passed the boundary, and if so, deletes the 
 * particle.
 * \param[in] p pointer to the particle which is to be checked
 * \param[out] pH the particle's ParticleHandler, from which
 * the particle is removed in case it has passed the boundary. 
 * \return TRUE if the particle has actually passed the boundary and is 
 * thus deleted.
 * \todo JMFT: The mass and volume counters need to be updated in the MPI code.
 */
bool DeletionBoundary::checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH)
{
    if (getDistance(p->getPosition()) < 0)
    {
        if (trackOutflow_) {
            auto dpm = getHandler()->getDPMBase();
            trackMassDeleted_ += p->getMass();
            if (!tracker.is_open()) {
                std::string name  = dpm->getName()
                        + helpers::toString(getIndex())
                        + ".out" + (NUMBER_OF_PROCESSORS==1?"":std::to_string(PROCESSOR_ID));
                logger(INFO, "Open file %", name);
                tracker.open(name);
                tracker << std::setw(13) << "Time "
                        << std::setw(8) << "Species "
                        << std::setw(8) << "ID "
                        << std::setw(13) << "Position_x "
                        << std::setw(13) << "Position_y "
                        << std::setw(13) << "Position_z "
                        << std::setw(13) << "Velocity_x "
                        << std::setw(13) << "Velocity_y "
                        << std::setw(13) << "Velocity_z "
                        << std::setw(13) << "Radius "
                        << std::setw(13) << "AnVelocity_x "
                        << std::setw(13) << "AnVelocity_y "
                        << std::setw(13) << "AnVelocity_z "
                        << std::setw(13) << "Mass "
                        << std::setw(13) << "TotalMass"
                        << '\n';
            }
            tracker << std::setw(12) << dpm->getTime() << ' '
                    << std::setw(8) << p->getIndSpecies() << ' '
                    << std::setw(8) << p->getId() << ' '
                    << std::setw(12) << p->getPosition().X << ' '
                    << std::setw(12) << p->getPosition().Y << ' '
                    << std::setw(12) << p->getPosition().Z << ' '
                    << std::setw(12) << p->getVelocity().X << ' '
                    << std::setw(12) << p->getVelocity().Y << ' '
                    << std::setw(12) << p->getVelocity().Z << ' '
                    << std::setw(12) << p->getRadius() << ' '
                    << std::setw(12) << p->getAngularVelocity().X << ' '
                    << std::setw(12) << p->getAngularVelocity().Y << ' '
                    << std::setw(12) << p->getAngularVelocity().Z << ' '
                    << std::setw(12) << p->getMass() << ' '
                    << std::setw(12) << trackMassDeleted_ << '\n';
        }
        else
            trackerCustom_(tracker, p);

        #ifdef MERCURYDPM_USE_MPI
            //Check if the particle is in the mpi communication zone
            if(p->isInMPIDomain())
            {
                //Add the particle to a list so we can flush it later on
                //The particles are not deleted yet. This is done in dpmBase after this function is called
                particlesToBeDeleted_.insert(p);
                return false;
            }
            else
            {
                pH.removeGhostObject(p->getIndex());
                return true;
            }
        #else
            numberOfParticlesDeleted_++;
            massDeleted_ += p->getMass();
            volumeDeleted_ += p->getVolume();
            pH.removeObject(p->getIndex());
            pH.computeLargestParticle();
            pH.computeSmallestParticle();
            return true;
        #endif
    }
    else
    {
        return false;
    }
}

void DeletionBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
    if (!isActivated_)
        return;

#ifdef MERCURYDPM_USE_MPI
    particlesToBeDeleted_.clear();
#endif
    for (unsigned int i = 0; i < pH.getSize(); i++)
    {
        //If the particle is deleted, change the iterator
        if (checkBoundaryAfterParticleMoved(pH.getObject(i), pH))
        {
            i--;
        }
    }
}

double DeletionBoundary::getMassOfParticlesDeleted() const
{
    return massDeleted_;
}

double DeletionBoundary::getVolumeOfParticlesDeleted() const
{
    return volumeDeleted_;
}

void DeletionBoundary::reset()
{
    numberOfParticlesDeleted_ = 0;
    massDeleted_ = 0;
    volumeDeleted_ = 0;
}

/*!
 * \details Reads a number of boundary properties from the given std::istream.
 * \param[in,out] is   the istream
 */
void DeletionBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
    std::string dummy;
    is >> dummy >> normal_
       >> dummy >> scaleFactor_
       >> dummy >> distance_
       >> dummy >> numberOfParticlesDeleted_
       >> dummy >> massDeleted_
       >> dummy >> volumeDeleted_
       >> dummy >> isActivated_
       >> dummy >> trackOutflow_;
}

/*!
 * \details the deprecated version of the read-method. Should not be used by new 
 * users!
 * \deprecated Should be gone by Mercury 2.0. Use DeletionBoundary::read() instead.
 */
void DeletionBoundary::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> normal_ >> dummy >> scaleFactor_ >> dummy >> distance_;
}

/*!
 * \details Writes the boundary properties to an std::ostream. 
 * \param[out] os   the ostream the properties are to be written to.
 */
void DeletionBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
    os << " normal " << normal_
       << " scaleFactor " << scaleFactor_
       << " distance " << distance_
       << " numberOfParticlesDeleted " << numberOfParticlesDeleted_
       << " massDeleted " << massDeleted_
       << " volumeDeleted " << volumeDeleted_
       << " isActivated " << isActivated_
       << " trackOutflow " << trackOutflow_;
}

/*!
 * \details Returns the object's class name (i.e. 'DeletionBoundary').
 * \return the object's class name
 */
std::string DeletionBoundary::getName() const
{
    return "DeletionBoundary";
}

/*!
 * \details Returns the number of particles deleted by this boundary.
 * TODO For some reason, this causes a segfault --- DO NOT USE.
 * valgrind complains that 'numberOfParticlesDeleted_' is not initialised.
 */
unsigned int DeletionBoundary::getNumberOfParticlesDeleted() const
{
    return numberOfParticlesDeleted_;
}

void DeletionBoundary::activate()
{
    isActivated_ = true;
}

void DeletionBoundary::deactivate()
{
    isActivated_ = false;
}
