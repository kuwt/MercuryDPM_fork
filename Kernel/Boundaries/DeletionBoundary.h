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

#ifndef DELETIONBOUNDARY_H
#define DELETIONBOUNDARY_H

#include "BaseBoundary.h"
#include "Math/Vector.h"

class ParticleHandler;

class BaseParticle;

/*!
 * \class DeletionBoundary
 * \brief Used for removing particles from the problem.
 * Inherits from BaseBoundary.
 * By default, a plane that deletes everything past it, but there are
 * derived classes such as CubeDeletionBoundary.
 */
class DeletionBoundary : public BaseBoundary
{
public:
    /*!
     * \brief default constructor
     */
    DeletionBoundary();

    /*!
     * \brief default copy constructor
     */
    DeletionBoundary(DeletionBoundary const &d)
        : normal_(d.normal_), scaleFactor_(d.scaleFactor_), distance_(d.distance_),
          numberOfParticlesDeleted_(d.numberOfParticlesDeleted_), massDeleted_(d.massDeleted_),
          volumeDeleted_(d.volumeDeleted_), isActivated_(d.isActivated_), trackOutflow_(d.trackOutflow_) {}

    /*!
     * \brief destructor
     */
    ~DeletionBoundary() override;
    
    /*!
     * \brief Copy method; creates copy on the heap and returns a pointer to it.
     */
    DeletionBoundary* copy() const override;
    
    /*!
     * \brief Sets boundary position based on a normal and distance.
     */
    virtual void set(const Vec3D& normal, Mdouble distance);
    
    /*!
     * \brief Sets the boundary's distance property to the given one.
     */
    void move(Mdouble position);
    
    /*!
     * \brief Returns a negative value if and only if the particle is inside the
     * boundary (and therefore to be deleted).
     */
    virtual Mdouble getDistance(const Vec3D& position) const;
    
    /*!
     * \brief Checks if particle is inside the boundary, and deletes the particle if so.
     */
    /// \todo: MX: update the above comment
    bool checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH);
    
    /*!
     * \todo MX: I changed the syntax to a more mpi compatible (and neater) code. Need to update the comments of the checkBoundaryAfterParticleMoved function
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!
     * \brief Gets the number of particles deleted by the boundary.
     */
    unsigned int getNumberOfParticlesDeleted() const;
    
    double getMassOfParticlesDeleted() const;
    
    double getVolumeOfParticlesDeleted() const;
    
    void reset();
    
    
    /*! 
     * \brief Turns on the DeletionBoundary.
     */
    void activate();
    
    /*!
     * \brief Turns off the DeletionBoundary.
     */
    void deactivate();

    /*!
     * \brief Turns on the outflow tracker.
     */
    void trackOutflow(bool trackOutflow = true) {trackOutflow_ = trackOutflow;}

    /*!
     * \brief Reads some boundary properties from an std::istream.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Deprecated read method. use DeletionBoundary::read() instead.
     */
    MERCURY_DEPRECATED
    void oldRead(std::istream& is);
    
    /*!
     * \brief Writes the boundary properties to an std::ostream.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;

private:
    /*!
     * \brief outward unit normal vector
     */
    Vec3D normal_;
    
    /*!
     * \brief This is the factor to rescale the given normal vector to a unit vectors. 
     * \details NB: Not only the normal vector is rescaled by this factor, also 
     * the 'distance' from the origin of the boundary is scaled by this factor! Also,
     * once the boundary position is set with DeletionBoundary::set(), the arguments 
     * of any reset of the distance_ property  (i.e. usage of DeletionBoundary::move()) 
     * will be rescaled by the same factor!
     */
    Mdouble scaleFactor_;
    
    /*!
     * \brief The boundary's distance from the origin.
     */
    Mdouble distance_;
    
    /*!
     * \brief Number of particles that have been deleted by this boundary.
     */
    unsigned int numberOfParticlesDeleted_;
    double massDeleted_;
    double trackMassDeleted_;
    double volumeDeleted_;
    
    /*!
     * \brief The DeletionBoundary is activated by default. If the
     * DeletionBoundary is deactivated, then it deletes no particles.
     */
    bool isActivated_;

    bool trackOutflow_;

    std::ofstream tracker;

};

#endif
