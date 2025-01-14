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

#ifndef LEESEDWARDSBOUNDARY_H
#define LEESEDWARDSBOUNDARY_H

#include "BaseBoundary.h"

#include <functional>

class ParticleHandler;

/*!
 * \brief Class which creates a boundary with Lees-Edwards type periodic boundary conditions. 
 * \details A LeesEdwardsBoundary is like a PeriodicBoundary, but when a
 * particle crosses one edge and is copied to the other side then the particle
 * is also shifted. This sort of boundary is useful for studying shear flows. 
 * \details See also Lees and Edwards (J. Phys. C 1921, 
 * <a href="http://dx.doi.org/1088/0022-3719/5/15/006">doi:1088/0022-3719/5/15/006</a>). 
 * Inherits from BaseBoundary.
 * \todo Add link to paper by Lees-Edwards in the documentation of this class.
 * \todo Is implemented for 2D only now. Needs extension to 3D.
 */
class LeesEdwardsBoundary : public BaseBoundary
{
public:
    
    LeesEdwardsBoundary()
    {
#ifdef MERCURYDPM_USE_MPI
        MPIContainer& communicator = MPIContainer::Instance();
        if (communicator.getNumberOfProcessors() > 1)
        {
            logger(WARN,"LeesEdwardsBoundaries are currently not implemented in parallel MercuryDPM");
        }
#endif
    }

    /*!
     * \brief Copy constructor 
     */
    LeesEdwardsBoundary(const LeesEdwardsBoundary& other);
    
    /*!
     * \brief Sets all boundary properties
     */
    void
    set(std::function<Mdouble(Mdouble)> shift, std::function<Mdouble(Mdouble)> velocity, Mdouble left, Mdouble right,
        Mdouble down, Mdouble up);
    
    void updateBoundaries(Mdouble left, Mdouble right, Mdouble down, Mdouble up);
    
    /*!
     * \brief Reads all boundary properties from a stream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Writes all boundary properties to a stream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    /*!
     * \brief Creates a copy of the object
     */
    LeesEdwardsBoundary* copy() const override;
    
    /*!
     * \brief Returns distance from given particle to the closest horizontal wall 
     */
    Mdouble getHorizontalDistance(BaseParticle& p, bool& positive);
    
    /*!
     * \brief Returns distance from given particle to the closest vertical wall
     */
    Mdouble getVerticalDistance(BaseParticle& p, bool& positive);
    
    /*!
     * \brief Applies a horizontal shift to the given particle
     */
    void shiftHorizontalPosition(BaseParticle* p, bool positive);
    
    /*!
     * \brief Applies a vertical shift to the given particle
     */
    void shiftVerticalPosition(BaseParticle* p, bool positive);
    
    /*!
     * \brief Checks if particle crossed a boundary wall and if so, applies periodic shift
     */
    void checkBoundaryAfterParticleMoved(BaseParticle* p);
    
    /*!
     * \brief Checks if particles need to be adjusted after their position has been updated
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    void createPeriodicParticle(BaseParticle* p, ParticleHandler& pH) override;
    
    /*!
     * \brief Creates horizontal and vertical periodic copies of given particle, if needed
     */
    void createPeriodicParticles(ParticleHandler& pH) override;
    
    /*!
     * \brief Creates horizontal periodic copies of given particle, if needed
     */
    void createHorizontalPeriodicParticle(BaseParticle* p, ParticleHandler& pH);
    
    /*!
     * \brief Creates vertical periodic copies of given particle, if needed
     */
    void createVerticalPeriodicParticle(BaseParticle* p, ParticleHandler& pH);
    
    Mdouble getCurrentShift();
    
    Mdouble getCurrentVelocity();
    
    void setShift(std::function<Mdouble(Mdouble)>);
    
    void setVelocity(std::function<Mdouble(Mdouble)>);

private:
    Mdouble left_;      ///(signed) Horizontal distance between the left wall and the origin
    Mdouble right_;     ///(signed) Horizontal distance between the right wall and the origin
    Mdouble down_;      ///(signed) Vertical distance between the bottom wall and the origin
    Mdouble up_;        ///(signed) Vertical distance between the top wall and the origin
    std::function<Mdouble(Mdouble)> shift_;
    std::function<Mdouble(Mdouble)> velocity_;  ///Velocity difference between the top and bottom wall
};

#endif
