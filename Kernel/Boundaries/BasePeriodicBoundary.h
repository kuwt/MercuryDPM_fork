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

#ifndef BasePeriodicBoundary_H
#define BasePeriodicBoundary_H

#include "BaseBoundary.h"
#include "Math/ExtendedMath.h"
#include "Side.h"

class PeriodicBoundaryHandler;

/*!
 * \class BasePeriodicBoundary
 * \brief 
 * \details Inherits from BaseObject
 */
class BasePeriodicBoundary : public BaseBoundary
{
public:
    /*!
     * \brief default constructor.
     */
    BasePeriodicBoundary();
    
    /*!
     * \brief copy constructor
     */
    BasePeriodicBoundary(const BasePeriodicBoundary& b);
    
    /*!
     * \brief destructor
     */
    ~BasePeriodicBoundary() override;
    
    /*!
     * \brief Reads the object's id_ from given istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Adds object's id_ to given ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Used to create a copy of the object
     * NB: purely virtual function
     */
    BasePeriodicBoundary* copy() const override = 0;
    
    /*!
     * \brief Sets the periodicBoundaryHandler, required for parallel periodic boundaries
     */
    void setPeriodicHandler(PeriodicBoundaryHandler* periodicHandler);
    
    /*!
     * \brief Returns the periodic boundary handler
     */
    PeriodicBoundaryHandler* getPeriodicHandler() const;
    
    /*!
     * \brief Returns the distance between a particle and the closest boundary, required for any periodic boundary
     * \param[in] particle The distance between this particle and the closest boundary is calculated
     */
    virtual Mdouble getDistance(const BaseParticle& particle) const = 0;
    
    /*!
     * \brief Returns the distance between a position and the closest boundary
     * \param[in] position The distance between this particle and the closest boundary is calculated
     */
    virtual Mdouble getDistance(const Vec3D& position) const = 0;
    
    /*!
     * \brief Returns true if it is closer to the left boundary than the right boundary
     * \details Computes if a certain position is close to the left boundary (true)
     * or if it is not close to the left boundary (false)
     * \param[in] position The position which is being checked
    */
    virtual bool isClosestToLeftBoundary(const Vec3D& position) const = 0;
    
    /*!
     * \brief Shifts the position (and velocity) of to the ghost particle
    * \details Shifts the position of a particle to the other boundary.
    * Note: In some cases it doesnt only shift the position, but also other quantities such
    * as velocity
    * \param[in] particle Pointer to the particle that will shift position
     */
    virtual void shiftPosition(BaseParticle* particle) const = 0;
    
    /*!
     * \brief Creates periodic ocpies of given particle in case of periodic boundaries in serial build
     */
    void createPeriodicParticles(ParticleHandler& pH) override;
    
    /*!
     * \brief Virtual function that does things to particles, each time step after particles have moved.
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!
     * \brief Modifies periodic complexity of a particle if necessary (i.e. maser boundary)
     */
    virtual void modifyPeriodicComplexity(std::vector<int>& complexity, int& totalPeriodicComplexity,
                                          BaseParticle* particle, int i) const;
    
    /*!
     * \brief Actions that need to be performed before adding new ghost particles
     */
    virtual void performActionsBeforeAddingParticles();

private:
    /*!
     * \brief pointer to the periodic boundary handler
     */
    PeriodicBoundaryHandler* periodicHandler_;
};

#endif
