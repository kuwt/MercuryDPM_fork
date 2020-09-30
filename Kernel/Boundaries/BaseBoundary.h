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

#ifndef BaseBoundary_H
#define BaseBoundary_H

#include <set>
#include "BaseObject.h"
#include "Math/ExtendedMath.h"

class BoundaryHandler;

class PeriodicBoundaryHandler;

class ParticleHandler;

class BaseParticle;

class DPMBase;

/*!
 * \class BaseBoundary
 * \brief 
 * \details Inherits from BaseObject
 */
class BaseBoundary : public BaseObject
{
public:
    /*!
     * \brief default constructor.
     */
    BaseBoundary();
    
    /*!
     * \brief copy constructor
     */
    BaseBoundary(const BaseBoundary& b);
    
    /*!
     * \brief destructor
     */
    ~BaseBoundary() override;
    
    /*!
     * \brief Used to create a copy of the object
     * NB: purely virtual function
     */
    virtual BaseBoundary* copy() const = 0;
    
    /*!
     * \brief Reads the object's id_ from given istream
     * NB: purely virtual function, overriding the version of BaseObject
     */
    void read(std::istream& is) override = 0;
    
    /*!
     * \brief Adds object's id_ to given ostream
     * NB: purely virtual function, overriding the version of BaseObject
     */
    void write(std::ostream& os) const override = 0;
    
    
    /*!
     * \brief Creates a periodic particle in case of periodic boundaries in serial build
     */
    virtual void createPeriodicParticle(BaseParticle* p UNUSED, ParticleHandler& pH UNUSED);
    
    /*!
     * \brief Creates periodic copies of given particle in case of periodic boundaries
     */
    virtual void createPeriodicParticles(ParticleHandler& pH UNUSED);
    
    /*!
     * \brief Virtual function that does things to particles, each time step after particles have moved.
     */
    virtual void checkBoundaryAfterParticlesMove(ParticleHandler& pH);
    
    /*!
     * \brief Virtual function that does things before each time step.
     */
    virtual void checkBoundaryBeforeTimeStep(DPMBase* md);
    
    /*!
     * \brief Virtual function that does something after DPMBase::setupInitialConditions but before the first time step.
     */
    virtual void actionsBeforeTimeLoop();
    
    virtual void modifyGhostAfterCreation(BaseParticle* particle, int i);

    virtual void writeVTK(std::fstream& file) {}

    /*!
     * \brief Sets the boundary's BoundaryHandler
     */
    void setHandler(BoundaryHandler* handler);
    
    /*!
     * \brief Returns the boundary's BoundaryHandler
     */
    BoundaryHandler* getHandler() const;

#ifdef MERCURY_USE_MPI
    /*!
     * \brief Returns a list of particles that need to be deleted
     * \details Particles can't just be deleted in parallel code - as they need to
     * be flushed out of the mpi communication zones first, in a systematic order
     * \return a vector of base particles that need to be deleted
     */
    std::set<BaseParticle*>& getParticlesToBeDeleted();
#endif

protected:
#ifdef MERCURY_USE_MPI
    /*!
     * \brief Storage for particles in mpi communication zone to be particlesToBeDeleted
     * \details Because the communication lists must remain synced, we must make sure
     * the particles are removed in a systemaic order. So first the particles that need
     * to be deleted are moved into this list. After all particles are checked these particles will
     * be flushed from the mpi communication lista and afterwards destroyed.
     */
    std::set<BaseParticle*> particlesToBeDeleted_;

#endif

private:
    
    /*!
     * \brief pointer to the boundary's BoundaryHandler
     */
    BoundaryHandler* handler_;
    
};

#endif
