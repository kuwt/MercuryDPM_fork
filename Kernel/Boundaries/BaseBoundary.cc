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

#include "BaseBoundary.h"
#include "BoundaryHandler.h"

/*!
 * \details Default constructor
 */
BaseBoundary::BaseBoundary()
{
    handler_ = nullptr;
    logger(DEBUG, "BaseBoundary::BaseBoundary() finished");
}

/*!
 * \details Copy constructor
 */
///Note: shallow copy! Otherwise the HGrid causes a stack overflow.
BaseBoundary::BaseBoundary(const BaseBoundary& b)
        : BaseObject(b)
{
    handler_ = b.handler_;
    logger(DEBUG, "BaseBoundary::BaseBoundary(const BaseBoundary &b) finished");
}

/*!
 * \details Destructor
 */
BaseBoundary::~BaseBoundary()
{
    logger(DEBUG, "BaseBoundary::~BaseBoundary() finished");
}

/*!
 * \details Reads the object's id_ from the given istream
 * \param[in,out] is    istream the id_ is read from
 */
void BaseBoundary::read(std::istream& is)
{
    BaseObject::read(is);
}

/*!
 * \details Adds the object's id_ to the given ostream
 * \param[in] os    ostream the id_ is added to
 */
void BaseBoundary::write(std::ostream& os) const
{
    BaseObject::write(os);
}

/*!
 * \details Used to create a single periodic copy of a particle in classes which implement periodic
 *   boundary conditions
 * \param[in] p     particle of which periodic copies are to be created
 * \param[in] pH   the particle handler
 */
void BaseBoundary::createPeriodicParticle(BaseParticle* p UNUSED, ParticleHandler& pH UNUSED)
{
}

/*!
 * \details Used to create periodic copies of particles in classes which 
 * implement periodic boundary conditions
 * NB: virtual function
 * \param[in] pH   the particle handler
 */
void BaseBoundary::createPeriodicParticles(ParticleHandler& pH UNUSED)
{
}

/*!
 * \details What this does depends on the type of boundary. 
 * For example, an InsertionBoundary introduces new particles (and how it does
 * that in turn depends on the type of InsertionBoundary).
 * \param[in] md    the problem's DPMBase object
 */
void BaseBoundary::checkBoundaryBeforeTimeStep(DPMBase* md)
{
}

/*!
 * \details checks whether given particle passed the boundary, and if so, does something
 * special with it. This is called after the particles moved, but before the force-computation.
 * NB: virtual function
 * \param[in] P    Particle checked
 * \param[out] pH  the particle handler.
 * \return returns TRUE if given particle actually did pass the boundary
 */
void BaseBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
}

void BaseBoundary::modifyGhostAfterCreation(BaseParticle* particle, int i)
{
}


/*!
 * \details Can be used to perform actions before the time loop, but after setupInitialConditions.
 */
void BaseBoundary::actionsBeforeTimeLoop()
{
    logger(VERBOSE, "In BaseBoundary::checkBoundaryBeforeTimeLoop\n");
}

/*!
 * \details Sets the pointer to the BounadaryHandler the boundary belongs to
 * \param[in] handler   pointer to the boundary handler
 */
void BaseBoundary::setHandler(BoundaryHandler* handler)
{
    handler_ = handler;
}

/*!
 * \details Returns the pointer to the BoundaryHandler the boundary belongs to
 * \return pointer to the handler
 */
BoundaryHandler* BaseBoundary::getHandler() const
{
    return handler_;
}

#ifdef MERCURY_USE_MPI
/*!
 * \brief Returns a list of particles that need to be deleted.
 */
std::set<BaseParticle*>& BaseBoundary::getParticlesToBeDeleted()
{
    return particlesToBeDeleted_;
}
#endif

