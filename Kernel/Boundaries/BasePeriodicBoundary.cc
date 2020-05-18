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

#include "BasePeriodicBoundary.h"
#include "PeriodicBoundaryHandler.h"
#include "ParticleHandler.h"
#include "MpiDataClass.h"
#include "DPMBase.h"

/*!
 * \details Default constructor
 */
BasePeriodicBoundary::BasePeriodicBoundary()
        : BaseBoundary()
{
    periodicHandler_ = nullptr;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"BasePeriodicBoundary::BasePeriodicBoundary() finished"<<std::endl;
#endif
}

/*!
 * \details Copy constructor
 */
///Note: shallow copy! Otherwise the HGrid causes a stack overflow.
BasePeriodicBoundary::BasePeriodicBoundary(const BasePeriodicBoundary& b)
        : BaseBoundary(b)
{
    periodicHandler_ = b.periodicHandler_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"BasePeriodicBoundary::BasePeriodicBoundary(const BasePeriodicBoundary&b) finished"<<std::endl;
#endif
}

/*!
 * \details Destructor
 */
BasePeriodicBoundary::~BasePeriodicBoundary()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout << "BasePeriodicBoundary::~BasePeriodicBoundary() finished"<<std::endl;
#endif
}

/*!
 * \details Reads the object's id_ from the given istream
 * \param[in,out] is    istream the id_ is read from
 */
void BasePeriodicBoundary::read(std::istream& is)
{
    BaseBoundary::read(is);
}

/*!
 * \details Adds the object's id_ to the given ostream
 * \param[in] os    ostream the id_ is added to
 */
void BasePeriodicBoundary::write(std::ostream& os) const
{
    BaseBoundary::write(os);
}


/*!
 * \details Sets the pointer to the BoundaryHandler the boundary belongs to
 * \param[in] handler   pointer to the boundary handler
 */
void BasePeriodicBoundary::setPeriodicHandler(PeriodicBoundaryHandler* periodicHandler)
{
    periodicHandler_ = periodicHandler;
}

/*!
 * \details Returns the pointer to the BoundaryHandler the boundary belongs to
 * \return pointer to the handler
 */
PeriodicBoundaryHandler* BasePeriodicBoundary::getPeriodicHandler() const
{
    return periodicHandler_;
}

/*!
 * \details checks whether particles have passed the boundary, and if so, does something
 * special with it. (i.e. deletion boundary, insertion boundary)
 * NB: virtual function
 * \param[out] pH  the particle handler.
 * \return returns TRUE if given particle actually did pass the boundary
 */
void BasePeriodicBoundary::checkBoundaryAfterParticlesMove(ParticleHandler& pH)
{
}

/*!
 * \details Checks the distance of given particle to the closest of both periodic 
 * walls, and creates a periodic copy of the particle if needed (i.e. if the particle
 * is closer to the periodic wall than the radius of the largest particle in the
 * system).
 * NOTE: This is only for a serial build - periodic particles work different in parallel
 * \param[in,out] pH    System's ParticleHandler, (1) from which the interaction radius
 *                      of its largest particle is retrieved to determine the maximum 
 *                      distance from the wall at which a particle should still have
 *                      a periodic copy created, and (2) to which a possible periodic
 *                      copy of the particle will be added
 */
void BasePeriodicBoundary::createPeriodicParticles(ParticleHandler& pH)
{
}

//TODO documentation
void BasePeriodicBoundary::modifyPeriodicComplexity(std::vector<int>& complexity, int& totalPeriodicComplexity,
                                                    BaseParticle* particle UNUSED, int i UNUSED) const
{
}

//TODO documentation
void BasePeriodicBoundary::performActionsBeforeAddingParticles()
{
}

