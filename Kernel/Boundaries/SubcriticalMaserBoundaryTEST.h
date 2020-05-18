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


#ifndef SUBCRITICALMASERBOUNDARYTEST_H
#define SUBCRITICALMASERBOUNDARYTEST_H

#include "Particles/BaseParticle.h"
#include "Boundaries/PeriodicBoundary.h"

class ParticleSpecies;

class SubcriticalMaserBoundaryTEST : public PeriodicBoundary
{
public:
    /*!
     * \brief MaserBoundary constructor
     */
    SubcriticalMaserBoundaryTEST();
    
    /*!
     * \brief destructor
     */
    ~SubcriticalMaserBoundaryTEST() override;
    
    /*!
     * \brief Creates a copy of this maser on the heap.
     */
    SubcriticalMaserBoundaryTEST* copy() const override;
    
    /*!
     * \brief reads boundary properties from istream
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief writes boundary properties to ostream
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;
    
    void actionsBeforeTimeLoop() override;
    
    /*!
     * \brief Creates periodic particles when the particle is a maser particle
     * and is sufficiently close to one of the boundary walls.
     */
    void createPeriodicParticle(BaseParticle* p, ParticleHandler& pH) override;
    
    /*!
     * \brief Shifts the particle to its 'periodic' position if it is a maser particle
     * and has crossed either of the walls. Creates a 'normal' particle at its current
     * position if it is a maser particle which crossed the RIGHT boundary wall.
     */
    bool checkBoundaryAfterParticleMoved(BaseParticle* p, ParticleHandler& pH) const;
    
    /*!
     * \brief Evaluates what the particles have to do after they have changed position
     */
    void checkBoundaryAfterParticlesMove(ParticleHandler& pH) override;
    
    /*!
     * \brief Activates the maser functionaly of this periodic boundary.
     */
    void activateMaser();

    /*!
     * \brief Stops copying particles, and act merely as a periodic domain
     */
    void deactivateMaser();

    /*!
     * \brief Returns whether the maser is activated or not.
     */
    bool isActivated() const;
    
    /*!
     * \brief sets the activate time of the maser
     */
    void setActivationTime(Mdouble time);
    
    /*!
     * \brief gets the distance to the closest wall if maser is inactive, otherwise distance to right wall
     */
    Mdouble getDistance(const Vec3D& position) const override;
    
    /*!
     * \brief returns the distance to the right wall
     */
    Mdouble getDistanceFromRight(const Vec3D& position) const;
    
    /*!
     * \brief modifies the periodic complexity to support a maser boundary
     */
    void modifyPeriodicComplexity(std::vector<int>& complexity, int& totalPeriodicComplexity,
                                  BaseParticle* particle, int i) const override;
    
    void modifyGhostAfterCreation(BaseParticle* particle, int i) override;
    
    
    /*!
     * \brief Checks before adding particles if the maser needs to be activated
     */
    void performActionsBeforeAddingParticles() override;
    
    void extendBottom() const;
    
    void copyExtraParticles() const;
    
    void setCopyFlowParticles(bool copyFlowParticles);

private:
    /*!
     * \brief Flag whether or not the gap is created and particles transformed already.
     */
    bool maserIsActivated_;
    
    /*!
     * \brief Time at which the maser opens
     */
    Mdouble activationTime_;
    
    /*!
     * \brief Flag for whether or not we copy a few blocks of flow particles in the front when activating the maser
     */
    bool copyFlowParticles_;
    
};

#endif
