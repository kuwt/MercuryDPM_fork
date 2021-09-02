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

#ifndef MERCURY_FIXEDCLUSTERINSERTIONBOUNDARY_H
#define MERCURY_FIXEDCLUSTERINSERTIONBOUNDARY_H
#include "BaseClusterInsertionBoundary.h"

/*!
 * RandomClusterInsertionBoundary
 * This class works exactly like CubeInsertionBoundary class but inserts clusters instead of particles, and the insertion
 *      process is not random, but determined by radii and positions specified by the user.
 * MPI: It does not work with MPI. One of the reasons is the need to adapt function DPMBase::importParticlesAs().
 *      In particular particles and interaction should be imported as MPIParticle and MPIInteraction.
 *      In addition to this, the same should be done inside BaseCluster::actionsAfterSolve() when particles are centred
 *      around the desired position.
 */
class FixedClusterInsertionBoundary : public BaseClusterInsertionBoundary {

public:

    /*!
     * \brief Constructor: inherits from BaseClusterInsertionBoundary constructor.
     */
    FixedClusterInsertionBoundary();

    /*!
     * \brief Copy constructor with deep copy.
     */
    FixedClusterInsertionBoundary(const FixedClusterInsertionBoundary& other);

    /*!
     * \brief Destructor: default destructor.
     */
    ~FixedClusterInsertionBoundary() override;

    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    FixedClusterInsertionBoundary* copy() const override;

    /*!
     * \brief Sets the properties of the ClusterInsertionBoundary
     */
    void set(BaseParticle *particleToCopy,
             std::vector<Vec3D> positions, std::vector<Mdouble> radii,
             Vec3D velMin, Vec3D velMax, Mdouble rMicroParticle);

    //!\brief this sets positions and radii of the desired clusters.
    void setPositionsAndRadii(std::vector<Vec3D> clusterPositions, std::vector<Mdouble> clusterRadii);

    //!\brief inserts cluster: differently from RandomClusterInsertionBoundary, here no check for interaction is computed.
    void checkBoundaryBeforeTimeStep(DPMBase* md) final;

    //!\brief Places particles according to vector clusterPositions_ and sets a random velocity, if required.
    void placeParticle(BaseParticle* p, RNG& random) final;

    //!\brief Sets cluster radii according to vector clusterRadii_.
    BaseParticle* generateParticle(RNG& random) final;

private:
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;


};
#endif //MERCURY_FIXEDCLUSTERINSERTIONBOUNDARY_H
