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

#ifndef MERCURY_RANDOMCLUSTERINSERTIONBOUNDARY_H
#define MERCURY_RANDOMCLUSTERINSERTIONBOUNDARY_H
#include "BaseClusterInsertionBoundary.h"

/*!
 * RandomClusterInsertionBoundary
 * This class works exactly like CubeInsertionBoundary class but inserts clusters instead of particles.
 * MPI: It does not work with MPI. One of the reasons is the need to adapt function DPMBase::importParticlesAs().
 *      In particular particles and interaction should be imported as MPIParticle and MPIInteraction.
 *      In addition to this, the same should be done inside BaseCluster::actionsAfterSolve() when particles are centred
 *      around the desired position.
 */

class RandomClusterInsertionBoundary : public BaseClusterInsertionBoundary {

public:

    /*!
     * \brief Constructor: inherits from BaseClusterInsertionBoundary constructor.
     */
    RandomClusterInsertionBoundary();

    /*!
     * \brief Copy constructor with deep copy.
     */
    RandomClusterInsertionBoundary(const RandomClusterInsertionBoundary& other);

    /*!
     * \brief Destructor: default destructor.
     */
    ~RandomClusterInsertionBoundary() override;

    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    RandomClusterInsertionBoundary* copy() const override;

    /*!
     * \brief Sets the properties of the ClusterInsertionBoundary
     */
    void set(BaseParticle *particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
             Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax, Mdouble rMicroParticle);

    void set(BaseParticle &particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
             Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax, Mdouble rMicroParticle);

    /*!
     * \brief Sets the properties of the ClusterInsertionBoundary
     */
    void set(BaseParticle *particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
            unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax);

    void set(BaseParticle &particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
            unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax);

    //!\brief sets the number of particles per cluster
    void setNumberOfParticlesPerCluster(unsigned int nParticlesPeCluster);

    //!\brief inserts cluster, if no interactions are detected.
    void checkBoundaryBeforeTimeStep(DPMBase* md) override;

    //!\brief sets random position and velocity for the cluster.
    void placeParticle(BaseParticle* p, RNG& random) override;

private:
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;


};
#endif //MERCURY_RANDOMCLUSTERINSERTIONBOUNDARY_H
