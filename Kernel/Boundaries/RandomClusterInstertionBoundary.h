//
// Created by paolo on 23-10-19.
//

#ifndef MERCURY_RANDOMCLUSTERINSTERTIONBOUNDARY_H
#define MERCURY_RANDOMCLUSTERINSTERTIONBOUNDARY_H

#endif //MERCURY_RANDOMCLUSTERINSTERTIONBOUNDARY_H


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


};