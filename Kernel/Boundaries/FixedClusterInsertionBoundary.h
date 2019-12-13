//
// Created by paolo on 23-10-19.
//

#ifndef MERCURY_FIXEDCLUSTERINSERTIONBOUNDARY_H
#define MERCURY_FIXEDCLUSTERINSERTIONBOUNDARY_H

#endif //MERCURY_FIXEDCLUSTERINSERTIONBOUNDARY_H


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
    BaseParticle* generateParticle(RNG &random) final;

private:
    /*!
     * \brief Returns the name of the object
     */
    std::string getName() const override;


};