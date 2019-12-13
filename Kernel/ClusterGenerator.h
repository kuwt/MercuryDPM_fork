//
// Created by paolo on 15/05/19.
//

#ifndef ClusterGenerator_H
#define ClusterGenerator_H

#endif //ClusterGenerator_H

#include "BaseCluster.h"

/*!
 * \class ClusterGenerator
 * \brief This class allows the user to create clusters of particles. All particles will be of LinearPlasticViscoelasticSpecies
 *          and will have a final overlap defined by the user.
 */
class ClusterGenerator
{
public:

    /*!
     * \brief This constructor initializes all variables to default values. After initialization all values are passed
     *          to ClusterDPM with setters.
     */
    ClusterGenerator();

    /*!
     * \brief Default destructor
     */
    ~ClusterGenerator();

    /*!
     *  \brief runs the simulation that creates the cluster.
     */
    void create();


    BaseCluster clusterProperties;

    ParticleHandler fakeParticleHandler;

    //InteractionHandler* interactionHandler;



private:

    /*
     * ----------------------------------------------------
     *                    VARIABLES
     * ----------------------------------------------------
     */

    // Position
    Vec3D position_;

    // TIME
    //!\brief Ratio between collision time and time step: should be at least 50.
    Mdouble collisionTimeOverTimeStep_;

    //Particles
    //!\brief Size dispersity of particles: must be between 0 and 1
    Mdouble sizeDispersityParticle_;
    //!\brief Total number of particles.
    int nParticles_;

    //Cluster
    //!\brief Total number of particles.
    unsigned int idCluster_;
    //!\brief safety factor for the initial size of the cluster: this must be greater than 1.
    Mdouble clusterSizeSafetyFactor_;

    // Central force
    //!\brief Value of damping modulus for velocity.
    Mdouble velocityDampingModulus_;

    // Data analysis
    //! \brief Number of points used for creating internal structure's grid.
    int internalStructureGridLength_;

    // Energy
    //! \brief Energy ratio threshold under wich the simulation can be considered static.
    Mdouble energyRatioTolerance_;

    // File
    //!\brief bool used to define whether or not cluster data output must be created.
    bool isCdatOutputOn_;
    //!\brief bool used to define whether or not overlap data output must be created.
    bool isOverlOutputOn_;
    //!\brief bool used to define whether or not adjacency matrix output must be created.
    bool isAmatOutputOn_;
    //!\brief bool used to define whether or not cluster internal structure output must be created.
    bool isIntStrucOutputOn_;
    //!\brief bool used to define whether or not vtk output must be created.
    bool isVtkOutputOn_;
    //!\brief bool used to define whether or not data output must be created.
    bool isDataOutputOn_;
    //!\brief bool used to define whether or not restart output must be created.
    bool isRestartOutputOn_;
    //!\brief bool used to define whether or not fStat output must be created.
    bool isFStatOutputOn_;
    //!\brief bool used to define whether or not eneOutput output must be created.
    bool isEneOutputOn_;

    // Restart self test
    //!\brief bool used to define whether or not the restart self test is being computed.
    bool restartSelfTest_;



};
