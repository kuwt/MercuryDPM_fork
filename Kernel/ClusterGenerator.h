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
