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

#include "ClusterGenerator.h"


/*!
 * \details This constructor initializes all variables to default values. After initialization all values are passed
 *          to ClusterDPM with setters.
 *          IMPORTANT: radiusParticle here is not initialized nor passed to ClusterDPM in order to
 *          force the user to set it, and same is for the species.
 */

ClusterGenerator::ClusterGenerator()
{
    /*
     * ----------------------------------------------------
     *      Initializing variables to default values
     * ----------------------------------------------------
     */

    //Position
    position_={0,0,0};

    // Time
    collisionTimeOverTimeStep_ = 51;

    // Particles
    sizeDispersityParticle_ = 0.0;
    nParticles_ = 100;

    // Cluster
    clusterSizeSafetyFactor_ = 2.0;
    idCluster_ = 0;

    // Central force
    velocityDampingModulus_ = 0.9;

    // Data analysis
    internalStructureGridLength_ = 300;

    // Energy
    energyRatioTolerance_ = 1.0e-8;

    // File
    isCdatOutputOn_ = true;//true;
    isOverlOutputOn_ = true;//true;
    isAmatOutputOn_ = false;//true;
    isIntStrucOutputOn_ = false; // Computational intensive
    isVtkOutputOn_ = false;
    isDataOutputOn_= true;
    isRestartOutputOn_ = false;//true;
    isFStatOutputOn_ = false;
    isEneOutputOn_ = false;

    // Restart self test
    restartSelfTest_ = false;


    /*
     * ----------------------------------------------------
     *           Passing values to BaseCluster
     * ----------------------------------------------------
     */

    clusterProperties.setCollisionTimeOverTimeStep(collisionTimeOverTimeStep_);

    clusterProperties.setSizeDispersityParticle(sizeDispersityParticle_);

    clusterProperties.setNumberOfParticles(nParticles_);

    clusterProperties.setClusterId(idCluster_);

    clusterProperties.setVelocityDampingModulus(velocityDampingModulus_);

    clusterProperties.setInternalStructureGridLength(internalStructureGridLength_);

    clusterProperties.setEnergyRatioTolerance(energyRatioTolerance_);

    clusterProperties.doCdatOutput(isCdatOutputOn_);

    clusterProperties.doOverlOutput(isOverlOutputOn_);

    clusterProperties.doAmatOutput(isAmatOutputOn_);

    clusterProperties.doIntStrucOutput(isIntStrucOutputOn_);

    clusterProperties.doVtkOutput(isVtkOutputOn_);

    clusterProperties.doRestartOutput(isRestartOutputOn_);

    clusterProperties.doFStatOutput(isFStatOutputOn_);

    clusterProperties.doEneOutput(isEneOutputOn_);

    logger(DEBUG, "ClusterGenerator::ClusterGenerator() finished");
}

/*!
 * \details Default destructor
 */
ClusterGenerator::~ClusterGenerator()
{
    logger(DEBUG, "ClusterGenerator::~ClusterGenerator() finished");
}

/*!
 *  \details runs the simulation that creates the cluster.
 */
void ClusterGenerator::create()
{

    clusterProperties.solve();
}
