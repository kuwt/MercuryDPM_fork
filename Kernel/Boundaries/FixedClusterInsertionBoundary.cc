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

#include "FixedClusterInsertionBoundary.h"

/*!
 * \details Default constructor: inherits from BaseClusterInsertionBoundary constructor.
 */
FixedClusterInsertionBoundary::FixedClusterInsertionBoundary() : BaseClusterInsertionBoundary()
{
    logger(DEBUG, "RandomClusterInsertionBoundary::RandomClusterInsertionBoundary() finished");
}

/*!
 * \details Copy constructor
 */
FixedClusterInsertionBoundary::FixedClusterInsertionBoundary(const FixedClusterInsertionBoundary& other)
        : BaseClusterInsertionBoundary(other) {

    logger(DEBUG, "RandomClusterInsertionBoundary::RandomClusterInsertionBoundary() finished");
}

/*!
 * \details Default Destructor.
 */
FixedClusterInsertionBoundary::~FixedClusterInsertionBoundary()
= default;

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 * \return      pointer to the copy on the heap
 */
FixedClusterInsertionBoundary* FixedClusterInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "BaseClusterInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new FixedClusterInsertionBoundary(*this);
}

/*!
 * \details Sets all the properties of the cuboidal insertion boundary.
 * \param[in] particleToCopy       Pointer to the BaseParticle which is used as a basis
 *                                 for clusters to be inserted
 * \param[in] positions            Vector containing all clusters positions
 * \param[in] radii                Vector containing all clusters radii
 * \param[in] velMin               Minimum velocity of inserted particles
 * \param[in] velMax               Maximum velocity of inserted particles
 * \param[in] rMicroParticle       Radius of the single particle composing the cluster.
 */
void FixedClusterInsertionBoundary::set(BaseParticle *particleToCopy,
                                        std::vector<Vec3D> positions, std::vector<Mdouble> radii,
                                        Vec3D velMin, Vec3D velMax, Mdouble rMicroParticle)
{
    setParticleToCopy(particleToCopy);

    setPositionsAndRadii(positions, radii);
    setVelocityRange(velMin, velMax);

    setRadiusMicroParticle(rMicroParticle);
}

//!\details After a few checks,  this sets positions and radii of the desired clusters.
void FixedClusterInsertionBoundary::setPositionsAndRadii(std::vector<Vec3D> clusterPositions, std::vector<Mdouble> clusterRadii)
{
    if (clusterPositions.size() != clusterRadii.size())
        logger(ERROR, "clusterPositions and clusterRadii have different size."
                      " clusterPositions.size() = %, clusterRadii.size() = %.",
               clusterPositions.size(), clusterRadii.size());
    else if (clusterPositions.empty())
        logger(ERROR, "clusterPositions is empty."
                      " clusterPositions.size() = %.", clusterPositions.size() );
    else if (clusterRadii.empty())
        logger(ERROR, "clusterRadii is empty."
                      " clusterRadii.size() = %.", clusterRadii.size());
    else
    {
        clusterPositions_ = clusterPositions;
        clusterRadii_ = clusterRadii;
    }
}

//!\details Here the insertion process takes place. Differently from RandomClusterInsertionBoundary,
//!         this is basically just a for loop.
void FixedClusterInsertionBoundary::checkBoundaryBeforeTimeStep(DPMBase* md)
{
    logger(VERBOSE, "In FixedClusterInsertionBoundary::checkBoundaryBeforeTimeStep\n");

    if (!isActivated_)
        return;

    /* Each time step, the InsertionBoundary attempts to fill up a region with
     * clusters.
     *
     * It first calls generateParticle() to get a randomised particle, subject
     * to a specified distribution over sizes and species. (The basic class
     * supports size dispersity only but PolydisperseInsertionBoundary will
     * support species dispersity.)
     * Then it repeatedly calls placeParticle(), which gives the cluster a
     * random location (and possibly velocity) in a specified,
     * geometry-dependent bound. Each time, it checks whether the new cluster
     * would have an interaction with another cluster or a wall.
     *
     * If it manages to do that within maxFailed_ tries, then:
     *   * the new cluster is created and inserted,
     *   * the failure counter is reset, and
     * the processes is repeated with a new generateParticle().
     *
     * Otherwise, the processes terminates for this time step.
     * */

    for (std::vector<Mdouble>::const_iterator i = clusterRadii_.begin(); i != clusterRadii_.end(); ++i) {


        auto p0 = generateParticle(md->random);

        placeParticle(p0, md->random);


        BaseCluster cluster;
        if(getRandomised())
            cluster.random.randomise();
        cluster.setCollisionTimeOverTimeStep(collisionTimeOverTimeStep_);
        cluster.setVelocityDampingModulus(velocityDampingModulus_);
        cluster.setNumberOfInternalStructurePoints(nInternalStructurePoints_);
        cluster.setEnergyRatioTolerance(energyRatioTolerance_);
        cluster.setSizeDispersityParticle(sizeDispersityParticle_);
        cluster.doCdatOutput(isCdatOutputOn_);
        cluster.doOverlOutput(isOverlOutputOn_);
        cluster.doAmatOutput(isAmatOutputOn_);
        cluster.doIntStrucOutput(isIntStrucOutputOn_);
        cluster.doVtkOutput(isVtkOutputOn_);
        cluster.doRestartOutput(isRestartOutputOn_);
        cluster.doFStatOutput(isFStatOutputOn_);
        cluster.doEneOutput(isEneOutputOn_);
        cluster.setClusterId(md->particleHandler.getNextGroupId());
        //Here, differently from RandomClusterInsertionBoundary the only possibility is to set radiusCluster and
        // radiusParticle.
        if (setRadiusParticleAndNotNumberOfParticles_)
            cluster.setRadiusParticle(radiusParticle_);
        cluster.setRadiusCluster(p0->getRadius());
        cluster.setPosition(p0->getPosition());
        cluster.setVelocity(p0->getVelocity());
        cluster.setParticleSpecies(
                dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(
                        md->speciesHandler.getObject(p0->getIndSpecies())));
        cluster.solve();

        md->importParticlesAs( cluster.particleHandler, cluster.interactionHandler, p0->getSpecies() );

        // Number of cluster inserted
        ++nClusterInserted_;
        // Total number of particles inserted
        numberOfParticlesInserted_ += cluster.particleHandler.getSize();
        // This is the total mass composed by every single particle (not cluster!) inserted
        massInserted_ += cluster.particleHandler.getMass();
        // This is the total volume composed by every single particle (not cluster!) inserted
        volumeInserted_ += cluster.particleHandler.getVolume();
        logger(VERBOSE, "Successfully inserted cluster %/%.", nClusterInserted_, clusterRadii_.size());

        delete p0;
    }

}

//!\brief Places particles according to vector clusterPositions_ and sets a random velocity, if required.
void FixedClusterInsertionBoundary::placeParticle(BaseParticle* p, RNG& random)
{
    if(getRandomised())
        random.randomise();
    Vec3D pos, vel;
    pos = clusterPositions_[nClusterInserted_];
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);
    p->setPosition(pos);
    p->setVelocity(vel);
}

//!\brief Sets cluster radii according to vector clusterRadii_.
BaseParticle* FixedClusterInsertionBoundary::generateParticle(RNG& random)
{
    BaseParticle* P = getParticleToCopy()->copy();
    P->setRadius(clusterRadii_[nClusterInserted_]);
    return P;
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'ClusterInsertionBoundary'
 */
std::string FixedClusterInsertionBoundary::getName() const
{
    return "FixedClusterInsertionBoundary";
}
