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

#include "RandomClusterInsertionBoundary.h"






/*!
 * \details Default constructor: inherits from BaseClusterInsertionBoundary constructor.
 */
RandomClusterInsertionBoundary::RandomClusterInsertionBoundary() : BaseClusterInsertionBoundary()
{
    logger(DEBUG, "RandomClusterInsertionBoundary::RandomClusterInsertionBoundary() finished");
}

/*!
 * \details Copy constructor
 */
RandomClusterInsertionBoundary::RandomClusterInsertionBoundary(const RandomClusterInsertionBoundary& other)
        : BaseClusterInsertionBoundary(other) {

    logger(DEBUG, "RandomClusterInsertionBoundary::RandomClusterInsertionBoundary() finished");
}

/*!
 * \details Default Destructor.
 */
RandomClusterInsertionBoundary::~RandomClusterInsertionBoundary()
= default;

/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 * \return      pointer to the copy on the heap
 */
RandomClusterInsertionBoundary* RandomClusterInsertionBoundary::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "BaseClusterInsertionBoundary::copy() const finished" << std::endl;
#endif
    return new RandomClusterInsertionBoundary(*this);
}

/*!
 * \details Sets all the properties of the cuboidal insertion boundary.
 * \param[in] particleToCopy    Pointer to the BaseParticle which is used as a basis
 *                              for clusters to be inserted
 * \param[in] maxFailed         The maximum number of times the insertion of a
 *                              particle may be tried and failed before the insertion
 *                              of particles is considered done
 *                              NB: this property is used in the parent's
 *                              InsertionBoundary::checkBoundaryBeforeTimeStep().
 * \param[in] posMin            First defining corner of cuboidal insertion boundary
 * \param[in] posMax            Second defining corner of cuboidal insertion boundary
 * \param[in] velMin            Minimum velocity of inserted particles
 * \param[in] velMax            Maximum velocity of inserted particles
 * \param[in] radMin            Minimum radius of inserted particles
 * \param[in] radMax            Maximum radius of inserted particles
 * \param[in] rMicroParticle    Radius of the single particle composing the cluster.
 */
void RandomClusterInsertionBoundary::set(BaseParticle *particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                       Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax,
                                       Mdouble rMicroParticle)
{
    setParticleToCopy(particleToCopy);
    setRadiusRange(radMin, radMax);

    setMaxFailed(maxFailed);
    setGeometry(posMin, posMax, velMin, velMax);

    setRadiusMicroParticle(rMicroParticle);
}

void RandomClusterInsertionBoundary::set(BaseParticle &particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                       Vec3D velMin, Vec3D velMax, Mdouble radMin, Mdouble radMax,
                                       Mdouble rMicroParticle)
{
    set(&particleToCopy, maxFailed, posMin, posMax, velMin, velMax, radMin, radMax, rMicroParticle);
}

/*!
 * \details Sets all the properties of the cuboidal insertion boundary.
 * \param[in] particleToCopy            Pointer to the BaseParticle which is used as a basis
 *                                      for clusters to be inserted
 * \param[in] maxFailed                 The maximum number of times the insertion of a
 *                                      particle may be tried and failed before the insertion
 *                                      of particles is considered done
 *                                      NB: this property is used in the parent's
 *                                      InsertionBoundary::checkBoundaryBeforeTimeStep().
 * \param[in] posMin                    First defining corner of cuboidal insertion boundary
 * \param[in] posMax                    Second defining corner of cuboidal insertion boundary
 * \param[in] nParticlesPerCluster      Number of particles composing the cluster.
 * \param[in] velMin                    Minimum velocity of inserted particles
 * \param[in] velMax                    Maximum velocity of inserted particles
 * \param[in] radMin                    Minimum radius of inserted particles
 * \param[in] radMax                    Maximum radius of inserted particles
 *
 * Important: this function differs from the class above because gives the possiblity to set the number of particles
 *              instead of the radius of the micro particle.
 */

void RandomClusterInsertionBoundary::set(BaseParticle *particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                         unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax,
                                         Mdouble radMin, Mdouble radMax)
{
    setParticleToCopy(particleToCopy);
    setRadiusRange(radMin, radMax);

    setMaxFailed(maxFailed);
    setGeometry(posMin, posMax, velMin, velMax);

    setNumberOfParticlesPerCluster(nParticlesPerCluster);
}

void RandomClusterInsertionBoundary::set(BaseParticle &particleToCopy, unsigned int maxFailed, Vec3D posMin, Vec3D posMax,
                                         unsigned int nParticlesPerCluster, Vec3D velMin, Vec3D velMax,
                                         Mdouble radMin, Mdouble radMax)
{
    set(&particleToCopy, maxFailed, posMin, posMax, nParticlesPerCluster, velMin, velMax, radMin, radMax);
}

void RandomClusterInsertionBoundary::setNumberOfParticlesPerCluster(unsigned int nParticlesPeCluster)
{
    if (nParticlesPeCluster <= 0)
        logger(ERROR, "The number of particles for a single cluster must be greater than zero. nParticlesPeCluster = %", nParticlesPeCluster);
    else {
        nParticles_ = nParticlesPeCluster;
        setRadiusParticleAndNotNumberOfParticles_ = false;
    }
}

void RandomClusterInsertionBoundary::checkBoundaryBeforeTimeStep(DPMBase* md)
{
    logger(VERBOSE, "In RandomClusterInsertionBoundary::checkBoundaryBeforeTimeStep\n");

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

    // Keep count of how many successive times we have failed to place a new
    // particle.
    unsigned int failed = 0;
    while (failed <= maxFailed_ && insertParticle(md->getNextTime())) // 'generating' loop
    {
        /* Generate random *intrinsic* properties for the new particle. */
        logger(VERBOSE, "about to call generateParticle\n");

        auto p0 = generateParticle(md->random);
        logger(VERBOSE, "generated a particle with intrinsics %", p0);

        while (true) // 'placing' loop
        {
            /* Generate extrinsic properties (position and velocity) for this
             * new particle. */

            placeParticle(p0, md->random);

            logger(VERBOSE, "attempting to place particle at %, vel %", p0->getPosition(), p0->getVelocity());

#ifdef MERCURY_USE_MPI
            /* Communicate the new particle's properties by setHandler (note
             * that this doesn't actually add the particle to the handler). */
            if (NUMBER_OF_PROCESSORS > 1)
            {
                MPIParticle particle;
                // 	//Every domain generates a particle (to get the species right etc)

                //Send particle data from root to other processors to sync the particle properties
                if (PROCESSOR_ID == 0)
                {
                    particle = copyDataFromParticleToMPIParticle(p0);
                }

                MPIContainer::Instance().broadcast(&particle,MercuryMPIType::PARTICLE);

                //Process the received data
                if (PROCESSOR_ID != 0)
                {
                    copyDataFromMPIParticleToParticle(&particle, p0, &(md->particleHandler));
                }
            }
#endif
            p0->setHandler(&md->particleHandler);

            /* Check whether the particle has any interactions. */
            if (md->checkParticleForInteraction(*p0))
            {
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
                if (setRadiusParticleAndNotNumberOfParticles_)
                    cluster.setRadiusParticle(radiusParticle_);
                else
                    cluster.setNumberOfParticles(nParticles_);
                cluster.setRadiusCluster(p0->getRadius());
                cluster.setPosition(p0->getPosition());
                cluster.setVelocity(p0->getVelocity());
                cluster.setParticleSpecies(
                        dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(
                                md->speciesHandler.getObject(p0->getIndSpecies())));
#ifdef MERCURY_USE_MPI
                cluster.setNumberOfDomains(md->getNumberOfDomains());
                cluster.setDomain(md->getMin(), md->getMax());
                cluster.speciesHandler.copyAndAddObject(md->speciesHandler.getObject(0));
#endif
                cluster.solve();

                md->importParticlesAs( cluster.particleHandler, cluster.interactionHandler, p0->getSpecies() );

                failed = 0;

                // Number of cluster inserted
                ++nClusterInserted_;
                // Total number of particles inserted
                numberOfParticlesInserted_ += cluster.particleHandler.getSize();
                // This is the total mass composed by every single particle (not cluster!) inserted
                massInserted_ += cluster.particleHandler.getMass();
                // This is the total volume composed by every single particle (not cluster!) inserted
                volumeInserted_ += cluster.particleHandler.getVolume();
                logger(VERBOSE, "successfully placed a particle %, with position: % after % fails.", p0,
                       p0->getPosition(), failed);

                /* JMFT: The generateParticle() routine allocates memory, so we should
                 * free it here. (Don't worry, the particle will have been copied to the
                 * particleHandler by this point iff we want it.) */
                delete p0;

                break; // out of the 'placing' loop
            }
            else
            {
                failed++;
                logger(VERBOSE, "failed to place a cluster; have failed % times", failed);
            }

            if (failed > maxFailed_)
            {
                logger(VERBOSE, "failed too many times; giving up");
                break; // out of the 'placing' loop (and will leave the 'generating' loop too
            }
        }
        logger(VERBOSE, "failed % times, so breaking out of InsertionBoundary loop for this timestep.", failed);
    }
    // logger(INFO, "volumeInserted_ = %", volumeInserted_);
}


//!\brief this sets position and velocity of the cluster.
void RandomClusterInsertionBoundary::placeParticle(BaseParticle* p, RNG& random)
{
    if(getRandomised())
        random.randomise();
    Vec3D pos, vel;
    pos.X = random.getRandomNumber(posMin_.X, posMax_.X);
    pos.Y = random.getRandomNumber(posMin_.Y, posMax_.Y);
    pos.Z = random.getRandomNumber(posMin_.Z, posMax_.Z);
    vel.X = random.getRandomNumber(velMin_.X, velMax_.X);
    vel.Y = random.getRandomNumber(velMin_.Y, velMax_.Y);
    vel.Z = random.getRandomNumber(velMin_.Z, velMax_.Z);
    p->setPosition(pos);
    p->setVelocity(vel);
}

/*!
 * \details Returns the name of the object class
 * \return      the object's class' name, i.e. 'ClusterInsertionBoundary'
 */
std::string RandomClusterInsertionBoundary::getName() const
{
    return "RandomClusterInsertionBoundary";
}
