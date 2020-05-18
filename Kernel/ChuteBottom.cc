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


#include <algorithm>
#include "ChuteBottom.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Species/FrictionForceSpecies/SlidingFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Logger.h"

/*!
 * \details Default constructor. Calls the constructor() method.
 */
ChuteBottom::ChuteBottom()
{
    constructor();
}

/*!
 * \details Copy constructor with a DPMBase object as an argument. This constructor
 * basically 'upgrades' a DPMBase object to one of the ChuteBottom class.
 * NB: The copy constructor of DPMBase has to be called because the link from DPMBase 
 * to MercuryBase is virtual.
 * \param[in] other     object of DPMBase class
 */
ChuteBottom::ChuteBottom(const DPMBase& other)
        : DPMBase(other), Chute(other)
{
    constructor();
}

/*!
 * \details Copy constructor with a MercuryBase object as an argument. This constructor
 * basically 'upgrades' a MercuryBase object to one of the ChuteBottom class.
 * NB: The copy constructor of DPMBase has to be called because the link from DPMBase 
 * to MercuryBase is virtual.
 * \param[in] other     object of MercuryBase class
 */
ChuteBottom::ChuteBottom(const MercuryBase& other)
        : DPMBase(other), Chute(other)
{
    constructor();
}

/*!
 * \details Copy constructor with a Mercury3D object as an argument. This constructor
 * basically 'upgrades' a Mercury3D object to one of the ChuteBottom class.
 * NB: The copy constructor of DPMBase has to be called because the link from DPMBase 
 * to MercuryBase is virtual.
 * \param[in] other     object of Mercury3D class
 */
ChuteBottom::ChuteBottom(const Mercury3D& other)
        : DPMBase(other), Chute(other)
{
    constructor();
}

/*!
 * \details Copy constructor with a Chute object as an argument. This constructor
 * basically 'upgrades' a Chute object to one of the ChuteBottom class.
 * NB: The copy constructor of DPMBase has to be called because the link from DPMBase 
 * to MercuryBase is virtual.
 * \param[in] other     object of Chute class
 */
ChuteBottom::ChuteBottom(const Chute& other)
        : DPMBase(other), Chute(other)
{
    constructor();
}

/*!
 * \details 'normal' copy constructor
 * \param[in] other     ChuteBottom object to be copied
 * \todo IFCD: why are the name and file-stuff not set here, in contrast to constructor()?
 */
ChuteBottom::ChuteBottom(const ChuteBottom& other)
        : DPMBase(other), Chute(other),
          thickness_(other.thickness_), isBottomPeriodic_(other.isBottomPeriodic_)
{
    
}

/*!
 * \details constructor METHOD, which sets all chute bottom properties to something sensible.
 */
void ChuteBottom::constructor()
{
    setName("roughbottom");
    fStatFile.setFileType(FileType::NO_FILE); //set to 0 for no data creation
    dataFile.setFileType(FileType::ONE_FILE);
    restartFile.setFileType(FileType::ONE_FILE);
    eneFile.setFileType(FileType::NO_FILE);
    interactionFile.setFileType(FileType::NO_FILE);
    setThickness(2.4);
    setIsBottomPeriodic(true);
}

/*!
 * \details Creates a multilayer rough bottom as follows:
 *    1. Sets up a mini simulation, with a horizontal chute and particles poring in
 *       on the left (between 0 < Z < inflowHeight_) 
 *    2. Runs the simulation for 2000 time steps 
 *    3. Removes all particles from the system, EXCEPT those with a Z-position s.t.:
 *          hmax = height - maxInflowParticleRadius_;
 *          hmax - (thickness_ * maxInflowParticleRadius_) <=  Z-position  <= hmax
 *    4. Moves all remaining particles to the bottom and fixes them in space
 *    5. Transfers all particles to the Chute argument
 * 
 * \param[out] chute    The Chute object for which the multilayer rough bottom is created
 * \todo Make sure that the post-conditions are:
 * - N fixed particles in particleHandler, with z-coordinate < 0
 * - 1 wall, below the rough bottom
 * - no boundaries
 */
void ChuteBottom::makeRoughBottom(Chute& chute)
{
    // set up mini-simulation with particle inflow on the left end and horizontal chute
    // set all parameters that should be different from the original chute
    setChuteAngle(0.0);
    setInflowHeight(25. * getInflowParticleRadius());
    //~ setInflowHeight(45.*getInflowParticleRadius());
    // note: Changing the Inflow height was an attempt to make the bottom density homogeneous, but it did not have the desired effect
    setRoughBottomType(MONOLAYER_DISORDERED);
    setFixedParticleRadius(getInflowParticleRadius());
    
    auto species = dynamic_cast<LinearViscoelasticNormalSpecies*>(speciesHandler.getObject(0));
    if (species != nullptr)
    {
        const Mdouble mass = speciesHandler.getObject(0)->getMassFromRadius(
                0.5 * (getMinInflowParticleRadius() + getMaxInflowParticleRadius()));
        //this number is chosen to be consistent with the old implementation (with contained a bug and was replaced); originally, a factor of 10.0 (and a restitution of 0.2) was chosen.
        const Mdouble collisionTime = 13.2863 * species->getCollisionTime(mass);
        species->setCollisionTimeAndRestitutionCoefficient(collisionTime, 0.343008, mass);
        setTimeStep(0.02 * 9.41823 * collisionTime);
        logger(INFO, "Time step: %", getTimeStep());
        //logger(INFO, "Species: %", *species);
    }
    else
    {
        ///\todo IFCD: Only a warning, or an error?
        logger(WARN, "[ChuteBottom::makeRoughBottom()] species type does not allow setting the parameters.");
    }
    
    auto species2 = dynamic_cast<SlidingFrictionSpecies*>(speciesHandler.getObject(0));
    if (species2 != nullptr)
        species2->setSlidingFrictionCoefficient(0);
    
    // set the simulation to run for 2000 time steps
    setTimeMax(getTimeStep() * 2000);
    setSaveCount(100);
    
    // run the simulation
    solve();
    
    //Find the Z-position of the highest particle in the system
    Mdouble height = 0;
    for (BaseParticle* const p : particleHandler)
    {
        height = std::max(height, p->getPosition().Z);
    }
    
    // Next, all particles are removed from the system, except those with a Z-position s.t.:
    // hmax = height - maxInflowParticleRadius_;
    // hmax - thickness_ * maxInflowParticleRadius_ <  Z-position  < hmax
    // note that hmax = height - rMax
    logger(INFO, "[ChuteBottom::makeRoughBottom()] Thickness: %", thickness_);
    
    for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
    {
        if ((*it)->getPosition().Z < height - (1.0 + thickness_) * getMaxInflowParticleRadius() ||
            (*it)->getPosition().Z > height - getMaxInflowParticleRadius())
        {
            //delete particles outside the given range
            particleHandler.removeObject((*it)->getIndex());
            --it;
        }
    }
    //The remaining particles must be moved downwards and fixed.
    // They also must have the species as specified by the chute.
    for (BaseParticle* p : particleHandler)
    {
        p->move(Vec3D(0.0, 0.0, getMaxInflowParticleRadius() - height));
        p->fixParticle();
        p->setSpecies(chute.speciesHandler.getObject(0));
    }
    
    //copy the rough bottom over
    logger(INFO, "[ChuteBottom::makeRoughBottom()] Chute bottom finished, consisting of % particles",
           particleHandler.getNumberOfObjects());
    chute.particleHandler = particleHandler;
    
}

/*!
 * \details Creates the chute bottom, (periodic or solid) side walls and fills the chute
 * with particles (without using an insertion boundary).
 */
void ChuteBottom::setupInitialConditions()
{
    
    particleHandler.setStorageCapacity(static_cast<unsigned int>(std::min(
            3.0 * getXMax() * getYMax() * getZMax() / mathsFunc::cubic(2.0 * getInflowParticleRadius()), 1e6)));
    
    createBottom();
    
    /*!
     * \todo The createBottom() function also creates some walls and boundaries, 
     * but at slightly different locations. In this version they are removed and reset, 
     * but this is (in my opinion (Dinant)) not the correct way.
     */
    wallHandler.clear();
    boundaryHandler.clear();
    if (isBottomPeriodic_)
    {
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin() - getInflowParticleRadius()));
        wallHandler.copyAndAddObject(w0);
        PeriodicBoundary b0;
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
    }
    else
    {
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin() - getInflowParticleRadius()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
    }
    
    // add particles
    /*!
     * \todo Particles are created without insertion boundary... this reeks of double work 
     * (see also Dinant's todo in this same function (ChuteBottom::setupInitialConditions()).
     * (BvdH)
     */
    hGridActionsBeforeTimeLoop();
    unsigned int failed = 0;
    const unsigned int max_failed = 500;
    
    SphericalParticle inflowParticle_;
    inflowParticle_.setSpecies(speciesHandler.getObject(0));
    inflowParticle_.setHandler(&particleHandler);
    inflowParticle_.setOrientation({1, 0, 0, 0});
    inflowParticle_.setAngularVelocity(Vec3D(0.0, 0.0, 0.0));
    
    // try max_failed times to find new insertable particle
    // every time an insertion succeeds the counter failed is set to 0 again.
    while (failed <= max_failed)
    {
        inflowParticle_.setRadius(getFixedParticleRadius());
        //inflowParticle_.computeMass();
        
        // The position components are first stored in a Vec3D, because if you pass 
        // them directly into setPosition the compiler is allowed to change the order 
        // in which the numbers are generated
        Vec3D position;
        position.X = random.getRandomNumber(inflowParticle_.getRadius(), getXMax() - inflowParticle_.getRadius());
        position.Y = random.getRandomNumber(inflowParticle_.getRadius(), getYMax() - inflowParticle_.getRadius());
        position.Z = random.getRandomNumber(2 * inflowParticle_.getRadius(), getZMax() - inflowParticle_.getRadius());
        inflowParticle_.setPosition(position);
        inflowParticle_.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        //if the volume we want to insert the particle is still free, insert the particle
        if (checkParticleForInteraction(inflowParticle_))
        {
            particleHandler.copyAndAddObject(inflowParticle_);
            failed = 0;
        }
        else
        {
            failed++;
        }
    }
    //set_Nmax(particleHandler.getNumberOfObjects());
    logger(INFO, "[ChuteBottom::setupInitialConditions()] Number of particles created: %",
           particleHandler.getNumberOfObjects());
    
    //fix hgrid (there is still an issue when particles are polydispersed)
    //assume 1-2 levels are optimal (which is the case for mono and bidispersed) and set the cell size to min and max
    // !this is not optimal for polydispersed
    Mdouble minCell = 2. * std::min(getFixedParticleRadius(), getMinInflowParticleRadius());
    Mdouble maxCell = 2. * std::max(getFixedParticleRadius(), getMaxInflowParticleRadius());
    if ((minCell == maxCell) | (minCell == 0.))
        setHGridMaxLevels(1);
    else
        setHGridMaxLevels(2);
    //set_HGRID_cell_to_cell_ratio (1.0000000001*maxCell/minCell);
    //optimize number of buckets
    //setHGridNumberOfBucketsToPower(particleHandler.getNumberOfObjects() * 1.5);
    //end: fix hgrid
    
    //~ write(std::cout,false);
}

/*!
 * \details Performs all necessary actions before the next time step, which are
 * none.
 */
void ChuteBottom::actionsBeforeTimeStep()
{
}

/*!
 * \details Returns the thickness_ of the multilayer rough chute bottom. See also 
 * the documentation of Chute::createBottom(). 
 * \return      the thickness of the multilayer rough chute bottom
 */
Mdouble ChuteBottom::getThickness() const
{
    return thickness_;
}

/*!
 * \details Sets the thickness of the multilayer rough chute bottom. See also the 
 * documentation of Chute::createBottom(). 
 * \param[in] new_      the thickness_ of the multilayer rough chute bottom to be set
 */
void ChuteBottom::setThickness(Mdouble new_)
{
    logger.assert_always(new_ > 0.0, "[ChuteBottom::setThickness()] thickness % not positive.", new_);
    thickness_ = new_;
}

/*!
 * \details Returns if the bottom is a periodic one.
 * \return      TRUE if the bottom is periodic, FALSE if not
 */
bool ChuteBottom::getIsBottomPeriodic() const
{
    return isBottomPeriodic_;
}

/*!
 * \details Sets whether the bottom should be periodic.
 * \param[in] isBottomPeriodic  TRUE if the bottom should be periodic, FALSE if not
 */
void ChuteBottom::setIsBottomPeriodic(bool isBottomPeriodic)
{
    isBottomPeriodic_ = isBottomPeriodic;
}
