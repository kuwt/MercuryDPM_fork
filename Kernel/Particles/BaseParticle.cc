//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "BaseParticle.h"
#include "DPMBase.h"

/*!
 * \details default constructor, creates an Particle at (0,0,0) with radius, 
 * mass and inertia equal to 1
 */
BaseParticle::BaseParticle()
        : BaseInteractable()
{
    handler_ = nullptr;
    displacement_.setZero();
    radius_ = 1.0;
    invMass_ = 1.0;
    invInertia_ = MatrixSymmetric3D(1,0,0,1,0,1);

    periodicFromParticle_ = nullptr;
    isMPIParticle_ = false;
    isInMPIDomain_ = false;
    isInPeriodicDomain_ = false;
    isPeriodicGhostParticle_ = false;
    isMaserParticle_ = false;
    communicationComplexity_ = 0;
    periodicComplexity_ = std::vector<int>(0);
    previousPeriodicComplexity_ = std::vector<int>(0);
#ifdef CONTACT_LIST_HGRID
    firstPossibleContact = nullptr;
#endif
    hGridNextObject_ = nullptr;
    hGridPrevObject_ = nullptr;

    hGridCell.setHGridLevel(99999);
    hGridCell.setHGridX(99999);
    hGridCell.setHGridY(99999);
    hGridCell.setHGridZ(99999);

    info_ = std::numeric_limits<double>::quiet_NaN();

    logger(DEBUG, "BaseParticle::BaseParticle() finished");
}

/*!
 * \details Constructor that copies most of the properties of the given particle.
 *          Please note that not everything is copied, for example the position 
 *          in the HGrid is not determined yet by the end of this constructor. 
 *          It also does not copy the interactions and the pointer to the handler
 *          that handles this particle. Use with care.
 * \param[in,out] p  Reference to the BaseParticle this one should become a copy of.
 */
BaseParticle::BaseParticle(const BaseParticle &p)
: BaseInteractable(p)
{
    handler_ = nullptr;
    displacement_ = p.displacement_;
    radius_ = p.radius_;
    invMass_ = p.getInvMass();
    invInertia_ = p.getInvInertia();

    hGridNextObject_ = nullptr;
    hGridPrevObject_ = nullptr;

    hGridCell.setHGridLevel(p.getHGridLevel());
    hGridCell.setHGridX(99999);
    hGridCell.setHGridY(99999);
    hGridCell.setHGridZ(99999);

    periodicFromParticle_ = p.periodicFromParticle_;
    isMPIParticle_ = p.isMPIParticle_;
    isInMPIDomain_ = p.isInMPIDomain_;
    isInPeriodicDomain_ = p.isInPeriodicDomain_;
    isMaserParticle_ = p.isMaserParticle_;
    isPeriodicGhostParticle_ = p.isPeriodicGhostParticle_;
    communicationComplexity_ = p.communicationComplexity_;
    periodicComplexity_ = p.periodicComplexity_;
    previousPeriodicComplexity_ = p.previousPeriodicComplexity_;
#ifdef CONTACT_LIST_HGRID
    firstPossibleContact = nullptr;
#endif

    info_ = p.info_;
    logger(DEBUG, "BaseParticle::BaseParticle(BaseParticle &p) finished");
}

BaseParticle::BaseParticle(const ParticleSpecies* s)
: BaseParticle()
{
    setSpecies(s);
#ifdef CONTACT_LIST_HGRID
    firstPossibleContact = nullptr;
#endif
    logger(DEBUG, "BaseParticle::BaseParticle(BaseSpecies &s) finished");
}

/*!
 * \details Destructor. It asks the ParticleHandler to check if this was the 
 *          smallest or largest particle and adjust itself accordingly.
 */
BaseParticle::~BaseParticle()
{

    if (getHandler() != nullptr)
    {
        getHandler()->checkExtremaOnDelete(this);
        if (isFixed())
            getHandler()->removedFixedParticle();
    }
    logger(DEBUG, "BaseParticle::~BaseParticle() of particle % finished.", getId());

}

/*!
 * \details Copy method. Uses copy constructor to create a copy on the heap. 
 *          Useful for polymorphism.
 * \return pointer to the particle's copy
 */
BaseParticle* BaseParticle::copy() const
{
    return new BaseParticle(*this);
}

/*!
 * \details Returns the volume of the BaseParticle, which is calculated using
 *          its number of dimensions and radius.
 * \return The actual volume of this BaseParticle.
 */
Mdouble BaseParticle::getVolume() const
{
    if (handler_ == nullptr)
    {
        logger(ERROR, "[BaseParticle::getVolume] no particle handler specified");
        return 0;
    }
    switch (getParticleDimensions())
    {
    case 3:
        return (4.0 / 3.0 * constants::pi * radius_ * radius_ * radius_);
    case 2:
        return (constants::pi * radius_ * radius_);
    case 1:
        return (2.0 * radius_);
    default:
        logger(ERROR, "[BaseParticle::getVolume] dimension of the particle is not set");
        return 0;
    }
}

/*!
 * \details Fixes a BaseParticle by setting its inverse mass and inertia and velocities to zero.
 */
void BaseParticle::fixParticle()
{
    invMass_ = 0.0;
    invInertia_ = MatrixSymmetric3D(0,0,0,0,0,0);
    setVelocity(Vec3D(0.0, 0.0, 0.0));
    setAngularVelocity(Vec3D(0.0, 0.0, 0.0));
    if (getHandler())
        getHandler()->addedFixedParticle();
}

/*!
 * \details Checks whether a BaseParticle is fixed or not, by checking its inverse Mass.
 * \return TRUE if particle is fixed, i.e. if the inverse mass (invMass_) is 0.
 */
bool BaseParticle::isFixed() const
{
    return (invMass_ == 0.0);
}

bool BaseParticle::isMPIParticle() const
{
    //make mpi-dependent so the compiler can optimise
#ifdef MERCURY_USE_MPI
    return isMPIParticle_;
#else
    return false;
#endif
}

void BaseParticle::setMPIParticle(bool flag)
{
    isMPIParticle_ = flag;
}

void BaseParticle::setCommunicationComplexity(unsigned complexity)
{
    communicationComplexity_ = complexity;
}

unsigned BaseParticle::getCommunicationComplexity()
{
    return communicationComplexity_;
}

void BaseParticle::setPeriodicComplexity(std::vector<int> complexity)
{
    periodicComplexity_ = complexity;
}



void BaseParticle::setPeriodicComplexity(int index, int value)
{
    //hack: generally you'd add particles after declaring the boundaries
    //but no official programming guildelines rules have been setup for that
    //So incase that doesnt happen we need to resize this periodicComplexity
    if (periodicComplexity_.size() == 0)
    {
        int numberOfPeriodicBoundaries = getHandler()->getDPMBase()->periodicBoundaryHandler.getSize();
        if (numberOfPeriodicBoundaries > 0 )
        {
            //First initialisation of the periodic complexity assumes the particle is completely
            //within the real domain
            periodicComplexity_ = std::vector<int>(numberOfPeriodicBoundaries,2);
        }
    }

    periodicComplexity_[index] = value;
}

const std::vector<int>& BaseParticle::getPeriodicComplexity()
{
    //TODO resolve this hack
    //hack: generally you'd add particles after declaring the boundaries
    //but no official programming guildelines rules have been setup for that
    //So incase that doesnt happen we need to resize this periodicComplexity
    if (periodicComplexity_.size() == 0)
    {
        int numberOfPeriodicBoundaries = getHandler()->getDPMBase()->periodicBoundaryHandler.getSize();
        if (numberOfPeriodicBoundaries > 0 )
        {
            periodicComplexity_.resize(numberOfPeriodicBoundaries,0);
        }
    }

    return periodicComplexity_;
}

int BaseParticle::getPeriodicComplexity(int index)
{
    //hack: generally you'd add particles after declaring the boundaries
    //but no official programming guildelines rules have been setup for that
    //So incase that doesnt happen we need to resize this periodicComplexity
    if (periodicComplexity_.size() == 0)
    {
        int numberOfPeriodicBoundaries = getHandler()->getDPMBase()->periodicBoundaryHandler.getSize();
        if (numberOfPeriodicBoundaries > 0 )
        {
            periodicComplexity_.resize(numberOfPeriodicBoundaries,0);
        }
    }

    return periodicComplexity_[index];
}

void BaseParticle::setPreviousPeriodicComplexity(std::vector<int> complexity)
{
    previousPeriodicComplexity_ = complexity;
}

const std::vector<int>& BaseParticle::getPreviousPeriodicComplexity() const
{
    return previousPeriodicComplexity_;
}

bool BaseParticle::isInMPIDomain()
{
    return isInMPIDomain_;
}

void BaseParticle::setInMPIDomain(bool flag)
{
    isInMPIDomain_ = flag;
}


bool BaseParticle::isInPeriodicDomain() const
{
    return isInPeriodicDomain_;
}

void BaseParticle::setInPeriodicDomain(bool flag)
{
    isInPeriodicDomain_ = flag;
}

bool BaseParticle::isPeriodicGhostParticle() const
{
    return isPeriodicGhostParticle_;
}

void BaseParticle::setPeriodicGhostParticle(bool flag)
{
    isPeriodicGhostParticle_ = flag;
}

bool BaseParticle::isMaserParticle() const
{
    return isMaserParticle_;
}

void BaseParticle::setMaserParticle(bool flag)
{
    isMaserParticle_ = flag;
}

/*!
 * \details Unfixes the particle by computing the Particles mass and inertia, using the 
 * species and radius.
 */

void BaseParticle::unfix()
{
    invMass_ = 1.0;
    getSpecies()->computeMass(this);
    if (getHandler())
        getHandler()->removedFixedParticle();
}

/*!
 * \details BaseParticle print method, which accepts an os std::ostream as 
 *          input. It prints human readable BaseParticle information to the 
 *          std::ostream.
 * \param[in,out] os    stream to which the info is written
 */
void BaseParticle::write(std::ostream& os) const
{
    BaseInteractable::write(os);
    os << " radius " << radius_
       << " invMass " << invMass_;
       //invMass_ is a computed value, but needs to be stored to see if a particle is fixed
}

/*!
 * \details Returns the name of the object; in this case 'BaseParticle'.
 * \return The object name.
 */
std::string BaseParticle::getName() const
{
    return "BaseParticle";
}

void BaseParticle::setInfo(Mdouble info)
{
    info_ = info; 
}

Mdouble BaseParticle::getInfo() const
{
    if (std::isnan(info_))
        return getSpecies()->getId();
    else
        return info_;
}

/*!
 * \details Particle read function. Has an std::istream as argument, from which 
 *          it extracts the radius_, invMass_ and invInertia_, respectively. 
 *          From these the mass and inertia are deduced. An additional set of
 *          properties is read through the call to the parent's method
 *          BaseInteractable::read().
 * \param[in,out] is    input stream with particle properties.
 */
void BaseParticle::read(std::istream& is)
{
    BaseInteractable::read(is);
    std::string dummy;
    is >> dummy >> radius_ >> dummy >> invMass_;// >> dummy >> invInertia_;
}

/*!
 * \details This is the previously used version of the read function. Now just kept
 *          for legacy purposes. 
 * \deprecated Should be gone in Mercury 2.0. Use BaseParticle::read() instead.
 */
void BaseParticle::oldRead(std::istream& is)
{
    logger(DEBUG, "reading particle old-style");
    static unsigned int id = 0;
    unsigned int indSpecies = 0;
    unsigned int numberOfContacts = 0;
    Vec3D orientation;
    Vec3D position;
    Vec3D velocity;
    Vec3D angularVelocity;
    double invInertiaScalar;
    double dummy = 0;
    is >> position >> velocity >> radius_ >> orientation >> angularVelocity;
    is >> invMass_ >> invInertiaScalar >> numberOfContacts;
    ///\todo incorporate contact information
    for (unsigned int i = 0; i < 12 * numberOfContacts; ++i)
    {
        is >> dummy;
    }
    is >> indSpecies;
    setPosition(position);
    setVelocity(velocity);
    setOrientation(Quaternion(orientation));
    setAngularVelocity(angularVelocity);
    invInertia_.XX = invInertiaScalar;
    invInertia_.YY = invInertiaScalar;
    invInertia_.ZZ = invInertiaScalar;
    setIndSpecies(indSpecies);
    setId(id);
    setIndex(id);
    id++;
}

/*!
 * \details Adds the particle's HGridLevel_ and HGRid x/y/z positions to an 
 *          std::ostream.
 * \param[in,out] os    the ostream which has the mentioned properties added.
 */
void BaseParticle::printHGrid(std::ostream& os) const
{
    os << "Particle( HGRID_Level:" << hGridCell.getHGridLevel()
        << ", HGRID_x:" << hGridCell.getHGridX()
        << ", HGRID_y:" << hGridCell.getHGridY()
        << ", HGRID_z:" << hGridCell.getHGridZ()
        << ")";
}

/*!
 * \details Returns the particle's HGridLevel_
 * \return the HGridLevel_
 */
unsigned int BaseParticle::getHGridLevel() const
{
    return hGridCell.getHGridLevel();
}

/*!
 * \details Returns the next object in the particle's HGrid cell
 * \return pointer to the next object in the particle's HGrid cell
 */
BaseParticle* BaseParticle::getHGridNextObject() const
{
    return hGridNextObject_;
}

/*!
 * \details Returns the previous object in the particle's HGrid cell
 * \return pointer to the previous object in the particle's HGrid
 * cell
 */
BaseParticle* BaseParticle::getHGridPrevObject() const
{
    return hGridPrevObject_;
}

#ifdef CONTACT_LIST_HGRID

/*!
 * \details 
 */
PossibleContact* BaseParticle::getFirstPossibleContact() const
{
    return firstPossibleContact;
}
#endif

/*!
 * \details Get the particle's HGrid cell's X-coordinate
 * \return the particle's HGrid cell's X-coordinate
 */
int BaseParticle::getHGridX() const
{
    return hGridCell.getHGridX();
}

/*!
 * \details Get the particle's HGrid cell's Y-coordinate
 * \return the particle's HGrid cell's Y-coordinate
 */
int BaseParticle::getHGridY() const
{
    return hGridCell.getHGridY();
}

/*!
 * \details Get the particle's HGrid cell's Z-coordinate
 * \return the particle's HGrid cell's Z-coordinate
 */
int BaseParticle::getHGridZ() const
{
    return hGridCell.getHGridZ();
}

/*!
 * \details Returns the particle's inverse inertia
 * \return the particles invInertia_
 */
MatrixSymmetric3D BaseParticle::getInvInertia() const
{
    return invInertia_;
}

/*!
 * \details Returns the particle's inverse mass
 * \return the particle's invMass_
 */
Mdouble BaseParticle::getInvMass() const
{
    return invMass_;
}

/*!
 * \details Calculates the particle's kinetic energy
 * \return the particle's kinetic energy
 */
Mdouble BaseParticle::getKineticEnergy() const
{
    if (isFixed())
        return 0.0;
    else
        return 0.5 * getMass() * getVelocity().getLengthSquared();
}

Mdouble BaseParticle::getRotationalEnergy() const
{
    if (isFixed())
        return 0.0;
    else
        return 0.5 * Vec3D::dot(getAngularVelocity(),getInertia() * getAngularVelocity());
}

/*!
 * \details Returns the mass of the particle
 * \return the mass of the particle
 */
Mdouble BaseParticle::getMass() const
{
    return 1.0/invMass_;
}

MatrixSymmetric3D BaseParticle::getInertia() const
{
    return MatrixSymmetric3D::inverse(invInertia_);
}

/*!
 * \details Returns a pointer to the 'original' particle if the current one is
 * a 'periodic copy' used for a periodic boundary implementation.
 * \return pointer to original particle
 */
BaseParticle* BaseParticle::getPeriodicFromParticle() const
{
    return periodicFromParticle_;
}

/*!
 * \details Returns the particle's radius
 * \return particle's radius
 */
Mdouble BaseParticle::getRadius() const
{
    return radius_;
}

/*!
 * \details Calculates the interaction radius of the particle (when it comes to 
 * interaction with other particles), including the effect of a possible additional
 * 'interaction distance' besides the 'normal' radius. The interaction radius
 * differs from the radius_, for example, when dealing with wet particles (i.e. 
 * particles with an additional liquid layer, which is dealt with in the particle's 
 * species).
 * \return the particle's interaction radius for particle-particle interaction
 */
Mdouble BaseParticle::getInteractionRadius() const
{
    return radius_ + getSpecies()->getInteractionDistance() * 0.5;
}

/*!
 * \details The interaction radius of the particle (when it comes to interaction
 * with walls). See also BaseParticle::getInteractionRadius().
 * \return the particle's interaction radius for particle-wall interaction
 */
Mdouble BaseParticle::getWallInteractionRadius() const
{
    return getInteractionRadius() + getSpecies()->getInteractionDistance() * 0.5;
}

/*!
 * \details Returns the particle's displacement_, which is the difference between
 * the current particle's position and its position in the previous time step.
 * \return (reference to) the particle displacement vector
 */
const Vec3D& BaseParticle::getDisplacement() const
{
    return displacement_;
}

/*!
 * \details Returns the particle's position in the previous time step.
 * \return (reference to) the previous position of the particle
 */
const Vec3D& BaseParticle::getPreviousPosition() const
{
    return previousPosition_;
}

/*!
 * \todo Rewrite, redefine (TW). Is only used in StatisticsVector.hcc, consider 
 * moving to that class.
 */
const Vec3D BaseParticle::getDisplacement2(Mdouble xmin, Mdouble xmax, Mdouble ymin, Mdouble ymax, Mdouble zmin, Mdouble zmax, Mdouble t) const
{
    Vec3D disp = getPosition() - getPreviousPosition();
    if (xmax > xmin && fabs(disp.X) > .5 * (xmax - xmin))
    {
        if (disp.X > 0)
            disp.X -= xmax - xmin;
        else
            disp.X += xmax - xmin;
    }
    if (ymax > ymin && fabs(disp.Y) > .5 * (ymax - ymin))
    {
        if (disp.Y > 0)
            disp.Y -= ymax - ymin;
        else
            disp.Y += ymax - ymin;
    }
    if (zmax > zmin && fabs(disp.Z) > .5 * (zmax - zmin))
    {
        if (disp.Z > 0)
            disp.Z -= zmax - zmin;
        else
            disp.Z += zmax - zmin;
    }
    disp /= t;
    return disp;
}

/*!
 * \details
 */
void BaseParticle::setInertia()
{

}

/*!
 * \details Sets the particle's inertia and invInertia_.
 * \param[in] newInertia  the new inertia to be set.
 */
void BaseParticle::setInertia(const MatrixSymmetric3D inertia)
{
    invInertia_ = MatrixSymmetric3D::inverse(inertia);
}

void BaseParticle::setInverseInertia(const MatrixSymmetric3D inverseInertia)
{
    invInertia_ = inverseInertia;
}

/*!
 * \details Sets the inertia to 1e20 and the invInertia_ (which is actually
 * used in the calculations) to 0.
 */
void BaseParticle::setInfiniteInertia()
{
    invInertia_.setZero();
} //> i.e. no rotations

/*!
 * \details Lets you set which particle this one is actually a periodic copy of
 * (used in periodic boundary condition implementations). 
 * \param[in] p  pointer to the 'original' particle this one is a periodic 
 * copy of.
 */
void BaseParticle::setPeriodicFromParticle(BaseParticle* p)
{
    periodicFromParticle_ = p;
}

/*!
 * \details Set the x-index of the particle's HGrid cell position
 * \param[in] x  x-index of particle's HGrid cell
 */
void BaseParticle::setHGridX(const int x)
{
    hGridCell.setHGridX(x);
}

/*!
 * \details Set the y-index of the particle's HGrid cell position
 * \param[in] y  y-index of particle's HGrid cell
 */
void BaseParticle::setHGridY(const int y)
{
    hGridCell.setHGridY(y);
}

/*!
 * \details Set the z-index of the particle's HGrid cell position
 * \param[in] z  z-index of particle's HGrid cell
 */
void BaseParticle::setHGridZ(const int z)
{
    hGridCell.setHGridZ(z);
}

/*!
 * \details Sets the particle's HGridLevel_
 * \param[in] level  the particle's HGridLevel_
 */
void BaseParticle::setHGridLevel(const unsigned int level)
{
    hGridCell.setHGridLevel(level);
}

/*!
 * \details Sets the pointer to the next object in the particle's HGrid cell
 * \param[in] p     pointer to the next object
 */
void BaseParticle::setHGridNextObject(BaseParticle* p)
{
    hGridNextObject_ = p;
}

/*!
 * \details Sets the pointer to the previous object in the particle's HGrid cell
 * \param[in] p     pointer to the previous object
 */
void BaseParticle::setHGridPrevObject(BaseParticle* p)
{
    hGridPrevObject_ = p;
}

/*!
 * \details 
 */
#ifdef CONTACT_LIST_HGRID

void BaseParticle::setFirstPossibleContact(PossibleContact* PC)
{
    firstPossibleContact = PC;
}
#endif

/*!
 * \details Sets the radius of the particle, and from that computes the new mass
 * (using its species) and checks whether it now is either the smallest or biggest 
 * particle in its ParticleHandler. 
 * \param[in] radius    the new radius
 */
void BaseParticle::setRadius(const Mdouble radius)
{
    radius_ = radius;
    if (getHandler())
    {
        getSpecies()->computeMass(this);
        getHandler()->checkExtrema(this);
    }
}

/*!
 * \details Sets the mass of the particle
 * \param[in] mass  the new particle's  mass
 */
void BaseParticle::setMass(const Mdouble mass)
{
    logger(WARN, "WARNING: Do not use particle->setMass, instead use "
           "particleSpecies->computeMass, since this function can cause "
           "inconsistencies between the mass, density and radius of this particle!");
    if (mass > 0.0 && !isFixed())
    {
        invMass_ = 1.0 / mass;
    } else {
        logger(ERROR, "Error in BaseParticle::setMass, the given mass to be set must be positive.");
    }
}

/*!
 * \details Sets the mass of the particle
 * \param[in] mass  the new particle's  mass
 */
void BaseParticle::setMassForP3Statistics(const Mdouble mass)
{
    if (mass > 0.0 && !isFixed())
    {
        invMass_ = 1.0 / mass;
    } else {
        logger(ERROR, "Error in BaseParticle::setMass, the given mass to be set must be positive.");
    }
} 

/*!
 * \details This is used to set the particle displacement_ 
 * \param[in] disp  the displacement vector
 */
void BaseParticle::setDisplacement(const Vec3D& disp)
{
    displacement_ = disp;
}

/*!
 * \details This is used to set the particle's previous position
 * \param[in] pos   the particle's previous position vector.
 */
void BaseParticle::setPreviousPosition(const Vec3D& pos)
{
    previousPosition_ = pos;
}

/*!
 * \details Lets you add a vector to the particle's previousPosition_ vector.
 * \param[in] posMove   the vector to be added to the current previousPosition_
 * vector.
 */
void BaseParticle::movePrevious(const Vec3D& posMove)
{
    previousPosition_ += posMove;
}

/*!
 * \details increases the the particle's velocity_ (BaseInteractable member) 
 * by adding the given vector.
 * \param[in] vel   vector to be added to the velocity_
 */
void BaseParticle::accelerate(const Vec3D& vel)
{
    addVelocity(vel);
}

/*!
 * \details increases the particle's angularVelocity_ (BaseInteractable member)
 * by adding the given vector.
 * \param[in] angVel    vector to be added to the angularVelocity_
 */
void BaseParticle::angularAccelerate(const Vec3D& angVel)
{
    addAngularVelocity(angVel);
}

/*!
 * \details Lets you add a vector to the particle's displacement_ vector.
 * \param[in] addDisp   vector to be added.
 */
void BaseParticle::addDisplacement(const Vec3D& addDisp)
{
    displacement_ += addDisp;
}

/*!
 * \details Assigns the particle to a ParticleHandler, and assigns a species to 
 * it based on the particles indSpecies_ (BaseInteractable data member).
 * \param[in] handler   pointer to the ParticleHandler
 */
void BaseParticle::setHandler(ParticleHandler* handler)
{
    handler_ = handler;
    setSpecies(getHandler()->getDPMBase()->speciesHandler.getObject(getIndSpecies()));
}

/*!
 * \details Returns the particle's ParticleHandler
 * \return pointer to the particle's ParticleHandler; null if handler not yet specified
 */
ParticleHandler* BaseParticle::getHandler() const
{
    return handler_;
}

/*!
 * \details Creates/updates a BaseInteraction object, treating the interaction between
 * this particle and a given one, in case there is an overlap between the two.
 * \param[in] P             particle to check the interaction with
 * \param[in] timeStamp     time stamp to be assigned to the interaction object (i.e., 
 * the current time) 
 * \param[in,out] interactionHandler    BaseInteraction container from where the 
 * interaction is retrieved, and to which it is assigned (if it is a new interaction). 
 * \return the pointer to the interaction object (if the particles overlap), or 0 
 * (if they don't overlap).
 */

std::vector<BaseInteraction*> BaseParticle::getInteractionWith(BaseParticle* const P, const unsigned timeStamp,
                                                  InteractionHandler* const interactionHandler)

{
    //get the normal (from P away from the contact)
    const Vec3D branchVector = P->getPosition() - getPosition();
    //Get the square of the distance between particle i and particle j
    const Mdouble distanceSquared = Vec3D::getLengthSquared(branchVector);
    const Mdouble sumOfInteractionRadii = P->getInteractionRadius() + getInteractionRadius();
    std::vector<BaseInteraction*> interactions;
    if (distanceSquared < (sumOfInteractionRadii * sumOfInteractionRadii))
    {
        BaseInteraction* const C = interactionHandler->getInteraction(P, this, timeStamp);
        const Mdouble distance = std::sqrt(distanceSquared);
        C->setNormal(branchVector / distance);
        C->setOverlap(P->getRadius() + getRadius() - distance);
        C->setDistance(distance);
        C->setContactPoint(P->getPosition() - (P->getRadius() - 0.5 * C->getOverlap()) * C->getNormal());
        ///\todo We should consider setting the contact point to \author weinhartt
        //Mdouble ratio=P->getRadius()/(getRadius()+P->getRadius());
        //C->setContactPoint(P->getPosition() - (P->getRadius() - ratio * C->getOverlap()) * C->getNormal());
        interactions.push_back(C);
    }
    return interactions;
}

/*!
 * \details First step of Velocity Verlet integration (see also 
 * http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).
 * \param[in] time          current time
 * \param[in] timeStep      current time step
 */
void BaseParticle::integrateBeforeForceComputation(double time, double timeStep)
{
    ///\todo If the position is described by the user, one should also call
    ///BaseInteractable::integrateBeforeForceComputation. To check if it works
    ///correctly, remove the p0.fixParticle() line from the DrivenParticleUnitTest
    ///\author irana
    if (getInvMass() == 0.0)
    {
        BaseInteractable::integrateBeforeForceComputation(time, timeStep);
    }
    else
    {
#ifdef MERCURY_USE_MPI
        //For periodic particles in paralle the previous position is required
        setPreviousPosition(getPosition());
#endif
        accelerate(getForce() * getInvMass() * 0.5 * timeStep);
        setDisplacement(getVelocity() * timeStep);
        move(getDisplacement());
        getHandler()->getDPMBase()->hGridUpdateMove(this, getDisplacement().getLength());
        if (getHandler()->getDPMBase()->getRotation())
        {
            angularAccelerate(getOrientation().rotateInverseInertiaTensor(getInvInertia())* getTorque() * 0.5 * timeStep);
            //apply to rotation quaternion q: q = normalise(q + \tilde{C}\omega*timeStep) (see Wouter's notes)
            rotate(getAngularVelocity() * timeStep);
        }
    }
}

/*!
 * \details Second step of Velocity Verlet integration (see also 
 * http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet).
 * \param[in] time      current time
 * \param[in] timeStep  current time step
 */
void BaseParticle::integrateAfterForceComputation(double time, double timeStep)
{
    if (getInvMass() == 0.0)
    {
	//Updates a baseParticle with a prescribed motion
        BaseInteractable::integrateAfterForceComputation(time, timeStep);
    }
    else
    {
        accelerate(getForce() * getInvMass() * 0.5 * timeStep);
        if (getHandler()->getDPMBase()->getRotation())
        {
            angularAccelerate(getOrientation().rotateInverseInertiaTensor(getInvInertia()) * getTorque() * 0.5 * timeStep);
        }
    }
}

/*!
 * \details Returns the amount of dimensions of the particle (2 or 3, basically)
 * \return the number of dimension of the particle
 */
unsigned int BaseParticle::getParticleDimensions() const
{
    return getHandler()->getDPMBase()->getParticleDimensions();
}

/*!
 * \details Set the particle's species and species' index. Logs a warning if
 * no ParticleHandler is assigned. 
 * \param[in] indSpecies    The index of the species in the SpeciesHandler.
 */
void BaseParticle::setIndSpecies(unsigned int indSpecies)
{
    if (handler_ != nullptr)
    {
        BaseInteractable::setIndSpecies(indSpecies);
	setSpecies(handler_->getDPMBase()->speciesHandler.getObject(indSpecies));
        ///\todo TW do we have to update the species stored in the interactions here?
    }
    else
    {
        BaseInteractable::setIndSpecies(indSpecies);
        logger(ERROR, "setIndSpecies called on a particle with no particle handler.\n"
               "Therefore I can't request the given species from the species handler.\n"
               " PartID = %", getId());
    }
}

/*!
 * \details Sets the particle's species. If this particle does not have a 
 *          handler yet, this function also assigns the ParticleHandler in the 
 *          same DPMBase as the SpeciesHandler of the given species as its handler.
 * \param[in] species   pointer to the ParticleSpecies object, to be set as the 
 *                      particle's species.
 */
void BaseParticle::setSpecies(const ParticleSpecies* species)
{
    BaseInteractable::setSpecies(species);
    ///\todo TW should we chaeck here if we have the right kind of species for the right kind of particle?
    //set pointer to the ParticleHandler handler_, which is needed to retrieve 
    //species information
    //\todo maybe these if statements should throw warnings
    if (handler_ == nullptr)
    {
        SpeciesHandler* sH = species->getHandler();
        DPMBase* dB = sH->getDPMBase();
        if (dB != nullptr)
        {
            setHandler(&dB->particleHandler);
        }
    }
}

unsigned BaseParticle::getNumberOfFieldsVTK() const {
    return 0;
}

std::string BaseParticle::getTypeVTK(unsigned i) const {
    return "";
}

std::string BaseParticle::getNameVTK(unsigned i) const {
    return "";
}

std::vector<Mdouble> BaseParticle::getFieldVTK(unsigned i) const {
    return std::vector<Mdouble>();
}

Vec3D BaseParticle::getAxes() {return Vec3D(0,0,0);}

double BaseParticle::getExponentEps1() {return 0;}

double BaseParticle::getExponentEps2() {return 0;}

bool BaseParticle::isInContactWith(const BaseParticle* const P) const
{
    if (P->getName() != "Superquadric")
    {
        return Vec3D::getDistanceSquared(getPosition(), P->getPosition()) <
               mathsFunc::square(getInteractionRadius() + P->getInteractionRadius());
    }
    return P->isInContactWith(this);
}

const HGridCell& BaseParticle::getHGridCell() const
{
    return hGridCell;
}
