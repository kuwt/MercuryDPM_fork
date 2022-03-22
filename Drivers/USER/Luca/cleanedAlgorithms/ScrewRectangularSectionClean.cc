#include "ScrewRectangularSectionClean.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include "WallHandler.h"
#include "DPMBase.h"

///\todo include this drawing in Doxygen comments
/*
 -- BLADE SECTION --
 / / / / / / / / / / / / / / / / / / / / / / / / / /
 / / / / / ............................. / / / / / /
 / / / / / :       :            :       :/ / / / / /
 / / / / / :   2   :      1     :   2   :/ / / / / /
 / / / / / :.......:____________:.......:/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :   3   |      |     |   3   :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 ..................|   E  |  E  |...................
           :       |      |     |       :
      4    :   5   |      |     |   5   :     4
 __________________/      |     \___________________
           :              |             :           
      E    :       E      |     E       :     E
           :              |             :
 
 REGION 1: BLADE SIDE COLLISION
 REGION 2: EDGE COLLISION
 REGION 3: BLADE SURFACE COLLISION
 REGION 4: SHAFT COLLISION
 REGION 5: CORNER COLLISION (ROUNDED EDGE)
 REGION E: UNACCESSIBLE REGION (COLLISION ERROR)
 ELSE: NO COLLISION
 */

///\details Default constructor: set all variables to sensible defaults
ScrewRectangularSectionClean::ScrewRectangularSectionClean()
{
    startPosition_.setZero();
    length_ = 1.0;
    rMax_ = 1.0;
    rMin_ = 0.5 * rMax_;
    numberOfTurns_ = 1.0;
    omega_ = 1.0;
    thickness_ = 0.0;
    rightHandedness_ = true;
    
    logger(DEBUG, "ScrewRectangularSectionClean() constructor finished.");
}

/*!
 * \details Copy constructor: copy all values of the given screw into the constructed screw
 * \param[in] other screw that should be copied
 */
ScrewRectangularSectionClean::ScrewRectangularSectionClean(const ScrewRectangularSectionClean& other)
        : BaseWall(other)
{
    // user defined variables
    startPosition_ = other.startPosition_;
    length_ = other.length_;
    rMax_ = other.rMax_;
    rMin_ = other.rMin_;
    numberOfTurns_ = other.numberOfTurns_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    rightHandedness_ = other.rightHandedness_;
    
    logger(DEBUG, "ScrewRectangularSectionClean(const ScrewRectangularSectionClean&) copy constructor finished.");
}

/*!
 * Constructor that takes parameters for all screw properties, checks if they are valid and then sets all screw
 * properties
 * \param[in] start The centre of the lower end of the screw (3D vector)
 * \param[in] length The length of the screw (Mdouble > 0)
 * \param[in] bladeRadius The outer screw blade radius (Mdouble > 0)
 * \param[in] shaftRadius The screw shaft radius (0 <= Mdouble < bladeRadius)
 * \param[in] numberOfTurns The number of turns the screw contains (Mdouble > 0)
 * \param[in] omega Rotation speed in rad per time unit (Mdouble)
 * \param[in] thickness The thickness of the screw blade (Mdouble >= 0)
 * \param[in] rightHandedness The right handedness of the screw, i.e. the direction of the screw-blade (bool)
 */
ScrewRectangularSectionClean::ScrewRectangularSectionClean(const Vec3D& start, Mdouble length, Mdouble bladeRadius,
                                                           Mdouble shaftRadius, Mdouble numberOfTurns,
                                                           Mdouble omega, Mdouble thickness, bool rightHandedness)
{
    // check if the quantities are correct
    logger.assert_always(length > 0, "\nERROR: INVALID LENGTH.");
    logger.assert_always(bladeRadius > 0, "\nERROR: INVALID BLADE RADIUS.");
    logger.assert_always(shaftRadius < bladeRadius, "\nERROR: SHAFT RADIUS LARGER THAN BLADE RADIUS.");
    logger.assert_always(shaftRadius >= 0, "\nERROR: SHAFT RADIUS NEGATIVE");
    logger.assert_always(numberOfTurns > 0, "\nERROR: INVALID NUMBER OF SCREW TURNS.");
    logger.assert_always(thickness < length && thickness >= 0, "\nERROR: INVALID THICKNESS.");

    // user defined variables
    startPosition_ = start;
    length_ = length;
    rMax_ = bladeRadius;
    rMin_ = shaftRadius;
    numberOfTurns_ = numberOfTurns;
    omega_ = omega;
    thickness_ = thickness;
    rightHandedness_ = rightHandedness;
    
    logger(INFO, "\n\nScrew parameters set.\n"
                   "User defined variables:\n"
                   "Start: %\n Length: %\n Blade radius: %\n Shaft radius: %\n, Number of turns: %\n"
                   "Angular velocity: %\n Thickness: %\n Initial angular offset: %\n Right handedness: %\n"
                   "Derived variables:\n Half thickness: %\n Pitch length: %",
           startPosition_, length_, rMax_, rMin_, numberOfTurns_, omega_, thickness_, getOrientation(),
           rightHandedness_, thickness_ / 2, length_ / numberOfTurns_, getOrientation().getAxis().Z);

    logger(DEBUG, "ScrewRectangularSectionClean(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, bool) "
            "constructor finished.");
}

/*!
 * \details Destructor, just prints a message if the logger is set to debug mode
 */
ScrewRectangularSectionClean::~ScrewRectangularSectionClean()
{
    logger(DEBUG, "~ScrewRectangularSectionClean() finished, destroyed the screw.");
}

/*!
 * \details Function that makes a copy of the current screw, by calling new. Note that delete should be called on the
 * returned screw later and outside this class.
 * \return A pointer to a copy of the current screw.
 */
ScrewRectangularSectionClean* ScrewRectangularSectionClean::copy() const
{
    return new ScrewRectangularSectionClean(*this);
}

// ----------   SET FUNCTIONS   ----------
// Function to set the screw parameters
/*!
 * \details Function that takes parameters for all screw properties, checks if they are valid and then sets all screw
 * properties
 * \param[in] start The centre of the lower end of the screw (3D vector)
 * \param[in] length The length of the screw (Mdouble > 0)
 * \param[in] bladeRadius The outer screw blade radius (Mdouble > 0)
 * \param[in] shaftRadius The screw shaft radius (0 <= Mdouble < bladeRadius)
 * \param[in] numberOfTurns The number of turns the screw contains (Mdouble > 0)
 * \param[in] omega Rotation speed in rad per time unit (Mdouble)
 * \param[in] thickness The thickness of the screw blade (0 <= Mdouble < length)
 * \param[in] rightHandedness The right handedness of the screw, i.e. the direction of the screw-blade (bool)
 */
void ScrewRectangularSectionClean::set(const Vec3D& start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius,
                                       Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandedness)
{
    // check if the quantities are correct
    logger.assert_always(length > 0, "\nERROR: INVALID LENGTH.");
    logger.assert_always(bladeRadius > 0, "\nERROR: INVALID BLADE RADIUS.");
    logger.assert_always(shaftRadius < bladeRadius, "\nERROR: SHAFT RADIUS LARGER THAN BLADE RADIUS.");
    logger.assert_always(shaftRadius >= 0, "\nERROR: SHAFT RADIUS NEGATIVE");
    logger.assert_always(numberOfTurns > 0, "\nERROR: INVALID NUMBER OF SCREW TURNS.");
    logger.assert_always(thickness < length && thickness >= 0, "\nERROR: INVALID THICKNESS.");
    
    // user defined variables
    startPosition_ = start;
    length_ = length;
    rMax_ = bladeRadius;
    rMin_ = shaftRadius;
    numberOfTurns_ = numberOfTurns;
    omega_ = omega;
    thickness_ = thickness;
    rightHandedness_ = rightHandedness;
    setPrescribedAngularVelocity([omega](Mdouble time)
                                 { return Vec3D(0, 0, omega); });
    
    logger(INFO, "\n\nScrew parameters set.\n"
                   "User defined variables:\n"
                   "Start: %\n Length: %\n Blade radius: %\n Shaft radius: %\n, Number of turns: %\n"
                   "Angular velocity: %\n Thickness: %\n Initial angular offset: %\n Right handedness: %\n"
                   "Derived variables:\n Half thickness: %\n Pitch length: %",
           startPosition_, length_, rMax_, rMin_, numberOfTurns_, omega_, thickness_, getOrientation(),
           rightHandedness_, thickness_ / 2, length_ / numberOfTurns_);
}

/*!
 * \details Function to change the blade radius (outer radius) of the screw
 * \param[in] bladeRadius new outer radius (Mdouble > 0)
 */
void ScrewRectangularSectionClean::setRadius(Mdouble bladeRadius)
{
    logger.assert_always(bladeRadius > 0, "\nERROR: INVALID BLADE RADIUS.");
    rMax_ = bladeRadius;
}

/*!
 * Function to change the blade-thickness
 * \param[in] thickness new blade thickness (0 <= Mdouble < length)
 */
void ScrewRectangularSectionClean::setThickness(Mdouble thickness)
{
    logger.assert_always(thickness < length_ && thickness >= 0, "\nERROR: INVALID THICKNESS.");
    thickness_ = thickness;
}

/*!
 * Function to change the screw angular velocity
 * \param[in] omega  new angular velocity
 */
void ScrewRectangularSectionClean::setOmega(Mdouble omega)
{
    omega_ = omega;
    setPrescribedAngularVelocity([omega](Mdouble time)
                                 { return Vec3D(0, 0, omega); });
}


// ----------   GET FUNCTIONS   ------------------------------------------------------
//
//
/*!
 * \details The contact detection algorithm, where it is checked whether the given particle is in contact with the
 * screw, and if so, computes the distance between the centre of the particle and the screw and the normal of this
 * contact. Note, that the distance and normal are returned as output parameters, not as return values. Furthermore,
 * all the quantities in cylindrical coordinates are referred to the screw axis
 * \param[in] p particle for which is checked if there is an interaction with this screw (BaseParticle)
 * \param[out] distance if there is an interaction, will contain the distance between the particle and screw (Mdouble)
 * \param[out] normal_return if there is an interaction, will contain the normal of the contact between the particle and
 * the screw (Vec3D)
 * @return Boolean to indicate whether there is an interaction: true if there is an interaction, false if there is no
 * interaction
 */
bool
ScrewRectangularSectionClean::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial position of the particle
    const Vec3D distancePToStart = p.getPosition() - startPosition_;
    Mdouble rho2 = distancePToStart.X * distancePToStart.X + distancePToStart.Y * distancePToStart.Y;
    
    // if the particle is outside the cylinder containing the screw there is no collision
    if (rho2 > (rMax_ + p.getWallInteractionRadius(this)) * (rMax_ + p.getWallInteractionRadius(this)) ||
        distancePToStart.Z > length_ + p.getWallInteractionRadius(this) ||
        distancePToStart.Z < -p.getWallInteractionRadius(this))
    {
        return false;
    }
    Vec3D normalVector;
    Vec3D radialVector;
    Mdouble deltaN;
    computeNormalRadialDeltaN(p, normalVector, radialVector, deltaN);
    
    // radial position of the particle
    const Mdouble rho = std::sqrt(rho2);
    // collision check
    if (rho >= rMax_) // check for collision with the blade edge
    {
        // radial component of the distance between the particle centre and the blade edge
        const Mdouble deltaR = rho - rMax_;
        
        if (deltaN <= 0.0) // region 1
        {
            logger(DEBUG, "region 1");

            // the distance between the collision point and the particle centre
            distance = deltaR;
            
            // if the distance is negative prints an error message
            logger.assert_debug(distance >= 0, "\nCOLLISION ERROR WITH THE SCREW SIDE: OVERLAP > PARTICLE RADIUS, ZONE 1");
            
            // the collision has no other component other than the radial vector
            normal_return = -radialVector;

            return true;
        }
        else // region 2
        {
            logger(DEBUG, "region 2");

            // distance between the particle centre and the and the edge
            distance = std::sqrt(deltaN * deltaN + deltaR * deltaR);
            
            // if the particle-blade_edge distance is higher than the particle radius there is no collision
            if (distance > p.getWallInteractionRadius(this))
                return false;
            
            // if the distance is negative prints an error message
            logger.assert_debug(distance >= 0,
                          "\nCOLLISION ERROR WITH THE SCREW EDGE: OVERLAP > PARTICLE RADIUS, ZONE 2");

            // the normal to the edge through the particle centre
            // IMPORTANT: it needs to be normalized!
            normal_return = deltaN * normalVector - deltaR * radialVector;
            normal_return.normalise();

            return true;
        }
    }
    else if (rho >= rMin_ + p.getWallInteractionRadius(this)) // check for collision with the blade side only, region 3
    {
        logger(DEBUG, "region 3");

        // if the particle-blade_surface distance is higher than the collision threshold there is no collision
        if (deltaN > p.getWallInteractionRadius(this))
            return false;
        
        // region 3
        // the distance between the collision point and the particle centre
        distance = deltaN;
        
        // if the distance is negative prints an error message
        logger.assert_debug(distance >= 0.0 || getHandler()->getDPMBase()->getTime() < 1e-2,
                      "\nCOLLISION ERROR WITH THE SCREW SURFACE: OVERLAP > PARTICLE RADIUS, particle ID = %, ZONE 3",
                      p.getId());
        
        // the collision has no other component other than the normal to the blade
        normal_return = normalVector;
        return true;
    }
    else if (rho >= rMin_) // collision with the shaft
    {
        // radial distance between the particle centre and the shaft
        const Mdouble deltaR = rho - rMin_;
        
        if (deltaN > p.getWallInteractionRadius(this)) // region 4
        {
            logger(DEBUG, "region 4");
            // the distance between the collision point and the particle centre
            distance = deltaR;
            
            // if the distance is negative prints an error message
            logger.assert_debug(distance >= 0, "\nCOLLISION ERROR WITH THE SCREW SHAFT: OVERLAP > PARTICLE RADIUS, ZONE 4");
            
            // the collision has no other component other than the radial vector
            normal_return = -radialVector;
            
            return true;
        }
        else // region 5
        {
            logger(DEBUG, "region 5");
            // distance between the particle centre and the and the point (rMin + rp, delta + rp)
            distance = std::sqrt((deltaN - p.getWallInteractionRadius(this)) * (deltaN - p.getWallInteractionRadius(this)) +
                                 (deltaR - p.getWallInteractionRadius(this)) * (deltaR - p.getWallInteractionRadius(this)));
            
            // transform the former to be the actual distance of the particle from the contact point
            // on the rounded corner
            distance = p.getWallInteractionRadius(this) - distance;
            
            // if the distance is negative prints an error message
            logger.assert_debug(distance >= 0.0 || getHandler()->getDPMBase()->getTime() < 1e-2,
                          "\nCOLLISION ERROR WITH THE SCREW CORNER: OVERLAP > PARTICLE RADIUS, particle ID = %, ZONE 5",
                          p.getId());
            
            // the normal to the edge through the particle centre
            // IMPORTANT: it needs to be normalized!
            normal_return = fabs(deltaN - p.getWallInteractionRadius(this)) * normalVector -
                            fabs(deltaR - p.getWallInteractionRadius(this)) * radialVector;
            normal_return.normalise();
            
            return true;
        }
    }
    else
    {
        // if the distance from the centre is smaller than the shaft radius prints an error
        logger(ERROR, "\nERROR: rho < rMin, particle centre inside the shaft");
        return false;
    }
}

/*!
 * Computes the normal and radial vector for the screw at the position of the particle. Furthermore, it also computes
 *  deltaN, which is component of the particle-blade centre distance normal to the blade surface
 * \param[in] p particle whose position is used to compute the vectors and distance to the blade
 * \param[out] normalVector the vector in (x,y,z)-direction that points perpendicular to the screw blade
 * \param[out] radialVector the vector in (x,y)-direction that points outward from shaft to particle
 * \param[out] deltaN component of the particle-blade_center distance normal to the blade surface
 */
void ScrewRectangularSectionClean::computeNormalRadialDeltaN(const BaseParticle& p, Vec3D& normalVector,
                                                             Vec3D& radialVector, Mdouble& deltaN) const
{
    Vec3D distanceParticleStart = p.getPosition() - startPosition_;
    // squared radial position of the particle
    const Mdouble rho2 = distanceParticleStart.X * distanceParticleStart.X
                         + distanceParticleStart.Y * distanceParticleStart.Y;

    // radial position of the particle
    const Mdouble rho = std::sqrt(rho2);

    // The rescaled length of the screw (length_/(2*pi*numberOfTurns_)).
    const Mdouble h = 0.5 * length_ / (constants::pi * numberOfTurns_);

    // The pitch length of the screw (length_/numberOfTurns_).
    const Mdouble pitch = length_ / numberOfTurns_;
    const Mdouble deltaZ = computeDeltaZ(distanceParticleStart, h, pitch);


    // trigonometric functions relative to the particle angle
    const Mdouble cosXi = (p.getPosition().X - startPosition_.X) / rho;
    const Mdouble sinXi = (p.getPosition().Y - startPosition_.Y) / rho;
    if (rightHandedness_)
    {
        normalVector.X = h * sinXi;
        normalVector.Y = -h * cosXi;
        normalVector.Z = rho;
    }
    else
    {
        normalVector.X = -h * cosXi;
        normalVector.Y = h * sinXi;
        normalVector.Z = rho;
    }

    // takes the right direction (+/-) of the vector and normalizes
    normalVector *= -deltaZ;
    normalVector /= normalVector.getLength();

    // radial vector at the particle position
    radialVector.X = cosXi;
    radialVector.Y = sinXi;
    radialVector.Z = 0.0;

    // The half-thickness of the screw.
    const Mdouble delta = 0.5 * thickness_;

    // cosine of the helix angle at the particle position
    const Mdouble cosEta = 1.0 / std::sqrt(1.0 + (h * h / rho / rho));

    // component of the particle-blade_center distance normal to the blade surface
    deltaN = fabs(deltaZ) * cosEta - delta;
}

/*!
 * Auxiliary function that computes the oriented axial distance between the particle's centre and the blade centre
 * \param[in] distanceParticleStart distance between the starting point of the screw, start_, and the particle position
 * (Vec3D)
 * \param[in] h The rescaled length of the screw (length_/(2*pi*numberOfTurns_)) (Mdouble > 0)
 * \param[in] pitch The pitch length of the screw (length_/numberOfTurns_) (Mdouble > 0)
 * @return oriented axial distance between the particle's centre and the blade centre
 */
Mdouble ScrewRectangularSectionClean::computeDeltaZ(const Vec3D& distanceParticleStart, Mdouble h, Mdouble pitch) const
{
    // oriented axial distance between the particle's centre and the blade centre
    Mdouble deltaZ;

    // normal to the blade at the particle position
    if (rightHandedness_) // right-handed thread
    {
        // angular coordinate of the particle
        // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
        Mdouble xi = atan2(distanceParticleStart.Y, distanceParticleStart.X);
        if (xi < 0.0)
            xi += 2.0 * constants::pi;

        deltaZ = fmod(distanceParticleStart.Z - h * (xi + getOrientation().getAxis().Z) -
                      static_cast<int> (distanceParticleStart.Z / pitch), pitch);
        logger(DEBUG, "xi: %, deltaZ: %", xi, deltaZ);
    }
    else // left-handed thread
    {
        // angular coordinate of the particle
        // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
        Mdouble xi = atan2(-distanceParticleStart.Y, distanceParticleStart.X);
        if (xi < 0.0)
            xi += 2.0 * constants::pi;
        xi += 0.5 * constants::pi;
        xi = fmod(xi, 2.0 * constants::pi);

        deltaZ = fmod(distanceParticleStart.Z - h * (xi + 0.5 * constants::pi - getOrientation().getAxis().Z) -
                      static_cast<int> (distanceParticleStart.Z / pitch), pitch);
        logger(DEBUG, "xi: %, deltaZ: %", xi, deltaZ);
    }

    if (deltaZ > 0.5 * pitch)
    {
        deltaZ -= pitch;
    }
    else if (deltaZ < -0.5 * pitch)
    {
        deltaZ += pitch;
    }
    return deltaZ;
}

/*!
 * Function that returns the current angular offset
 * \return current angular offset (Mdouble)
 */
Mdouble ScrewRectangularSectionClean::getOffset() const
{
    return getOrientation().getAxis().Z;
}

//
/*!
 * Function that returns the screw's righthandedness
 * \return  rightHandedness_ (boolean)
 */
bool ScrewRectangularSectionClean::getRightHandedness() const
{
    return rightHandedness_;
}

// ----------   ROTATION FUNCTIONS   ------------------------------------------------------
/*!
 * Function to make the screw rotate by adding \omega*dt to the angular offset
 * \param[in] dt time duration for which the screw should be turned (Mdouble)
 */
void ScrewRectangularSectionClean::rotate(Mdouble dt)
{
    incrementOffset(dt*omega_);
}

//
/*!
 * Function to make the screw rotate by incrementing the angular offset
 * \param[in] off The amount by which the screw should be rotated (Mdouble)
 */
void ScrewRectangularSectionClean::incrementOffset(Mdouble off)
{
    setOrientation(getOrientation() + Quaternion(Vec3D(0, 0, off)));
}

// ------------INPUT/OUTPUT   -------------------------------------------------------------
/*!
 * Read the screw properties from an input-stream, for example a restart file
 * \param[in|out] is inputstream from which the screw is read (std::istream)
 */
void ScrewRectangularSectionClean::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> startPosition_
       >> dummy >> length_
       >> dummy >> rMax_
       >> dummy >> rMin_
       >> dummy >> numberOfTurns_
       >> dummy >> omega_
       >> dummy >> thickness_
       >> dummy >> rightHandedness_;
}

/*!
* Read the screw properties from an input-stream, for example a restart file
* \param[in|out] is input-stream from which the screw is read (std::istream)
*/
void ScrewRectangularSectionClean::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> startPosition_
       >> dummy >> length_
       >> dummy >> rMax_
       >> dummy >> rMin_
       >> dummy >> numberOfTurns_
       >> dummy >> omega_
       >> dummy >> rightHandedness_;
}

/*!
 * Write the screw properties to an output-stream, for example a restart file, can also be used to write to the console
 * by using std::cout as a parameter.
 * \param[in|out] os output-stream to which the screw properties are written (std::ostream)
 */
void ScrewRectangularSectionClean::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << startPosition_
       << " Length " << length_
       << " Blade radius " << rMax_
       << " Shaft radius " << rMin_
       << " Turns " << numberOfTurns_
       << " Omega " << omega_
       << " Thickness " << thickness_
       << " Right handedness " << rightHandedness_;
}

/*!
 * Returns the string "ScrewRectangularSectionClean"
 * @return "ScrewRectangularSectionClean" (std::string)
 */
std::string ScrewRectangularSectionClean::getName() const
{
    return "ScrewRectangularSectionClean";
}

// Checks for the interaction between a particle p at a time timeStamp
// In case of interaction returns a pointer to the BaseInteraction happened between the Screw and the BaseParticle at time timeStamp
BaseInteraction* ScrewRectangularSectionClean::getInteractionWith(BaseParticle* p, unsigned timeStamp,
                                                                               InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    if (getDistanceAndNormal(*p, distance, normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        return c;
    }
    return nullptr;
}
