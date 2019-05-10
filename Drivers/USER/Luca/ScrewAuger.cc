#include "ScrewAuger.h"
#include "array"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

///\todo implement the finite screw
///\todo implement the right/left handed screw
///\todo check the implementation for the left handed screw, especially the angle and deltaZ computation
///\todo add multi-threads option
///\todo add proper Mercury headers
///\todo rewrite the details above the .h
///\todo make this compliant with Trunk

/*
 Last update 07.05.19

 ----- BLADE SECTION -----
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

/*!
 * \details Make a Screw which is centered in the origin, has a length of 1, one
 * revolution, a blade radius of 1, a shaft radius of 0.5, turns with 1 revolution per second,
 * has 0 thickness, starts at its normal initial point and is right-handed.
 */
ScrewAuger::ScrewAuger()
{
    // explicit variables
    screwOrigin_.setZero();
    length_ = 1.0;
    bladeRadius_ = 1.0;
    shaftRadius_ = 0.5;
    numberOfTurns_ = 1.0;
    omega_ = 1.0;
    thickness_ = 0.0;
    rightHandedness_ = true;
    setOrientationViaNormal({0, 0, 1});

    // implicit variables
    angularOffset_ = 0.0;
    delta_ = 0.0;
    pitch_ = 1.0;
    rescaledPitch_ = 0.5/constants::pi;

    logger(DEBUG, "ScrewAuger() constructor finished.");
}

/*!
 * \details Copy constructor: copy all values of the given screw into the constructed screw
 * \param[in] other screw that should be copied
 */
ScrewAuger::ScrewAuger(const ScrewAuger& other)
    : BaseWall(other)
{
    // explicit variables
    screwOrigin_ = other.screwOrigin_;
    length_ = other.length_;
    bladeRadius_ = other.bladeRadius_;
    shaftRadius_ = other.shaftRadius_;
    numberOfTurns_ = other.numberOfTurns_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    angularOffset_ = other.angularOffset_;
    rightHandedness_ = other.rightHandedness_;

    // implicit variables
    delta_ = 0.5*other.thickness_;
    pitch_ = other.length_/other.numberOfTurns_;
    rescaledPitch_ = 0.5*other.length_/(constants::pi * other.numberOfTurns_);

    logger(DEBUG, "ScrewAuger(const ScrewAuger&) copy constructor finished.");
}

/*!
 * \details Constructor that takes parameters for all screw properties, checks if they are valid and then sets all screw properties
 * \param[in] origin The centre of the lower end of the screw (3D vector)
 * \param[in] length The length of the screw (Mdouble > 0)
 * \param[in] bladeRadius The screw blade radius (Mdouble > 0)
 * \param[in] shaftRadius The screw shaft radius (0 <= Mdouble < bladeRadius)
 * \param[in] numberOfTurns The number of turns the screw contains (Mdouble > 0)
 * \param[in] omega Screw rotational velocity around its axix in rad/s (Mdouble)
 * \param[in] thickness The thickness of the screw blade (0 <= Mdouble < length)
 * \param[in] rightHandedness The right handedness of the screw, i.e. the direction of the screw-blade (bool)
 */
ScrewAuger::ScrewAuger(Vec3D origin, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandedness)
{
    // check if the quantities are correct
    logger.assert_always(length > 0, "\nERROR: INVALID LENGTH.");
    logger.assert_always(bladeRadius > 0, "\nERROR: INVALID BLADE RADIUS.");
    logger.assert_always(shaftRadius < bladeRadius, "\nERROR: SHAFT RADIUS LARGER THAN BLADE RADIUS.");
    logger.assert_always(shaftRadius > 0, "\nERROR: SHAFT RADIUS NEGATIVE");
    logger.assert_always(numberOfTurns > 0, "\nERROR: INVALID NUMBER OF SCREW TURNS.");
    logger.assert_always(thickness < length && thickness >= 0, "\nERROR: INVALID THICKNESS.");

    // explicit variables
    screwOrigin_ = origin;
    length_ = length;
    bladeRadius_ = bladeRadius;
    shaftRadius_ = shaftRadius;
    numberOfTurns_ = numberOfTurns;
    omega_ = omega;
    thickness_ = thickness;
    rightHandedness_ = rightHandedness;

    // implicit variables
    angularOffset_ = 0.0;
    delta_ = 0.5*thickness_;
    pitch_ = length_/numberOfTurns_;
    rescaledPitch_ = 0.5 * pitch_/constants::pi;

    logger(INFO, "\n\nScrew parameters set.\n"
           "User defined variables:\n"
           "Origin: %\n Length: %\n Blade radius: %\n Shaft radius: %\n, Number of turns: %\n"
           "Angular velocity: %\n Thickness: %\n Right handedness: %\n"
           "Derived variables:\n"
           "Initial angular offset: %\n Half thickness: %\n Pitch length: %\n Rescaled pitch length : %",
           screwOrigin_, length_, bladeRadius_, shaftRadius_, numberOfTurns_, omega_, thickness_, rightHandedness_,
           angularOffset_, delta_, pitch_, rescaledPitch_);

    logger(DEBUG, "ScrewAuger(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, bool) constructor finished.");
}

/*!
 * \details Destructor, just prints a message if the logger is set to debug mode
 */
ScrewAuger::~ScrewAuger()
{
    logger(DEBUG, "~ScrewAuger() finished, destroyed the screw.");
}

/*!
 * \details Function that makes a copy of the current screw, by calling new. Note that delete should be called on the
 * returned screw later and outside this class.
 * \return A pointer to a copy of the current screw.
 */
ScrewAuger* ScrewAuger::copy() const
{
    return new ScrewAuger(*this);
}

// --------------------   SET FUNCTIONS   --------------------
// Function to set the screw parameters
/*!
 * \details Function that takes parameters for all screw properties, checks if they are valid and then sets all screw properties
 * \param[in] origin The centre of the lower end of the screw (3D vector)
 * \param[in] length The length of the screw (Mdouble > 0)
 * \param[in] bladeRadius The screw blade radius (Mdouble > 0)
 * \param[in] shaftRadius The screw shaft radius (0 <= Mdouble < bladeRadius)
 * \param[in] numberOfTurns The number of turns the screw contains (Mdouble > 0)
 * \param[in] omega Screw rotational velocity around its axix in rad/s (Mdouble)
 * \param[in] thickness The thickness of the screw blade (0 <= Mdouble < length)
 * \param[in] rightHandedness The right handedness of the screw, i.e. the direction of the screw-blade (bool)
 */
void ScrewAuger::set(Vec3D origin, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandedness)
{
    // check if the quantities are correct
    logger.assert_always(length > 0, "\nERROR: INVALID LENGTH.");
    logger.assert_always(bladeRadius > 0, "\nERROR: INVALID BLADE RADIUS.");
    logger.assert_always(shaftRadius < bladeRadius, "\nERROR: SHAFT RADIUS LARGER THAN BLADE RADIUS.");
    logger.assert_always(shaftRadius > 0, "\nERROR: SHAFT RADIUS NEGATIVE");
    logger.assert_always(numberOfTurns > 0, "\nERROR: INVALID NUMBER OF SCREW TURNS.");
    logger.assert_always(thickness < length && thickness >= 0, "\nERROR: INVALID THICKNESS.");

    // explicit variables
    screwOrigin_ = origin;
    length_ = length;
    bladeRadius_ = bladeRadius;
    shaftRadius_ = shaftRadius;
    numberOfTurns_ = numberOfTurns;
    omega_ = omega;
    thickness_ = thickness;
    rightHandedness_ = rightHandedness;

    // implicit variables
    angularOffset_ = 0.0;
    delta_ = 0.5 * thickness_;
    pitch_ = length_/numberOfTurns_;
    rescaledPitch_ = 0.5 * pitch_/constants::pi;

    logger(INFO, "\n\nScrew parameters set.\n"
           "User defined variables:\n"
           "Origin: %\n Length: %\n Blade radius: %\n Shaft radius: %\n, Number of turns: %\n"
           "Angular velocity: %\n Thickness: %\n Right handedness: %\n"
           "Derived variables:\n"
           "Initial angular offset: %\n Half thickness: %\n Pitch length: %\n Rescaled pitch length : %",
           screwOrigin_, length_, bladeRadius_, shaftRadius_, numberOfTurns_, omega_, thickness_, rightHandedness_,
           angularOffset_, delta_, pitch_, rescaledPitch_);
}

/*!
 * \details Function to change the length of the screw
 * \param[in] length new screw length (Mdouble > 0 && Mdouble > thickness)
 */
void ScrewAuger::setLength(Mdouble length)
{
    logger.assert_always(length > 0 && length > thickness_, "\nERROR: INVALID LENGTH.");
    length_ = length;
}

/*!
 * \details Function to change the blade radius of the screw
 * \param[in] bladeRadius new blade radius (Mdouble > 0 && Mdouble > shaftRadius)
 */
void ScrewAuger::setBladeRadius(Mdouble bladeRadius)
{
    logger.assert_always(bladeRadius > 0 && bladeRadius > shaftRadius_, "\nERROR: INVALID BLADE RADIUS.");
    bladeRadius_ = bladeRadius;
}

/*!
 * \details Function to change the shaft radius of the screw
 * \param[in] shaftRadius new shaft radius (Mdouble > 0 && Mdouble < bladeRadius)
 */
void ScrewAuger::setShaftRadius(Mdouble shaftRadius)
{
    logger.assert_always(shaftRadius > 0 && shaftRadius < bladeRadius_, "\nERROR: INVALID SHAFT RADIUS.");
    shaftRadius_ = shaftRadius;
}

/*!
 * \details Function to change the blade thickness of the screw
 * \param[in] thickness new blade thickness (0 < Mdouble < length)
 */
void ScrewAuger::setThickness(Mdouble thickness)
{
    logger.assert_always(thickness > 0 && thickness < length_, "\nERROR: INVALID BLADE THICKNESS.");
    thickness_ = thickness;
}

/*!
 * \details Function to change the screw angular velocity
 * \param[in] omega new angular velocity
 */
void ScrewAuger::setOmega(Mdouble omega)
{
    omega_ = omega;
    setPrescribedAngularVelocity([omega](Mdouble time){ return Vec3D(0, 0, omega); });
}

/*!
 * \details Function to change the screw angle
 * \param[in] angle new angular offset
 */
void ScrewAuger::setOffset(Mdouble angularOffset)
{
    angularOffset_ = angularOffset;
}

// --------------------   GET FUNCTIONS   --------------------
// Function to get the screw parameters or related quantities
/*!
 * \details Function that returns the particle coordinate in cylindrical frame of reference, with respect to the screw origin
 * \param[in] particle position in euclidean coordinates with respect to the screw (Vec3D)
 * \return particle position in cylindrical coordinates with respect to the screw (Vec3D)
 */
Vec3D ScrewAuger::getCylindricalCoordinate(const Vec3D pInScrewReference) const
{
    Vec3D pInCylindricalFrame;
    pInCylindricalFrame.X = std::sqrt(pInScrewReference.Y * pInScrewReference.Y + pInScrewReference.Z * pInScrewReference.Z);
    pInCylindricalFrame.Y = atan2(pInScrewReference.Z, pInScrewReference.Y);
    if (pInCylindricalFrame.Y < 0.0) pInCylindricalFrame.Y += 2.0 * constants::pi;
    pInCylindricalFrame.Z = pInScrewReference.X;

    return pInCylindricalFrame;
}

/*!
 * \details Function that returns the oriented axial distance (can be negative!) between the particle and the blade centre
 * \param[in] particle position in cylindrical coordinates with respect to the screw (Vec3D)
 * \return oriented axial distance between particle centre and blade centre (Mdouble)
 */
Mdouble ScrewAuger::getDeltaX(const Vec3D pInCylindricalFrame) const
{
    Mdouble deltaX;

    if (rightHandedness_) // right-handed thread
    {
        deltaX = fmod(pInCylindricalFrame.Z - rescaledPitch_*(pInCylindricalFrame.Y + 2.0*constants::pi*static_cast<int>(pInCylindricalFrame.Z/pitch_)) - rescaledPitch_*angularOffset_, pitch_);
    }
    else // left-handed thread
    {
        deltaX = fmod(pInCylindricalFrame.Z - rescaledPitch_*(pInCylindricalFrame.Y + 2.0*constants::pi*static_cast<int>(pInCylindricalFrame.Z/pitch_)) - rescaledPitch_*angularOffset_, pitch_);
    }

    if (deltaX > 0.5*pitch_)
    {
        deltaX -= pitch_;
    }
    else if (deltaX < -0.5*pitch_)
    {
        deltaX += pitch_;
    }

    return deltaX;
}

/*!
 * \details Function that returns the normal vector with respect to the blade surface at the particle's position
 * \param[in] particle position in euclidean coordinates with respect to the screw (Vec3D)
 * \param[in] particle position in cylindrical coordinates with respect to the screw (Vec3D)
 * \return vector through the particle's centre perpendicular to the screw surface (Vec3D)
 */
Vec3D ScrewAuger::getNormalVector(const Vec3D pInScrewReference, const Vec3D pInCylindricalFrame) const
{
    Vec3D normalVector;

    if (rightHandedness_) // right-handed thread
    {
        normalVector.X = pInCylindricalFrame.X;
        normalVector.Y = rescaledPitch_ * pInScrewReference.Z/pInCylindricalFrame.X;
        normalVector.Z = - rescaledPitch_ * pInScrewReference.Y/pInCylindricalFrame.X;
    }
    else // left-handed thread
    {
        normalVector.X = pInCylindricalFrame.X;
        normalVector.Y = - rescaledPitch_ * pInScrewReference.Y/pInCylindricalFrame.X;
        normalVector.Z = rescaledPitch_ * pInScrewReference.Z/pInCylindricalFrame.X;
    }

    // takes the right direction of the vector
    normalVector *= - getDeltaX(pInCylindricalFrame);
    normalVector /= normalVector.getLength();

    return normalVector;
}

/*!
 * \details Function that returns the radial vector with respect to the blade axis at the particle's position
 * \param[in] particle position in euclidean coordinates with respect to the screw (Vec3D)
 * \param[in] particle position in cylindrical coordinates with respect to the screw (Vec3D)
 * \return vector through the particle's centre perpendicular to the screw axis (Vec3D)
 */
Vec3D ScrewAuger::getRadialVector(const Vec3D pInScrewReference, const Vec3D pInCylindricalFrame) const
{
    Vec3D radialVector;

    radialVector.X = 0.0;
    radialVector.Y = - pInScrewReference.Y/pInCylindricalFrame.X;
    radialVector.Z = - pInScrewReference.Z/pInCylindricalFrame.X;

    return radialVector;
}

/*!
 * \details The contact detection algorithm, where it is checked whether the given particle is in contact with the
 * screw, and if so, computes the distance between the centre of the particle and the screw and the normal of this
 * contact. Note, that the distance and normal are returned as output parameters, not as return values. Furthermore,
 * all the quantities in cylindrical coordinates are referred to the screw axis
 * \param[in] p particle for which is checked if there is an interaction with this screw (BaseParticle)
 * \param[out] distance if there is an interaction, will contain the distance between the particle and screw (Mdouble)
 * \param[out] normal_return if there is an interaction, will contain the normal of the contact between the particle and
 * the screw (Vec3D)
 * \return Boolean to indicate whether there is an interaction: true if there is an interaction, false if there is no
 * interaction (bool)
 */
bool ScrewAuger::getDistanceAndNormalLabCoordinates(Vec3D position, Mdouble wallInteractionRadius, Mdouble& distance,
                                                    Vec3D& normal_return) const
{
    // position of the particle with respect to the screw origin
    const Vec3D pInScrewReference = position - screwOrigin_;

    // squared radial distance
    Mdouble rho2 = pInScrewReference.Y * pInScrewReference.Y + pInScrewReference.Z * pInScrewReference.Z;

    // if the particle is outside the cylinder containing the screw there is no collision
    if (rho2 > (bladeRadius_ + wallInteractionRadius) * (bladeRadius_ + wallInteractionRadius) ||
        pInScrewReference.X > length_ + wallInteractionRadius ||
        pInScrewReference.X < - wallInteractionRadius)
    {
        return false;
    }

    // position of the particle with respect to the screw origin in cylindrical coordinates
    const Vec3D pInCylindricalScrewReference = getCylindricalCoordinate(pInScrewReference);

    // component of the particle-blade_center distance normal to the blade surface
    const Mdouble deltaN = fabs(getDeltaX(pInCylindricalScrewReference)) * pInCylindricalScrewReference.X / std::sqrt(pInCylindricalScrewReference.X * pInCylindricalScrewReference.X + rescaledPitch_ * rescaledPitch_) - delta_;

    // collision check
    if (pInCylindricalScrewReference.X >= bladeRadius_) // check for collision with the blade edge
    {
        // radial component of the distance between the particle centre and the blade edge
        const Mdouble deltaR = pInCylindricalScrewReference.X - bladeRadius_;

        if (deltaN <= 0.0) // region 1
        {
            logger(DEBUG, "region 1");

            // distance between the particle centre and the blade flat side
            distance = deltaR;

            // if the distance is negative prints an error message
            logger.assert(distance >= 0.0, "\nCOLLISION ERROR WITH THE SCREW SIDE: OVERLAP > PARTICLE RADIUS, REGION 1");

            // the collision has no other component other than the radial vector
            normal_return = getRadialVector(pInScrewReference, pInCylindricalScrewReference);

            return true;
        }
        else // region 2
        {
            logger(DEBUG, "region 2");

            // distance between the particle centre and the edge
            distance = std::sqrt(deltaN * deltaN + deltaR * deltaR);

            // if the particle-blade_edge distance is higher than the particle radius there is no collision
            if (distance > wallInteractionRadius)
                return false;

            // if the distance is negative prints an error message
            logger.assert(distance >= 0.0, "\nCOLLISION ERROR WITH THE SCREW EDGE: OVERLAP > PARTICLE RADIUS, REGION 2");

            // the collision has mixed radial and normal components
            normal_return = deltaN * getNormalVector(pInScrewReference, pInCylindricalScrewReference) + deltaR * getRadialVector(pInScrewReference, pInCylindricalScrewReference);
            normal_return /= normal_return.getLength();

            return true;
        }
    }
    else if (pInCylindricalScrewReference.X >= shaftRadius_ + wallInteractionRadius) // region 3
    {
        logger(DEBUG, "region 3");

        // if the particle-blade_surface distance is higher than the collision threshold there is no collision
        if (deltaN > wallInteractionRadius)
            return false;

        // distance between the particle centre and the blade surface
        distance = deltaN;

        // if the distance is negative prints an error message
        logger.assert(distance >= 0.0, "\nCOLLISION ERROR WITH THE SCREW BLADE SURFACE: OVERLAP > PARTICLE RADIUS, REGION 3");

        // the collision has no other component other than the normal to the blade
        normal_return = getNormalVector(pInScrewReference, pInCylindricalScrewReference);

        return true;
    }
    else if (pInCylindricalScrewReference.X >= shaftRadius_) // check for collision with the shaft
    {
        // radial component of the distance between the particle centre and the screw shaft
        const Mdouble deltaR = pInCylindricalScrewReference.X - shaftRadius_;

        if (deltaN > wallInteractionRadius) // region 4
        {
            logger(DEBUG, "region 4");

            // distance between the particle centre and the shaft
            distance = deltaR;

            // if the distance is negative prints an error message
            logger.assert(distance >= 0.0, "\nCOLLISION ERROR WITH THE SCREW SHAFT: OVERLAP > PARTICLE RADIUS, REGION 4");

            // the collision has no other component other than the radial vector
            normal_return = getRadialVector(pInScrewReference, pInCylindricalScrewReference);

            return true;
        }
        else // region 5
        {
            // distance between the particle centre and the and the point (shaftRadius + rp, delta + rp)
            distance = std::sqrt((deltaN - wallInteractionRadius) * (deltaN - wallInteractionRadius) + (deltaR - wallInteractionRadius) * (deltaR - wallInteractionRadius));

            // transform the former distance to be the actual distance of the particle from the contact point on the rounded corner
            distance = wallInteractionRadius - distance;

            // if the distance is negative prints an error message
            logger.assert(distance >= 0.0, "\nCOLLISION ERROR WITH THE SCREW CORNER: OVERLAP > PARTICLE RADIUS, REGION 5");

            // the collision has mixed radial and normal components
            normal_return = fabs(deltaN - wallInteractionRadius) * getNormalVector(pInScrewReference, pInCylindricalScrewReference) + fabs(deltaR - wallInteractionRadius) * getRadialVector(pInScrewReference, pInCylindricalScrewReference);
            normal_return /= normal_return.getLength();

            return true;
        }
    }
    else
    {
        // if the distance from the centre is smaller than the shaft radius prints an error
        logger(ERROR, "\nERROR: AXIAL DISTANCE < SHAFT RADIUS: PARTICLE CENTRE INSIDE THE SCREW SHAFT");
        return false;
    }
}

/*!
 * \details Function that transforms the coordinates into the position-orientation frame.
 * It also takes into account the mixed species of the particle-wall system, therefore using the correct interaction distance automatically.
 * \return Boolean to indicate whether there is an interaction: true if there is an interaction, false if there is no
 * interaction (bool)
 */
bool ScrewAuger::getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const
{
    // transform coordinates into position-orientation frame
    Vec3D position = P.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    
    // automatically uses the correct interaction distance obtained by the mixed species
    BaseSpecies* s = getHandler()->getDPMBase()->speciesHandler.getMixedObject(P.getSpecies(), getSpecies());
    
    // if there is a collision in the lab frame, then switches back to the normal frame
    if (getDistanceAndNormalLabCoordinates(position, P.getRadius() + s->getInteractionDistance(), distance,
                                           normal_return))
    {
        getOrientation().rotate(normal_return);
        return true;
    }
    else
    {
        return false;
    }
}

/*!
 * \details Function that returns the current angular offset
 * \return current angular offset (Mdouble)
 */
Mdouble ScrewAuger::getAngularOffset() const
{
    return angularOffset_;
}

/*!
 * \details Function that returns the current screw angular velocity
 * \return current angular velocity (Mdouble)
 */
Mdouble ScrewAuger::getOmega() const
{
    return omega_;
}

/*!
 * \details Function that returns the screw righthandedness
 * \return Boolean that is true if the screw is right-handed (bool)
 */
bool ScrewAuger::getRightHandedness() const
{
    return rightHandedness_;
}

// --------------------   ROTATION FUNCTIONS   --------------------
// Function to rotate the screw
/*!
 * \details Function to make the screw rotate by incrementing the angular offset
 * \param[in] offset the amount by which the screw should be rotated (Mdouble)
 */
void ScrewAuger::incrementOffset(Mdouble angularOffsetIncrement)
{
    angularOffset_ += angularOffsetIncrement;
}

/*!
 * \details Function to make the screw rotate by incrementing by dt*omega_the angular offset
 * \param[in] dt time duration for which the screw should be turned (Mdouble)
 */
void ScrewAuger::rotate(Mdouble dt)
{
    BaseInteractable::rotate({dt*omega_, 0.0, 0.0});
    incrementOffset(dt*omega_);
}

// --------------------   INPUT/OUTPUT FUNCTIONS   --------------------
// Function designated for input/output purpouses
/*!
 * \details Read the screw properties from an input-stream, for example a restart file
 * \param[in|out] is inputstream from which the screw is read (std::istream)
 */
void ScrewAuger::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> screwOrigin_
        >> dummy >> length_
        >> dummy >> bladeRadius_
        >> dummy >> shaftRadius_
        >> dummy >> numberOfTurns_
        >> dummy >> omega_
        >> dummy >> thickness_
        >> dummy >> angularOffset_
        >> dummy >> rightHandedness_;
}

/*!
 * \details Read the screw properties from an input-stream, for example a restart file
 * \param[in|out] is input-stream from which the screw is read (std::istream)
 */
void ScrewAuger::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> screwOrigin_
        >> dummy >> length_
        >> dummy >> bladeRadius_
        >> dummy >> shaftRadius_
        >> dummy >> numberOfTurns_
        >> dummy >> omega_
        >> dummy >> angularOffset_
        >> dummy >> rightHandedness_;
}

/*!
 * \details Write the screw properties to an output-stream, for example a restart file, and can also be used to write to the console by using std::cout as a parameter.
 * \param[in|out] os output-stream to which the screw properties are written (std::ostream)
 */
void ScrewAuger::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " Origin " << screwOrigin_
        << " Length " << length_
        << " Blade radius " << bladeRadius_
        << " Shaft radius " << shaftRadius_
        << " Turns " << numberOfTurns_
        << " Omega " << omega_
        << " Thickness " << thickness_
        << " Angular offset " << angularOffset_
        << " Right handeness " << rightHandedness_;
}

/*!
 * \details Returns the string "ScrewAuger"
 * \return the std string "ScrewAuger" (std::string)
 */
std::string ScrewAuger::getName() const
{
    return "ScrewAuger";
}

void ScrewAuger::writeVTK(VTKContainer& vtk) const
{
    //number of points in radial direction (for one side of the screw)
    unsigned nr = 5;
    //number of points in axial direction
    unsigned nz = static_cast<unsigned int>(99 * fabs(numberOfTurns_));
    
    unsigned long nPoints = vtk.points.size();
    vtk.points.reserve(nPoints + 2 * nr * nz);
    Vec3D contactPoint;
    // either one or two helices
    for (Mdouble offset = angularOffset_; offset<=angularOffset_; offset+=0.5)
    {
        for (unsigned iz = 0; iz < nz; iz++)
        {
            for (unsigned ir = 0; ir < nr; ir++)
            {
                double q = (double) iz / nz;
                double r = (double) ir / (nr - 1) * bladeRadius_;
                contactPoint.Y = screwOrigin_.Y - r * cos(2.0 * constants::pi * (offset + numberOfTurns_ * q));
                contactPoint.Z = screwOrigin_.Z - r * sin(2.0 * constants::pi * (offset + numberOfTurns_ * q));
                contactPoint.X = screwOrigin_.X + q * length_ - thickness_;
                getOrientation().rotate(contactPoint);
                contactPoint += getPosition();
                vtk.points.push_back(contactPoint);
            }
            for (unsigned ir = 0; ir < nr; ir++)
            {
                double q = (double) iz / nz;
                double r = (double) (nr - 1 - ir) / (nr - 1) * bladeRadius_;
                contactPoint.Y = screwOrigin_.Y - r * cos(2.0 * constants::pi * (offset + numberOfTurns_ * q));
                contactPoint.Z = screwOrigin_.Z - r * sin(2.0 * constants::pi * (offset + numberOfTurns_ * q));
                contactPoint.X = screwOrigin_.X + q * length_ + thickness_;
                getOrientation().rotate(contactPoint);
                contactPoint += getPosition();
                vtk.points.push_back(contactPoint);
            }
        }
    }
    
    unsigned long nCells = vtk.triangleStrips.size();
    vtk.triangleStrips.reserve(nCells + (nz - 1));
    for (unsigned iz = 0; iz < nz-1; iz++)
    {
        //skip step that would connect the two screw parts
        if (iz==nz-1) ++iz;
        std::vector<double> cell;
        cell.reserve(2 * nr);
        for (unsigned ir = 0; ir < 2*nr; ir++)
        {
            cell.push_back(nPoints + ir + 2*iz * nr);
            cell.push_back(nPoints + ir + 2*(iz + 1) * nr);
        }
        vtk.triangleStrips.push_back(cell);
    }
}

void ScrewAuger::writeVTK(std::string filename) const
{
    VTKContainer vtk;
    writeVTK(vtk);
    
    std::stringstream file;
    file << "# vtk DataFile Version 2.0\n"
    << getName() << "\n"
    "ASCII\n"
    "DATASET UNSTRUCTURED_GRID\n"
    "POINTS " << vtk.points.size() << " double\n";
    for (const auto& vertex : vtk.points)
        file << vertex << '\n';
    file << "\nCELLS " << vtk.triangleStrips.size() << ' ' << 4 * vtk.triangleStrips.size() << "\n";
    for (const auto& face : vtk.triangleStrips)
        file << "3 " << face[0] << ' ' << face[1] << ' ' << face[2] << '\n';
    file << "\nCELL_TYPES " << vtk.triangleStrips.size() << "\n";
    for (const auto& face : vtk.triangleStrips)
        file << "5\n";
    helpers::writeToFile(filename, file.str());
}
