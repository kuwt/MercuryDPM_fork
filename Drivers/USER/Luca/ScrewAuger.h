// This function defines an Archimedes' screw in the z-direction from a (constant) starting point,
// a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s).

// The screw H=H(r,theta) is a 2D surface embedded in E3, which parametric equations are:
// xH = r*Cos(theta)
// yH = r*Sin(theta)
// zH = L*theta/(2*pi*N)
// where r and theta are the parameters with r\in[shaftRadius;bladeRadius] and theta\in[0;2*pi*N[.

#ifndef SCREWAUGER_H
#define SCREWAUGER_H

#include "Walls/BaseWall.h"

class ScrewAuger : public BaseWall
{
public:

    ///\details Default constructor: set all variables to sensible defaults
    ScrewAuger();

    /*!
     * \details Copy constructor: copy all values of the given screw into the constructed screw
     * \param[in] other screw that should be copied
     */
    ScrewAuger(const ScrewAuger& other);

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
     ScrewAuger(Vec3D origin, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandeness);

     /*!
      * \details Destructor, just prints a message if the logger is set to debug mode
      */
    ~ScrewAuger();

    /*!
     * \details Function that makes a copy of the current screw, by calling new. Note that delete should be called on the
     * returned screw later and outside this class.
     * \return A pointer to a copy of the current screw.
     */
    ScrewAuger* copy() const final;

    // --------------------   SET FUNCTIONS   --------------------
    // Function to set the screw parameters
    // Sets the screw parameters
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
    void set(Vec3D origin, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandeness);

    /*!
     * \details Function to change the length of the screw
     * \param[in] length new screw length (Mdouble > 0 && Mdouble > thickness)
     */
    void setLength(Mdouble length);

    /*!
     * \details Function to change the blade radius of the screw
     * \param[in] bladeRadius new blade radius (Mdouble > 0 && Mdouble > shaftRadius)
     */
    void setBladeRadius(Mdouble bladeRadius);

    /*!
     * \details Function to change the shaft radius of the screw
     * \param[in] shaftRadius new shaft radius (Mdouble > 0 && Mdouble < bladeRadius)
     */
    void setShaftRadius(Mdouble shaftRadius);

    /*!
     * \details Function to change the blade thickness of the screw
     * \param[in] thickness new blade thickness (0 < Mdouble < length)
     */
    void setThickness(Mdouble thickness);

    /*!
     * \details Function to change the screw angular velocity
     * \param[in] omega new angular velocity
     */
    void setOmega(Mdouble omega);

    /*!
     * \details Function to change the screw angle
     * \param[in] angle new angular offset
     */
    void setOffset(Mdouble angularOffset);

    // --------------------   GET FUNCTIONS   --------------------
    // Function to get the screw parameters or related quantities
    /*!
     * \details Function that returns the particle coordinate in cylindrical frame of reference, with respect to the screw origin
     * \param[in] particle position in euclidean coordinates with respect to the screw (Vec3D)
     * \return particle position in cylindrical coordinates with respect to the screw (Vec3D)
     */
    Vec3D getCylindricalCoordinate(const Vec3D pInScrewReference) const;

    /*!
     * \details Function that returns the oriented axial distance (can be negative!) between the particle and the blade centre
     * \param[in] particle position in cylindrical coordinates with respect to the screw (Vec3D)
     * \return oriented axial distance between particle centre and blade centre (Mdouble)
     */
    Mdouble getDeltaZ(const Vec3D pInCylindricalFrame) const;

    /*!
     * \details Function that returns the normal vector with respect to the blade surface at the particle's position
     * \param[in] particle position in euclidean coordinates with respect to the screw (Vec3D)
     * \param[in] particle position in cylindrical coordinates with respect to the screw (Vec3D)
     * \return vector through the particle's centre perpendicular to the screw surface (Vec3D)
     */
    Vec3D getNormalVector(const Vec3D pInScrewReference, const Vec3D pInCylindricalFrame) const;

    /*!
     * \details Function that returns the radial vector with respect to the blade axis at the particle's position
     * \param[in] particle position in euclidean coordinates with respect to the screw (Vec3D)
     * \param[in] particle position in cylindrical coordinates with respect to the screw (Vec3D)
     * \return vector through the particle's centre perpendicular to the screw axis (Vec3D)
     */
    Vec3D getRadialVector(const Vec3D pInScrewReference, const Vec3D pInCylindricalFrame) const;

    /*!
     * \details The contact detection algorithm, where it is checked whether the given particle is in contact with the
     * screw, and if so, computes the distance between the centre of the particle and the screw and the normal of this
     * contact. Note, that the distance and normal are returned as output parameters, not as return values. Furthermore,
     * all the quantities in cylindrical coordinates are referred to the screw axis
     * \param[in] p particle for which is checked if there is an interaction with this screw (BaseParticle)
     * \param[out] distance if there is an interaction, will contain the distance between the particle and screw (Mdouble)
     * \param[out] normalength_return if there is an interaction, will contain the normal of the contact between the particle and
     * the screw (Vec3D)
     * \return Boolean to indicate whether there is an interaction: true if there is an interaction, false if there is no
     * interaction (bool)
     */
    bool getDistanceAndNormalLabCoordinates(Vec3D position, Mdouble wallInteractionRadius, Mdouble& distance, Vec3D& normal_return) const;
    
    /*!
     * \details Function that transforms the coordinates into the position-orientation frame.
     * It also takes into account the mixed species of the particle-wall system, therefore using the correct interaction distance automatically.
     * \return Boolean to indicate whether there is an interaction: true if there is an interaction, false if there is no
     * interaction (bool)
     */
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    /*!
     * \details Function that returns the current angular offset
     * \return current angular offset (Mdouble)
     */
    Mdouble getAngularOffset() const;

    /*!
     * \details Function that returns the current screw angular velocity
     * \return current angular velocity (Mdouble)
     */
    Mdouble getOmega() const;

    /*!
     * \details Function that returns the screw righthandedness
     * \return Boolean that is true if the screw is right-handed (bool)
     */
    bool getRightHandedness() const;

    // --------------------   ROTATION FUNCTIONS   --------------------
    // Function to rotate the screw
    /*!
     * \details Function to make the screw rotate by incrementing the angular offset
     * \param[in] offset the amount by which the screw should be rotated (Mdouble)
     */
    void incrementOffset(Mdouble angularOffsetIncrement);

    /*!
     * \details Function to make the screw rotate by incrementing by dt*omega_the angular offset
     * \param[in] dt time duration for which the screw should be turned (Mdouble)
     */
    void rotate(Mdouble dt);

    // --------------------   INPUT/OUTPUT FUNCTIONS   --------------------
    // Function designated for input/output purpouses
    /*!
     * \details Read the screw properties from an input-stream, for example a restart file
     * \param[in|out] is inputstream from which the screw is read (std::istream)
     */
    void read(std::istream& is) override;

    /*!
     * \details Read the screw properties from an input-stream, for example a restart file
     * \param[in|out] is input-stream from which the screw is read (std::istream)
     */
    void oldRead(std::istream& is);

    /*!
     * \details Write the screw properties to an output-stream, for example a restart file, and can also be used to write to the console by using std::cout as a parameter.
     * \param[in|out] os output-stream to which the screw properties are written (std::ostream)
     */
    void write(std::ostream& os) const override;

    /*!
     * \details Returns the string "ScrewAuger"
     * \return the std string "ScrewAuger" (std::string)
     */
    std::string getName() const final;
    
    void writeVTK(VTKContainer& vtk) const override;
    
    void writeVTK(std::string filename) const;

//    /*!
//     * \details Function that checks if a particle p is interacting with the screw at a given time stamp. If getDistanceAndNormal returns true, the particle is interacting with
//     * \param[in] p particle for which is checked if there is an interaction with this screw (BaseParticle)
//     * \param[in] timeStamp time stamp at which the interaction between particle and screw is checked (Mdouble)
//     * \param[in] interactionHandler pointer to the interaction handler (InteractionHandler*)
//     * \param[out] an interaction is added to the interaction handler, containing the particle-screw interaction detected by getDistanceAndNormal (Mdouble)
//     * \return interactions an updated list of all the interactions of particle p at time timeStamp, with the actual screw-particle interaction at the bottom of the list (std::vector<BaseInteraction*>)
//     */
//    std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler) final;
    
private:
   // explicit variables
    ///\details The centre of the lower end of the screw
    Vec3D screwOrigin_;

    ///\details The length of the screw
    Mdouble length_;

    ///\details The screw blade radius
    Mdouble bladeRadius_;

    ///\details The screw shaft radius
    Mdouble shaftRadius_;

    ///\details The number of turns.
    Mdouble numberOfTurns_;

    ///\details Screw rotational velocity around its axix in rad/s
    Mdouble omega_;

    ///\details The thickness of the screw blade
    Mdouble thickness_;

    ///\details The right handedness of the screw, i.e. the direction of the screw-blade
    bool rightHandedness_;

    // implicit variables
    ///\details Angular offset of the screw, i.e. the angle the screw has rotated
    Mdouble angularOffset_;

    ///\details The half thickness of the screw blade
    Mdouble delta_;

    ///\details The pitch length of the screw ( pitch_ = length_/numberOfTurns_)
    Mdouble pitch_;

    ///\details The rescaled pitch length of the screw (rescaledPitch_ = pitch_/2*pi)
    Mdouble rescaledPitch_;

};

#endif
