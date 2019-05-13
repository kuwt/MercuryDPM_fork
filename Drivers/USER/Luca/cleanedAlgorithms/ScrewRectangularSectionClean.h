#ifndef ScrewRectangularSectionClean_H
#define ScrewRectangularSectionClean_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

/*! \brief This class defines an Archimedes' screw in the z-direction from a (constant) starting point,
 a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s).

 \details The screw H=H(r,theta) is a 2D surface embedded in E3, which parametric equations are:
 xH = r*Cos(theta)
 yH = r*Sin(theta)
 zH = L*theta/(2*pi*N)
 where r and theta are the parameters with r\in[rMin;rMax] and theta\in[0;2*pi*N[.
 For a schematic drawing, see ScrewRectangularSectionClean.cc
*/
class ScrewRectangularSectionClean : public BaseWall
{
public:
    
    ///\brief Default constructor: set all variables to sensible defaults
    ScrewRectangularSectionClean();
    
    ///\brief Copy constructor
    ScrewRectangularSectionClean(const ScrewRectangularSectionClean& other);
    
    ///\brief Constructor that takes relevant parameters and initialises the screw values
    // Parameters in: initial starting point of the screw (bottom centre), length, blade radius, number of turns,
    // angular velocity, thickness of the blade.
    ScrewRectangularSectionClean(const Vec3D& start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius,
                                 Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandeness);
    
    ///\brief Destructor
    ~ScrewRectangularSectionClean();
    
    ///\brief Copy function: returns a pointer to a new copy of this screw
    ScrewRectangularSectionClean* copy() const override;
    
    // ----------   SET FUNCTIONS   ----------
    ///\brief Function that sets the screw parameters
    void set(const Vec3D& start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns,
             Mdouble omega, Mdouble thickness, bool rightHandeness);
    
    ///\brief Function that sets the screw max radius, which is the outer radius of the blade
    void setRadius(Mdouble bladeRadius);
    
    ///\brief Function that sets the thickness of the screw blade
    void setThickness(Mdouble thickness);
    
    ///\brief Function that sets the angular velocity of the screw
    void setOmega(Mdouble omega);
    
    ///\brief Function that sets the angular offset, i.e. the rotation in z-direction
    void setOffset(Mdouble offset);
    
    // ----------   GET FUNCTIONS   ------------------------------------------------------
    ///\brief Function that computes the distance of a BaseParticle to this wall and the normal of this wall if there is a collision.
    // Returns true if there is a collision
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const override;
    
    ///\brief Function that returns the current angular offset, i.e. the rotation in z-direction
    Mdouble getOffset() const;
    
    ///\brief Function that returns the right handeness of the screw
    bool getRightHandedness() const;
    
    // ----------   ROTATION FUNCTIONS   ------------------------------------------------------
    
    ///\brief Function that rotates the screw of a period dt, so that the orientation in z-direction changes with omega_*dt.
    void rotate(Mdouble dt);
    
    ///\brief Increments the angular offset, de facto substituting the move_time function
    void incrementOffset(Mdouble off);
    
    // ----------------------------------------------------------------------------------------
    
    ///\brief Reads a screw from an input stream, for example a restart file.
    void read(std::istream& is) override;
    
    ///\brief Reads a screw in the old style from an input stream, for example a restart file old style.
    void oldRead(std::istream& is);
    
    ///\brief Writes this screw to an output stream, for example a restart file.
    void write(std::ostream& os) const override;
    
    ///\brief Returns the name of the object, here the string "ScrewRectangularSectionClean".
    std::string getName() const final;

    // Gets the interaction between this screw and a given BaseParticle at a given time.
    BaseInteraction* getInteractionWith(BaseParticle* p, unsigned timeStamp,
                                                     InteractionHandler* interactionHandler) final;

private:
    ///\brief Auxiliary function for getDistanceAndNormal, this computes the normal vector, radial vector and distance in normal direction to the screw
    void
    computeNormalRadialDeltaN(const BaseParticle& p, Vec3D& normalVector, Vec3D& radialVector, Mdouble& deltaN) const;

    ///\brief Auxiliary function for computeNormalRadialDeltaN, computes the oriented axial distance between the particle's centre and the blade centre
    Mdouble computeDeltaZ(const Vec3D& distanceParticleStart, Mdouble h, Mdouble pitch) const;

    ///\brief The centre of the lower end of the screw.
    Vec3D startPosition_;

    ///\brief The length of the screw.
    Mdouble length_;

    ///\brief The outer screw blade radius.
    Mdouble rMax_;

    ///\brief The screw shaft radius.
    Mdouble rMin_;

    ///\brief The number of turns.
    Mdouble numberOfTurns_;

    ///\brief Rotation speed in rad per time unit.
    ///\todo Can probably removed? \author irana
    Mdouble omega_;

    ///\brief The thickness of the screw blade.
    Mdouble thickness_;

    ///\brief The right handedness of the screw, i.e. the direction of the screw-blade.
    bool rightHandedness_;
};

#endif
