// This function defines an Archimedes' screw in the z-direction from a (constant) starting point,
// a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s).

// The screw H=H(r,theta) is a 2D surface embedded in E3, which parametric equations are:
// xH = r*Cos(theta)
// yH = r*Sin(theta)
// zH = L*theta/(2*pi*N)
// where r and theta are the parameters with r\in[rMin;rMax] and theta\in[0;2*pi*N[.

#ifndef ScrewRectangularSection_H
#define ScrewRectangularSection_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class ScrewRectangularSection : public BaseWall
{
public:
    
    // Default constructor
    ScrewRectangularSection();
    
    // Copy constructor
    ScrewRectangularSection(const ScrewRectangularSection& other);

    // Initialized constructor
    // Parameters in: initial starting point of the screw, length, blade radius, number of turns, angular velocity, thickness of the blade.
    ScrewRectangularSection(Vec3D start, Mdouble l, Mdouble R, Mdouble n, Mdouble omega, Mdouble thickness);
    
    // Destructor
    ~ScrewRectangularSection();
    
    // Pointer to a copy of the screw
    ScrewRectangularSection* copy() const final;
    
    // Sets the screw parameters
    void set(Vec3D Start, Mdouble length, Mdouble bladeRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness);
    
    // Sets the screw max radius
    void setRadius(Mdouble bladeRadius);
    
    // Sets the screw thickness
    void setThickness(Mdouble thickness);
    
    // Sets the screw angular velocity
    void setOmega(Mdouble omega);

    // Determines if there is a collision between a BaseParticle P and the screw.
    // If so, computes the distance between particle centre and collision point and the normal to the latter and returns true:
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    // Rotates the screw of a period dt, so that the offset_ changes with omega_*dt.
    void rotate(Mdouble dt);
    
    // Increments the angular offset, de facto substituting the move_time function
    void incrementOffset(Mdouble off);
    
    // Returns the current angular offset
    Mdouble getOffset();
    
    // Sets the angular offset
    void setOffset(Mdouble off);

    // Reads a screw from an input stream, for example a restart file.
    void read(std::istream& is) override;

    //Reads a screw in the old style from an input stream, for example a restart file old style.
    void oldRead(std::istream& is);

    // Writes this screw to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "ScrewRectangularSection".
    std::string getName() const final;

    // Gets the interaction between this screw and a given BaseParticle at a given time.
    BaseInteraction* getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) final;
    
private:
    // The centre of the lower end of the screw.
    Vec3D start_;
    // The length of the screw.
    Mdouble l_;
    // The outer radius of the screw.
    Mdouble rMax_;
    // The number of turns.
    Mdouble n_;
    // Rotation speed in rad/s.
    Mdouble omega_;
    // The angle that describes how much the screw has turned.
    Mdouble offset_;
    // The thickness of the screw.
    Mdouble thickness_;
    // The half-thickness of the screw.
    Mdouble delta_;
    // The pitch length of the screw (l_/n_).
    Mdouble pitch_;
    // The rescaled length of the screw (l_/2*pi*n_).
    Mdouble h_;
};

#endif
