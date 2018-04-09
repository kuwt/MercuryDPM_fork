// This function defines an Archimedes' Helicoid in the z-direction from a (constant) starting point,
// a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s).

// The helicoid H=H(r,theta) is a 2D surface embedded in E3, which parametric equations are:
// xH = r*Cos(theta)
// yH = r*Sin(theta)
// zH = L*theta/(2*pi*N)
// where r and theta are the parameters with r\in[rMin;rMax] and theta\in[0;2*pi*N[.

#ifndef Helicoid03bis_H
#define Helicoid03bis_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class Helicoid03bis : public BaseWall
{
public:
    
    // Default constructor
    Helicoid03bis();
    
    // Copy constructor
    Helicoid03bis(const Helicoid03bis& other);

    // Good constructor
    Helicoid03bis(Vec3D start, Mdouble l, Mdouble R, Mdouble n, Mdouble omega, Mdouble thickness);
    
    // Destroyer
    ~Helicoid03bis();
    
    // Pointer to a copy of the Helicoid
    Helicoid03bis* copy() const final;
    
    // Sets the Helicoid parameters
    void set(Vec3D Start, Mdouble length, Mdouble bladeRadius, Mdouble numberOfRevelations, Mdouble omega, Mdouble thickness);
    
    // Sets the Helicoid max radius
    void setRadius(Mdouble bladeRadius);
    
    // Sets the Helicoid thickness
    void setThickness(Mdouble thickness);
    
    // Sets the Helicoid angular velocity
    void setOmega(Mdouble omega);

    // Determines if there is a collision between a BaseParticle P and the Helicoid.
    // If so, computes the distance between particle centre and collision point and the normal to the latter and returns true:
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    // Rotates the Helicoid of a period dt, so that the offset_ changes with omega_*dt.
    void move_time(Mdouble dt);
    
    // Increments the angular offset, de facto substituting the move_time function
    void incrementOffset(Mdouble off);
    
    Mdouble getOffset();
    
    void setOffset(Mdouble off);

    // Reads a Helicoid from an input stream, for example a restart file.
    void read(std::istream& is) override;

    //Reads a Helicoid in the old style from an input stream, for example a restart file old style.
    void oldRead(std::istream& is);

    // Writes this Helicoid to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "Helicoid".
    std::string getName() const final;

    // Gets the interaction between this Helicoid and a given BaseParticle at a given time.
    std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) final;

    // The angle that describes how much the Helicoid has turned.
    Mdouble offset_;
    
private:
    // The centre of the lower end of the Helicoid.
    Vec3D start_;
    // The length of the Helicoid.
    Mdouble l_;
    // The outer radius of the Helicoid.
    Mdouble maxR_;
    // The number of revelations.
    Mdouble n_;
    // Rotation speed in rad/s.
    Mdouble omega_;
    
    // The thickness of the Helicoid.
    Mdouble thickness_;
    // The half-thickness of the Helicoid.
    Mdouble delta_;
    // The pitch length of the Helicoid (l_/n_).
    Mdouble pitch_;
    // The rescaled length of the Helicoid (l_/2*pi*n_).
    Mdouble h_;
};

#endif
