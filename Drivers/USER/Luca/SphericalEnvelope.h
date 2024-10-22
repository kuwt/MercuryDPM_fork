
// This function defines an Archimedes' screw in the z-direction from a (constant) starting point,
// a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s).

// The screw H=H(r,theta) is a 2D surface embedded in E3, which parametric equations are:
// xH = r*Cos(theta)
// yH = r*Sin(theta)
// zH = L*theta/(2*pi*N)
// where r and theta are the parameters with r\in[rMin;rMax] and theta\in[0;2*pi*N[.

// In this version the screw has a triangular section, rather than a rectangular one.

#ifndef SphericalEnvelope_H
#define SphericalEnvelope_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

// Last update 22.03.19
// - creation

class SphericalEnvelope : public BaseWall
{
public:

    // Default constructor
    SphericalEnvelope();

    // Copy constructor
    SphericalEnvelope(const SphericalEnvelope& other);

    // Initialized constructor
    SphericalEnvelope(Vec3D start, Mdouble radius);

    // Destructor
    ~SphericalEnvelope();

    // Pointer to a copy of the envelope
    SphericalEnvelope* copy() const final;

    // Sets the envelope parameters
    void set(Vec3D start, Mdouble radius);

    // Function to set the radius
    void setRadius(Mdouble radius);

    // Function to get the radius
    double getRadius(Mdouble radius);

    // Reads a screw from an input stream, for example a restart file.
    void read(std::istream& is) override;

    // Reads a spherical shell in the old style from an input stream, for example a restart file old style.
    void oldRead(std::istream& is);

    // Writes this envelope to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "SphericalEnvelope".
    std::string getName() const final;

    // Determines if there is a collision between a BaseParticle P and the envelope.
    // If so, computes the distance between particle centre and collision point and the normal to the latter and returns true, false otherwise
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    // Gets the interaction between this envelope and a given BaseParticle at a given time.
    std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler);

private:
    // The centre of the spherical envelope.
    Vec3D start_;

    // The radius of the sphere.
    Mdouble r_;
};

#endif
