
// This function defines an Archimedes' screw in the z-direction from a (constant) starting point,
// a (constant) length L, a (constant) radius r, a (constant) number or revelations N and a (constant) rotation speed (rev/s).

// The screw H=H(r,theta) is a 2D surface embedded in E3, which parametric equations are:
// xH = r*Cos(theta)
// yH = r*Sin(theta)
// zH = L*theta/(2*pi*N)
// where r and theta are the parameters with r\in[rMin;rMax] and theta\in[0;2*pi*N[.

// In this version the screw has a triangular section, rather than a rectangular one.

#ifndef DualFeederCasing_BOT_H
#define DualFeederCasing_BOT_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

// Last update 19.03.19
// - creation

class DualFeederCasing_BOT : public BaseWall
{
public:

    // Default constructor
    DualFeederCasing_BOT();

    // Copy constructor
    DualFeederCasing_BOT(const DualFeederCasing_BOT& other);

    // Initialized constructor
    DualFeederCasing_BOT(Vec3D start, Mdouble length, Mdouble radius, Mdouble distance);

    // Destructor
    ~DualFeederCasing_BOT();

    // Pointer to a copy of the casing
    DualFeederCasing_BOT* copy() const final;

    // Sets the casing parameters
    void set(Vec3D start, Mdouble length, Mdouble radius, Mdouble distance);

    // Function to set the radius
    void setRadius(Mdouble radius);

    // Function to set the distance
    void setDistance(Mdouble distance);

    // Reads a screw from an input stream, for example a restart file.
    void read(std::istream& is) override;

    //Reads a screw in the old style from an input stream, for example a restart file old style.
    void oldRead(std::istream& is);

    // Writes this screw to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "screw".
    std::string getName() const final;

    // Determines if there is a collision between a BaseParticle P and the screw.
    // If so, computes the distance between particle centre and collision point and the normal to the latter and returns true, false otherwise
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const final;

    // Gets the interaction between this screw and a given BaseParticle at a given time.
    std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler);

private:
    // The mid point between the screw axes.
    Vec3D start_;

    // The length of the casing.
    Mdouble l_;

    // The radius of the two cylindrical sections of the casing.
    Mdouble r_;

    // The distance between the screw axes.
    Mdouble delta_;
};

#endif
