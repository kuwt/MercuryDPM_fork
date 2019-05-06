#ifndef CompressionPistonSurface_H
#define CompressionPistonSurface_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class CompressionPistonSurface : public BaseWall
{
public:

    // Default constructor
    CompressionPistonSurface();

    // Copy constructor
    CompressionPistonSurface(const CompressionPistonSurface& other);

    // Good constructor
    CompressionPistonSurface(Mdouble position, Mdouble height, Mdouble thickness, Mdouble angle);

    // Destroyer
    ~CompressionPistonSurface();

    // Pointer to a copy of the piston
    CompressionPistonSurface* copy() const final;

    // Sets the piston parameters
    void set(Mdouble position, Mdouble height, Mdouble thickness, Mdouble angle);

    // Moves the piston along the axis
    void move(Mdouble dh);

    // Rotates the piston
    void rotate(Mdouble dTheta);

    // Returns the actual piston position
    double getPosition() const;

    // Returns the actual piston angle
    double getAngle() const;

    // Reads a compressing piston from an input stream, for example a restart file.
    void read(std::istream& is) override;

    // Writes this compressing piston to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "CompressionPistonSurface".
    std::string getName() const final;

    // Determines if there is a collision between a BaseParticle P and the piston.
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const final;

    // Gets the interaction between the piston and a given BaseParticle at a given time.
    std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, Mdouble timeStamp, InteractionHandler* interactionHandler);

private:
   // The piston position.
   Mdouble position_;

   // The piston surface height.
   Mdouble height_;

   // The piston thickness.
   Mdouble thickness_;

   // The piston angle.
   Mdouble angle_;
};

#endif
