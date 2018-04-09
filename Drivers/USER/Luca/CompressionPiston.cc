#include "CompressionPiston.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// This is an infinite wall used to compress particles inside a cylindrical tube.

// Last update 2.05.17

// Default constructor
CompressionPiston::CompressionPiston()
{
    height_ = 1.0;
    radius_ = 1.0;
    velocity_ = 0.0;
    maxPressure_ = 0.0;
    
    logger(DEBUG, "CompressionPiston() constructor finished.");
}

// Copy constructor
CompressionPiston::CompressionPiston(const CompressionPiston& other)
: BaseWall(other)
{
    height_ = other.height_;
    radius_ = other.radius_;
    velocity_ = other.velocity_;
    maxPressure_ = other.maxPressure_;
    
    logger(DEBUG, "CompressionPiston(const CompressionPiston&) copy constructor finished.");
}

// The good constructor
CompressionPiston::CompressionPiston(Mdouble height, Mdouble radius, Mdouble velocity, Mdouble maxPressure)
{
    height_ = height;
    radius_ = radius;
    velocity_ = velocity;
    maxPressure_ = maxPressure;
    
    logger(DEBUG, "CompressionPiston(Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
CompressionPiston::~CompressionPiston()
{
    logger(DEBUG, "~CompressionPiston() finished, destroyed the Compression Piston.");
}

// A pointer to the copy of the piston
CompressionPiston* CompressionPiston::copy() const
{
    return new CompressionPiston(*this);
}

// Function to set the piston parameters
void CompressionPiston::set(Mdouble height, Mdouble radius, Mdouble velocity, Mdouble maxPressure)
{
    height_ = height;
    radius_ = radius;
    velocity_ = velocity;
    maxPressure_ = maxPressure;

    std::cout << "\nCompression piston parameters set.\n";
    std::cout << "Initial height: " << height_ << "\n";
    std::cout << "Radius: " << radius_ << "\n";
    std::cout << "Initial velocity: " << velocity_ << "\n";
    std::cout << "Maximum pressure: " << maxPressure_ << "\n";
}

// Sets the piston height
void CompressionPiston::setHeight(Mdouble height)
{
    height_ = height;
}

// Sets the piston velocity
void CompressionPiston::setVelocity(Mdouble velocity)
{
    velocity_ = velocity;
}

// Sets the piston maximum pressure
void CompressionPiston::setMaxPressure(Mdouble maxPressure)
{
    maxPressure_ = maxPressure;
}

// Sets the piston radius
void CompressionPiston::setRadius(Mdouble radius)
{
    radius_ = radius;
}

// Moves the piston.
void CompressionPiston::movePiston(Mdouble timeStep)
{
    height_ += velocity_*timeStep;
}

// Returns the actual piston height
double CompressionPiston::getHeight() const
{
    return height_;
}

// Returns the actual piston velocity
double CompressionPiston::getVelocity()
{
    return velocity_;
}

// Reads a compressing piston from an input stream, for example a restart file.
void CompressionPiston::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> height_
    >> dummy >> radius_
    >> dummy >> velocity_
    >> dummy >> maxPressure_;
}

// Writes this compressing piston to an output stream, for example a restart file.
void CompressionPiston::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " Height " << height_
    << " Radius " << radius_
    << " Velocity " << velocity_
    << " Maximum pressure " << maxPressure_;
}

// Returns the name of the driver
std::string CompressionPiston::getName() const
{
    return "CompressionPiston";
}

// The contact detection algorithm
// Takes as input a pointer to a particle and evaluates if there is a collision.
// If so updates the distance between the particle centre and the collision point, the collision normal and returns true.
// Returns false if there is no collision.
bool CompressionPiston::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // The particle centre-piston distance
    distance = height_ - p.getPosition().Z;
    
    // Checks if the particle is in contact with the wall
    if (distance > p.getWallInteractionRadius()) return false;
    
    // If the compression is too high throws a warning but doesn't exit
    if (distance <= 0.0)
    {
        std::cout << "\nWARNING: OVER-COMPRESSION!\n";
    }
    
    // The normal is in the Z direction
    normal_return.X = 0.0;
    normal_return.Y = 0.0;
    normal_return.Z = 1.0;
    
    return true;
}

// Checks for the interaction between a particle p at a time timeStamp.
// In case of interaction returns a pointer to the BaseInteraction happened between the inter-well wall and the
// BaseParticle at time timeStamp
std::vector<BaseInteraction *> CompressionPiston::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    std::vector<BaseInteraction*> interactions;
    if (getDistanceAndNormal(*p,distance,normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        interactions.push_back(c);
    }
    return interactions;
}









