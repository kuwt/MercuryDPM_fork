#ifndef CompressionPiston_H
#define CompressionPiston_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class CompressionPiston : public BaseWall
{
public:
    
    // Default constructor
    CompressionPiston();
    
    // Copy constructor
    CompressionPiston(const CompressionPiston& other);

    // Good constructor
    CompressionPiston(Mdouble height, Mdouble radius, Mdouble velocity, Mdouble maxPressure);
    
    // Destroyer
    ~CompressionPiston();
    
    // Pointer to a copy of the piston
    CompressionPiston* copy() const final;
    
    // Sets the piston parameters
    void set(Mdouble height, Mdouble radius, Mdouble velocity, Mdouble maxPressure);
    
    // Sets the piston height
    void setHeight(Mdouble height);
    
    // Sets the piston velocity
    void setVelocity(Mdouble velocity);
    
    // Sets the piston maximum pressure
    void setMaxPressure(Mdouble maxPressure);
    
    // Sets the piston radius
    void setRadius(Mdouble radius);
    
    // Moves the piston by an amount heightIncrement.
    void movePiston(Mdouble heightIncrement);
    
    // Returns the actual piston height
    double getHeight() const;
    
    // Returns the actual piston velocity
    double getVelocity();
    
    // Reads a compressing piston from an input stream, for example a restart file.
    void read(std::istream& is) override;

    // Writes this compressing piston to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "CompressionPiston".
    std::string getName() const final;
    
    // Determines if there is a collision between a BaseParticle P and the piston.
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const final;
    
    // Gets the interaction between the piston and a given BaseParticle at a given time.
    BaseInteraction* getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) final;
    
private:
    // The piston height.
    Mdouble height_;
    
    // The piston radius.
    Mdouble radius_;
    
    // The piston velocity.
    Mdouble velocity_;
    
    // The piston maximum pressure.
    Mdouble maxPressure_;
    
    // A parameter used to change and reverse the velocity automatically.
    bool velocityInverted_;
};

#endif
