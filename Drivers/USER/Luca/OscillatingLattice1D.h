#ifndef OscillatingLattice1D_H
#define OscillatingLattice1D_H

#include "Walls/BaseWall.h"
#include "Math/Vector.h"
#include "Math/ExtendedMath.h"

class OscillatingLattice1D : public BaseWall
{
public:
    
    // Default constructor
    OscillatingLattice1D();
    
    // Copy constructor
    OscillatingLattice1D(const OscillatingLattice1D& other);

    // Good constructor
    OscillatingLattice1D(int numberOfwells, Mdouble wellLength, Mdouble wellWidth, Mdouble wellHeight, Mdouble wallThickness, Mdouble initialHeight, Mdouble oscillationAmplitude, Mdouble oscillationFrequency);
    
    // Destroyer
    ~OscillatingLattice1D();
    
    // Pointer to a copy of the Helicoid
    OscillatingLattice1D* copy() const final;
    
    // Sets the oscillating lattice parameters
    void set(int numberOfwells, Mdouble wellLength, Mdouble wellWidth, Mdouble wellHeight, Mdouble wallThickness, Mdouble initialHeight, Mdouble oscillationAmplitude, Mdouble oscillationFrequency);
    
    // Sets the oscillating lattice geometry
    void setGeometry(Mdouble wellLength, Mdouble wellWidth, Mdouble wellHeight, Mdouble wallThickness);
    
    // Sets the well height
    void setWellHeight(Mdouble wellHeight);
    
    // Set the oscillation parameters
    void setOscillation(Mdouble oscillationAmplitude, Mdouble oscillationFrequency);
    
    // Set the initial height of the wells
    void setInitialHeight(Mdouble initialHeight);
    
    // Moves the wells
    void setOscillation(Mdouble t);
    
    // Returns the actual well height
    double getWellPosition();
    
    // Reads an oscillating lattice from an input stream, for example a restart file.
    void read(std::istream& is) override;

    // Writes this oscillating lattice to an output stream, for example a restart file.
    void write(std::ostream& os) const override;

    // Returns the name of the object, here the string "OscillatingLattice1D".
    std::string getName() const final;
    
    // Determines if there is a collision between a BaseParticle P and one of the inter-well walls.
    bool getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const final;
    
    // Gets the interaction between the inter-well walls and a given BaseParticle at a given time.
    std::vector<BaseInteraction *> getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) final;

    // The time-dependent height of the well.
    Mdouble z_;
    
private:
    // The number of wells in the lattice.
    int n_;
    // The well length.
    Mdouble l_;
    // The well width.
    Mdouble w_;
    // The well height.
    Mdouble h_;
    // The inter-well wall thickness.
    Mdouble d_;
    // The initial well height.
    Mdouble z0_;
    
    // The oscillation amplitude.
    Mdouble A_;
    // The oscillation frequency.
    Mdouble f_;
    
    // The section length.
    Mdouble L_;
    // The oscillation rescaled frequency (angular velocity).
    Mdouble omega_;
};

#endif
