#include "OscillatingLattice1D.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// Last update 21.11.16
// ToDo: check if it is better the setOscillation(t) implementation or the setOscillation(t, t0) one
// ToDo: the well width is NEVER USED!

/*
 -- WELL WALL SECTION --
 / / / / / / / / / / / / / / / / / / / / / / / / / /
 / / / / / ............................. / / / / / /
 / / / / / :       :            :       :/ / / / / /
 / / / / / :   3   :      2     :   3   :/ / / / / /
 / / / / / :.......:____________:.......:/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :   1   |   4  |  4  |   1   :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 ___________________________________________________
 
 REGION 1: X-DIRECTION COLLISION
 REGION 2: Z-DIRECTION COLLISION
 REGION 3: EDGE (XZ PLANE DIRECTION) COLLISION
 REGION 4: UNACCESSIBLE REGION (COLLISION ERROR)
 */

// Default constructor
OscillatingLattice1D::OscillatingLattice1D()
{
    n_ = 1.0;
    l_ = 1.0;
    w_ = 1.0;
    h_ = 0.0;
    d_ = 0.0;
    z0_ = 0.0;
    
    A_ = 0.0;
    f_ = 0.0;
    
    L_ = 1.0;
    z_ = 0.0;
    omega_ = 0.0;
    
    logger(DEBUG, "OscillatingLattice1D() constructor finished.");
}

// Copy constructor
OscillatingLattice1D::OscillatingLattice1D(const OscillatingLattice1D& other)
: BaseWall(other)
{
    n_ = other.n_;
    l_ = other.l_;
    w_ = other.w_;
    h_ = other.h_;
    d_ = other.d_;
    z0_ = other.z0_;
    
    A_ = other.A_;
    f_ = other.f_;
    
    L_ = other.l_ + other.d_;
    z_ = other.z0_;
    omega_ = 2.0*constants::pi*other.f_;
    
    logger(DEBUG, "OscillatingLattice1D(const OscillatingLattice1D&) copy constructor finished.");
}

// The good constructor
// Parameters in: number of wells, well length, well width, well height, inter-well wall thickness, initial well height, oscillation amplitude, oscillation frequency
OscillatingLattice1D::OscillatingLattice1D(int numberOfwells, Mdouble wellLength, Mdouble wellWidth, Mdouble wellHeight, Mdouble wallThickness, Mdouble initialHeight, Mdouble oscillationAmplitude, Mdouble oscillationFrequency)
{
    n_ = numberOfwells;
    l_ = wellLength;
    w_ = wellWidth;
    h_ = wellHeight;
    d_ = wallThickness;
    z0_ = initialHeight;
    
    A_ = oscillationAmplitude;
    f_ = oscillationFrequency;
    
    L_ = wellLength + wallThickness;
    z_ = z0_;
    omega_ = 2.0*constants::pi*oscillationFrequency;
    
    logger(DEBUG, "OscillatingLattice1D(int, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
OscillatingLattice1D::~OscillatingLattice1D()
{
    logger(DEBUG, "~OscillatingLattice1D() finished, destroyed the Oscillating Lattice 1D.");
}

// A pointer to the copy of the oscillating lattice
OscillatingLattice1D* OscillatingLattice1D::copy() const
{
    return new OscillatingLattice1D(*this);
}

// Function to set the oscillating lattice parameters
void OscillatingLattice1D::set(int numberOfwells, Mdouble wellLength, Mdouble wellWidth, Mdouble wellHeight, Mdouble wallThickness, Mdouble initialHeight, Mdouble oscillationAmplitude, Mdouble oscillationFrequency)
{
    n_ = numberOfwells;
    l_ = wellLength;
    w_ = wellWidth;
    h_ = wellHeight;
    d_ = wallThickness;
    z0_ = initialHeight;
    
    A_ = oscillationAmplitude;
    f_ = oscillationFrequency;
    
    L_ = wellLength + wallThickness;
    z_ = z0_;
    omega_ = 2.0*constants::pi*oscillationFrequency;

    std::cout << "\nOscillating lattice parameters set.\n";
    std::cout << "Number of wells: " << n_ << "\n";
    std::cout << "Length of the wells: " << l_ << "\n";
    std::cout << "Width of the wells: " << w_ << "\n";
    std::cout << "Height of the wells: " << h_ << "\n";
    std::cout << "Inter-wells wall thickness: " << d_ << "\n";
    std::cout << "Initial well height: " << z0_ << "\n";
    std::cout << "Well section length: " << L_ << "\n";
    std::cout << "Oscillation amplitude: " << A_ << "\n";
    std::cout << "Oscillation frequency: " << f_ << "\n";
}

// Function to change the well geometry
void OscillatingLattice1D::setGeometry(Mdouble wellLength, Mdouble wellWidth, Mdouble wellHeight, Mdouble wallThickness)
{
    l_ = wellLength;
    w_ = wellWidth;
    h_ = wellHeight;
    d_ = wallThickness;
}

// Function to change the well height
void OscillatingLattice1D::setWellHeight(Mdouble wellHeight)
{
    h_ = wellHeight;
}

// Function to change the oscillation parameters
void OscillatingLattice1D::setOscillation(Mdouble oscillationAmplitude, Mdouble oscillationFrequency)
{
    A_ = oscillationAmplitude;
    f_ = oscillationFrequency;
    
    omega_ = 2.0*constants::pi*oscillationFrequency;
}

// Function to change the initial well height
void OscillatingLattice1D::setInitialHeight(Mdouble initialHeight)
{
    z0_ = initialHeight;
}

// Function to make the lattice oscillate
void OscillatingLattice1D::setOscillation(Mdouble t)
{
    z_ = z0_ + A_*sin(omega_*t);
}

// Returns the actual well height
double OscillatingLattice1D::getWellPosition()
{
    return z_;
}

// Input stream
void OscillatingLattice1D::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> n_
    >> dummy >> l_
    >> dummy >> w_
    >> dummy >> h_
    >> dummy >> d_
    >> dummy >> z0_
    >> dummy >> A_
    >> dummy >> f_;
}

// Output stream
void OscillatingLattice1D::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " N. of wells " << n_
    << " Well length " << l_
    << " Well width " << w_
    << " Well height " << h_
    << " Wall thickness " << d_
    << " Init. height " << z0_
    << " Osc. amplitude " << A_
    << " Osc. frequency " << f_;
}

// Returns the name of the driver
std::string OscillatingLattice1D::getName() const
{
    return "OscillatingLattice1D";
}

// The contact detection algorithm
// Takes as input a pointer to a particle and evaluates if there is a collision.
// If so updates the distance between the particle centre and the collision point, the collision normal and returns true.
// Returns false if there is no collision.
bool OscillatingLattice1D::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // Checks if the particle is too high to touch the walls
    if (p.getPosition().Z > z_ + h_ + p.getWallInteractionRadius(this)) return false;
    
    // Check if the particles are close enough to the walls from the side
    if (fabs(fmod(p.getPosition().X,L_) - 0.5*L_) < 0.5*d_ + p.getWallInteractionRadius(this)) return false;
    
    // collision checker, the space being divided in 4 regions
    if (p.getPosition().Z <= z_ + h_ && fabs(fmod(p.getPosition().X,L_) - 0.5*L_) >= 0.5*d_)  // REGION 1
    {
        if (fmod(p.getPosition().X,L_) >= 0.5*L_)
        {
            distance = L_ - fmod(p.getPosition().X,L_) - 0.5*d_;
            
            normal_return.X = 1.0;
            normal_return.Y = 0.0;
            normal_return.Z = 0.0;
        }
        else
        {
            distance = fmod(p.getPosition().X,L_) - 0.5*d_;
            
            normal_return.X = -1.0;
            normal_return.Y = 0.0;
            normal_return.Z = 0.0;
        }
        
        return true;
    }
    else if (p.getPosition().Z >= z_ + h_ && fabs(fmod(p.getPosition().X,L_) - 0.5*L_) <= 0.5*d_) // REGION 2
    {
        distance = p.getPosition().Z - (z_ + h_);
        
        normal_return.X = 0.0;
        normal_return.Y = 0.0;
        normal_return.Z = 1.0;
        
        return true;
    }
    else if (p.getPosition().Z > z_ + h_ && fabs(fmod(p.getPosition().X,L_) - 0.5*L_) > 0.5*d_) // REGION 3
    {
        if (fmod(p.getPosition().X,L_) >= 0.5*L_)
        {
            distance = sqrt(pow(fmod(p.getPosition().X,L_) - 0.5*d_,2) + pow(p.getPosition().Z - (z_ + h_),2));
            
            if (distance > p.getWallInteractionRadius(this)) return false;
            
            normal_return.X = (fmod(p.getPosition().X,L_) - 0.5*d_)/distance;
            normal_return.Y = 0.0;
            normal_return.Z = (p.getPosition().Z - (z_ + h_))/distance;
        }
        else
        {
            distance = sqrt(pow(L_ - fmod(p.getPosition().X,L_) - 0.5*d_,2) + pow(p.getPosition().Z - (z_ + h_),2));
            
            if (distance > p.getWallInteractionRadius(this)) return false;
            
            normal_return.X = (L_ - fmod(p.getPosition().X,L_) - 0.5*d_)/distance;
            normal_return.Y = 0.0;
            normal_return.Z = (p.getPosition().Z - (z_ + h_))/distance;
        }
        
        return true;
    }
    else    // REGION 4
    {
        std::cout << "\nTHE PARTICLE IS IN THE FORBIDDEN REGION!\n";
        
        return false;
        //        exit(0);
    }
}

// Checks for the interaction between a particle p at a time timeStamp.
// In case of interaction returns a pointer to the BaseInteraction happened between the inter-well wall and the
// BaseParticle at time timeStamp
BaseInteraction* OscillatingLattice1D::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    if (getDistanceAndNormal(*p,distance,normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        c->setContactPoint(p->getPosition() - (p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        return c;
    }
    return nullptr;
}







