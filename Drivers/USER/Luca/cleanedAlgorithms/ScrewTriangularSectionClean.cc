
#include "ScrewTriangularSectionClean.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

// Last update 12.01.17
// Cleaned and optimized
// Needs debugging

// Default constructor
ScrewTriangularSectionClean::ScrewTriangularSectionClean()
{
    start_.setZero();
    l_ = 1.0;
    rMax_ = 1.0;
    rMin_ = 0.0;
    n_ = 1.0;
    omega_ = 1.0;
    thickness_ = 0.0;
    delta_ = 0.0;
    pitch_ = 1.0;
    h_ = 0.5 / constants::pi;
    offset_ = 0.0;
    
    logger(DEBUG, "ScrewTriangularSectionClean() constructor finished.");
}

// Copy constructor
ScrewTriangularSectionClean::ScrewTriangularSectionClean(const ScrewTriangularSectionClean& other)
        : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    rMax_ = other.rMax_;
    rMin_ = other.rMin_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    delta_ = 0.5 * other.thickness_;
    pitch_ = other.l_ / other.n_;
    h_ = 0.5 * other.l_ / (constants::pi * other.n_);
    offset_ = other.offset_;
    
    logger(DEBUG, "ScrewTriangularSectionClean(const ScrewTriangularSectionClean&) copy constructor finished.");
}

// Initialized constructor
// Parameters in: initial starting point of the screw axis, length, blade radius, shaft radius, number of turns, angular velocity, thickness of the blade.
ScrewTriangularSectionClean::ScrewTriangularSectionClean(Vec3D start, Mdouble l, Mdouble R, Mdouble r, Mdouble n,
                                                         Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = l;
    rMax_ = R;
    rMin_ = r;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
    pitch_ = l_ / n_;
    h_ = 0.5 * l_ / (constants::pi * n_);
    offset_ = 0.0;
    
    logger(DEBUG,
           "ScrewTriangularSectionClean(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
ScrewTriangularSectionClean::~ScrewTriangularSectionClean()
{
    logger(DEBUG, "~ScrewTriangularSectionClean() finished, destroyed the screw.");
}

// A pointer to the copy of the screw
ScrewTriangularSectionClean* ScrewTriangularSectionClean::copy() const
{
    return new ScrewTriangularSectionClean(*this);
}

// Function to set the screw parameters
void ScrewTriangularSectionClean::set(Vec3D start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius,
                                      Mdouble numberOfTurns, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = length;
    rMax_ = bladeRadius;
    rMin_ = shaftRadius;
    n_ = numberOfTurns;
    omega_ = omega;
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
    pitch_ = l_ / n_;
    h_ = 0.5 * l_ / (constants::pi * n_);
    offset_ = 0.0;
    
    std::cout << "\n\n screw parameters set.\n";
    std::cout << "Start : " << start_ << "\n";
    std::cout << "Length : " << l_ << "\n";
    std::cout << "Blade radius : " << rMax_ << "\n";
    std::cout << "Shaft radius : " << rMin_ << "\n";
    std::cout << "Number of turns : " << n_ << "\n";
    std::cout << "Angular velocity : " << omega_ << "\n";
    std::cout << "Thickness : " << thickness_ << "\n";
    std::cout << "Half thickness : " << delta_ << "\n";
    std::cout << "Pitch length : " << pitch_ << "\n";
    std::cout << "Rescaled length : " << h_ << "\n";
}

// Function to change the screw radius
void ScrewTriangularSectionClean::setRadius(Mdouble bladeRadius)
{
    rMax_ = bladeRadius;
}

// Function to change the thickness (and the half-thickness as well)
void ScrewTriangularSectionClean::setThickness(Mdouble thickness)
{
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
}

// Function to change the screw rotation velocity
void ScrewTriangularSectionClean::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// Function to make the screw rotate by adding \omega*t to the angular offset
void ScrewTriangularSectionClean::rotate(Mdouble dt)
{
    offset_ += omega_ * dt;
}

// Input stream
void ScrewTriangularSectionClean::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> start_
       >> dummy >> l_
       >> dummy >> rMax_
       >> dummy >> rMin_
       >> dummy >> n_
       >> dummy >> omega_
       >> dummy >> thickness_
       >> dummy >> offset_;
}

// Old input stream (no idea what this is)
void ScrewTriangularSectionClean::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> start_
       >> dummy >> l_
       >> dummy >> rMax_
       >> dummy >> rMin_
       >> dummy >> n_
       >> dummy >> omega_
       >> dummy >> offset_;
}

// Writer
void ScrewTriangularSectionClean::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << start_
       << " Length " << l_
       << " Blade Radius " << rMax_
       << " Shaft Radius " << rMin_
       << " Revolutions " << n_
       << " Omega " << omega_
       << " Thickness " << thickness_
       << " Offset " << offset_;
}

// Returns the string "ScrewTriangularSectionClean"
std::string ScrewTriangularSectionClean::getName() const
{
    return "ScrewTriangularSectionClean";
}

// The contact detection algorithm
// All the quantities in cylindrical coordinates are referred to the screw axis
bool
ScrewTriangularSectionClean::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial position of the particle
    Mdouble rho2 = pow(p.getPosition().X - start_.X, 2) + pow(p.getPosition().Y - start_.Y, 2);
    
    // if the particle is outside the cylinder containing the screw there is no collision
    if (rho2 > pow(rMax_ + p.getWallInteractionRadius(this), 2)) return false;
    if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius(this)) return false;
    if (p.getPosition().Z < start_.Z - p.getWallInteractionRadius(this)) return false;
    
    // radial position of the particle
    Mdouble rho = sqrt(rho2);
    
    // cosine of the helix angle at the particle position
    Mdouble cosEta = 1.0 / sqrt(1.0 + pow(h_ / rho, 2));
    
    // angular coordinate of the particle
    // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
    Mdouble xi = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
    if (xi < 0.0) xi += 2.0 * constants::pi;
    
    // oriented axial distance between the particle's centre and the blade centre
    Mdouble deltaZ = fmod((p.getPosition().Z - start_.Z) -
                          h_ * (xi + offset_ + 2.0 * constants::pi * ((int) ((p.getPosition().Z - start_.Z) / pitch_))),
                          pitch_);
    if (deltaZ > 0.5 * pitch_)
    {
        deltaZ -= pitch_;
    }
    else if (deltaZ < -0.5 * pitch_)
    {
        deltaZ += pitch_;
    }
    
    // component of the particle-blade_center distance normal to the blade center surface
    Mdouble deltaN = fabs(deltaZ) * cosEta;
    
    // radial component of the distance between the particle centre and the edge
    Mdouble deltaR = rho - rMax_;
    
    // ratio between the blade half thickness and the blade length (it is also equal to tanGamma)
    Mdouble lambda = delta_ / (rMax_ - rMin_);
    
    // cosine of the blade section angle
    Mdouble cosGamma = 1.0 / sqrt(1.0 + pow(lambda, 2));
    
    // if the particle-blade_surface distance is higher than the collision threshold there is no collision
//    if ((deltaN + lambda*rho)*cosGamma > p.getWallInteractionRadius(this)) return false;
    if ((deltaN + lambda * deltaR) * cosGamma > p.getWallInteractionRadius(this))
        return false;
    
    
    // trigonometric functions relative to the particle angle
    Mdouble cosXi = (p.getPosition().X - start_.X) / rho;
    Mdouble sinXi = (p.getPosition().Y - start_.Y) / rho;
    
    // normal to the blade center at the particle position
    Vec3D normalVector;
    normalVector.X = h_ * sinXi;
    normalVector.Y = -h_ * cosXi;
    normalVector.Z = rho;
    
    // takes the right direction of the vector and normalizes
    normalVector *= -deltaZ;
    normalVector /= normalVector.getLength();
    
    // radial vector at the particle position
    Vec3D radialVector;
    radialVector.X = cosXi;
    radialVector.Y = sinXi;
    radialVector.Z = 0.0;
    
    // normal to the blade surface obtained by adding the radial vector properly scaled with tanGamma and normalizing
    normalVector -= lambda * radialVector;
    normalVector /= normalVector.getLength();
    
    // check which type of collision it is
    if (deltaR - lambda * deltaN < rMax_) // collision with the blade surface
    {
        // distance between the contact point and the particle's centre
        distance = (deltaN + lambda * deltaR) * cosGamma;
        
        // if the distance is negative prints an error message
        if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SURFACE: OVERLAP > PARTICLE RADIUS\n";
        
        // the collision has no other component other than the normal to the blade surface
        normal_return = normalVector;
        
        return true;
    }
    else // collision with the blade edge
    {
        // distance between the particle centre and the and the edge
        distance = sqrt(pow(deltaN, 2) + pow(deltaR, 2));
        
        // if the particle-blade_edge distance is higher than the particle radius there is no collision
        if (distance > p.getWallInteractionRadius(this)) return false;
        
        // if the distance is negative prints an error message
        if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW EDGE: OVERLAP > PARTICLE RADIUS\n";
        
        // the normal to the edge through the particle centre
        // IMPORTANT: it needs to be normalized!
        normal_return = deltaN * normalVector - deltaR * radialVector;
        normal_return /= normal_return.getLength();
        
        return true;
    }
    
    // error message if some case has not been accounted for
    std::cout << "\nERROR: THE SCREW COLLISION FUNCTION DID NOT RETURN PROPERLY\n";
    
    return false;
}







