#include "Helicoid03.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// Last update 19.11.16
// ToDo: check for different boundaries depending on maximum possible overlap and particle radius
// ToDo: put checks on quantities (e.g. distances must be positive, etc)
// ToDo: check if oldRead can be deleted (I guess so...)

/*
                 -- BLADE SECTION --
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
 
 REGION 1: NORMAL COLLISION
 REGION 2: RADIAL COLLISION
 REGION 3: EDGE (MIXED) COLLISION
 REGION 4: UNACCESSIBLE REGION (COLLISION ERROR)
 */

// Default constructor
Helicoid03::Helicoid03()
{
    start_.setZero();
    l_ = 1.0;
    maxR_ = 1.0;
    n_ = 1.0;
    omega_ = 1.0;
    thickness_ = 0.0;
    delta_ = 0.0;
    pitch_ = 1.0;
    h_ = 0.5/constants::pi;
    offset_ = 0.0;
    
    logger(DEBUG, "Helicoid03() constructor finished.");
}

// Copy constructor
Helicoid03::Helicoid03(const Helicoid03& other)
    : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    maxR_ = other.maxR_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    delta_ = 0.5*other.thickness_;
    pitch_ = other.l_/other.n_;
    h_ = 0.5*other.l_/(constants::pi * other.n_);
    offset_ = other.offset_;
    
    logger(DEBUG, "Helicoid03(const Helicoid03&) copy constructor finished.");
}

// The good constructor
// Parameters in: the initial starting point of the screw axis, the length, the radius, the number of revelations,
// the angular velocity, the thickness of the blade.
Helicoid03::Helicoid03(Vec3D start, Mdouble l, Mdouble R, Mdouble n, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = l;
    maxR_ = R;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    delta_ = 0.5*thickness_;
    pitch_ = l_/n_;
    h_ = 0.5*l_/(constants::pi * n_);
    offset_ = 0.0;
    
    logger(DEBUG, "Helicoid03(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
Helicoid03::~Helicoid03()
{
    logger(DEBUG, "~Helicoid03() finished, destroyed the Helicoid.");
}

// A pointer to the copy of the helicoid
Helicoid03* Helicoid03::copy() const
{
    return new Helicoid03(*this);
}

// Function to set the Helicoid parameters
void Helicoid03::set(Vec3D start, Mdouble length, Mdouble bladeRadius, Mdouble numberOfRevelations, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = length;
    maxR_ = bladeRadius;
    n_ = numberOfRevelations;
    omega_ = omega;
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
    pitch_ = l_/n_;
    h_ = 0.5 * l_/(constants::pi * n_);
    offset_ = 0.0;
    
    std::cout << "\n\n Helicoid parameters set.\n";
    std::cout << "Start : " << start_ << "\n";
    std::cout << "Length : " << l_ << "\n";
    std::cout << "Blade radius : " << maxR_ << "\n";
    std::cout << "Number of turns : " << n_ << "\n";
    std::cout << "Angular velocity : " << omega_ << "\n";
    std::cout << "Thickness : " << thickness_ << "\n";
    std::cout << "Half thickness : " << delta_ << "\n";
    std::cout << "Pitch length : " << pitch_ << "\n";
    std::cout << "Rescaled length : " << h_ << "\n";
    std::cout << "Initial angular offset : " << offset_ << "\n";
}

// Function to change the helicoid radius
void Helicoid03::setRadius(Mdouble bladeRadius)
{
    maxR_ = bladeRadius;
}

// Function to change the thickness (and the half-thickness as well)
void Helicoid03::setThickness(Mdouble thickness)
{
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
}

// Function to change the helicoid rotation velocity
void Helicoid03::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// Function to make the screw rotate by adding \omega*t to the angular offset
void Helicoid03::move_time(Mdouble dt)
{
    offset_ += omega_ * dt;
}

Mdouble Helicoid03::getOffset()
{
    return offset_;
}

void Helicoid03::setOffset(Mdouble off)
{
    offset_ = off;
}

// Function to make the screw rotate by incrementing the angular offset (should substitute the move_time function)
void Helicoid03::incrementOffset(Mdouble off)
{
    offset_ += off;
}

// Helicoid Input stream
void Helicoid03::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> start_
    >> dummy >> l_
    >> dummy >> maxR_
    >> dummy >> n_
    >> dummy >> omega_
    >> dummy >> thickness_
    >> dummy >> offset_;
}

// I think I can kill this since it refers to the old implementation I have no idea about
void Helicoid03::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> start_
    >> dummy >> l_
    >> dummy >> maxR_
    >> dummy >> n_
    >> dummy >> omega_
    >> dummy >> offset_;
}

// Helicoid writer
void Helicoid03::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << start_
    << " Length " << l_
    << " Radius " << maxR_
    << " Revolutions " << n_
    << " Omega " << omega_
    << " Thickness " << thickness_
    << " Offset " << offset_;
}

// Returns the string "Helicoid"
std::string Helicoid03::getName() const
{
    return "Helicoid";
}

// The contact detection algorithm
// Takes as input a pointer to a particle and evaluates if there is a collision.
// If so updates the distance between the particle centre and the collision point, the collision normal and returns true.
// Returns false if there is no collision.
bool Helicoid03::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial distance of the particle from the helicoid axis
    Mdouble rho2 = pow(p.getPosition().X - start_.X, 2) + pow(p.getPosition().Y - start_.Y, 2);
    
    // if the particle is outside of the cylinder that contains the helicoid returns false
    if (rho2 > pow(maxR_ + p.getWallInteractionRadius(this), 2)) return false;
    if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius(this) + delta_) return false;
    if (p.getPosition().Z < start_.Z - p.getWallInteractionRadius(this) - delta_) return false;
    
    // radial distance of the particle from the helicoid axis
    Mdouble rho = sqrt(rho2);
    
    // cosine of the helix angle
    Mdouble cosEta = 1.0/sqrt(1.0+pow(h_/rho,2));
    
    // angular coordinate of the particle in the cylindrical frame of reference centered on the helicoid axis
    // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
    Mdouble xi = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
    if (xi < 0.0) xi += 2.0*constants::pi;
    
    // finding the "oriented" axial distance between the particle's centre and the blade centre
    int k = (int)((p.getPosition().Z - start_.Z)/pitch_);
    
    Mdouble deltaZ = fmod((p.getPosition().Z - start_.Z) - h_*(xi + 2.0*constants::pi*k) - h_*offset_,pitch_);
    if (deltaZ > 0.5*pitch_)
    {
        deltaZ -= pitch_;
    }
    else if (deltaZ < -0.5*pitch_)
    {
        deltaZ += pitch_;
    }
    
    // normal collision checkcheck for possible collision
    if (fabs(deltaZ)*cosEta > p.getWallInteractionRadius(this) + delta_) return false;
    
    // trigonometric functions relative to the particle position
    Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
    Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
    
    // the component normal to the blade and the radial one at the collision point
    Vec3D normalVector;
    Vec3D radialVector;
    
    // the two components of the "oriented" distance between the particle centre and the contact point on the blade
    Mdouble normalDistance;
    Mdouble radialDistance;
    
    // trigonometric functions relative to the particle position with respect to the edge of the blade
    Mdouble cosGamma;
    Mdouble sinGamma;
    
    // collision checker, the space being divided in 4 regions
    if (rho <= maxR_ && fabs(deltaZ)*cosEta >= delta_)  // REGION 1
    {
        // the components of the distance
        normalDistance = fabs(deltaZ)*cosEta - delta_;
        radialDistance = 0.0;
        
        // the normal to the blade at the particle position
        normalVector.X = h_*sinXi;
        normalVector.Y = -h_*cosXi;
        normalVector.Z = rho;
        
        // takes the right direction of the vector and normalizes
        normalVector *= - deltaZ;
        normalVector /= normalVector.getLength();
        
        // the particle centre distance from the blade
        distance = normalDistance;
        
        // the normal at the collision point
        normal_return = normalVector;
        
        return true;
    }
    else if (rho >= maxR_ && fabs(deltaZ)*cosEta <= delta_) // REGION 2
    {
        // the components of the distance
        normalDistance = 0.0;
        radialDistance = rho - maxR_;
        
        // the radial component of the normal to the blade at the particle position
        radialVector.X = cosXi;
        radialVector.Y = sinXi;
        radialVector.Z = 0.0;
        
        // the particle centre distance from the blade
        distance = radialDistance;
        
        // the normal at the collision point
        normal_return = -radialVector;
        
        return true;
    }
    else if (rho >= maxR_ && fabs(deltaZ)*cosEta >= delta_) // REGION 3
    {
        // the components of the distance
        normalDistance = fabs(deltaZ)*cosEta - delta_;
        radialDistance = rho - maxR_;
        
        cosGamma = 1.0/sqrt(1.0+pow(radialDistance/normalDistance,2));
        sinGamma = 1.0/sqrt(1.0+pow(normalDistance/radialDistance,2));
        
        if (normalDistance > p.getWallInteractionRadius(this)*cosGamma) return false;
        
        // the normal to the blade at the particle position
        normalVector.X = h_*sinXi;
        normalVector.Y = -h_*cosXi;
        normalVector.Z = rho;
        
        // the radial component of the normal to the blade at the particle position
        radialVector.X = cosXi;
        radialVector.Y = sinXi;
        radialVector.Z = 0.0;
        
        // takes the right direction of the vector and normalizes
        normalVector *= - deltaZ;
        normalVector /= normalVector.getLength();
        
        // the particle centre distance from the blade
        distance = sqrt(pow(normalDistance,2) + pow(radialDistance,2));
        
        // the normal at the collision point
        normal_return = normalVector*cosGamma - radialVector*sinGamma;
        normal_return /= normal_return.getLength();
        
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
// In case of interaction returns a pointer to the BaseInteraction happened between the Helicoid and the
// BaseParticle at time timeStamp
BaseInteraction* Helicoid03::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
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









