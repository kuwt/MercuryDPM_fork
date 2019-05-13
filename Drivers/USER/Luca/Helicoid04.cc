
#include "Helicoid04.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// Last update 05.01.17
// ToDo: check for different boundaries depending on maximum possible overlap and particle radius
// ToDo: i believe it can be re-written to have a unique loop just defining if (...) cosBeta = whatever else cosBeta = 1.0;
// ToDo: put checks on quantities (e.g. distances must be positive, etc)
// ToDo: check if oldRead can be deleted (I guess so...)

// ToDo: implement the minR value throughout this

// Default (bad) constructor
Helicoid04::Helicoid04()
{
    start_.setZero();
    l_ = 1.0;
    maxR_ = 1.0;
    minR_ = 0.0;
    n_ = 1.0;
    omega_ = 1.0;
    thickness_ = 0.0;
    delta_ = 0.0;
    pitch_ = 1.0;
    h_ = 0.5/constants::pi;
    offset_ = 0.0;
    
    logger(DEBUG, "Helicoid04() constructor finished.");
}

// Copy constructor
Helicoid04::Helicoid04(const Helicoid04& other)
    : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    maxR_ = other.maxR_;
    minR_ = other.minR_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    delta_ = 0.5*other.thickness_;
    pitch_ = other.l_/other.n_;
    h_ = 0.5*other.l_/(constants::pi * other.n_);
    offset_ = other.offset_;
    
    logger(DEBUG, "Helicoid04(const Helicoid04&) copy constructor finished.");
}

// The good constructor
// Parameters in: the initial starting point of the screw axis, the length, the blade radius, the shaft radius, the number of revelations,
// the angular velocity, the thickness of the blade.
Helicoid04::Helicoid04(Vec3D start, Mdouble l, Mdouble R, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = l;
    maxR_ = R;
    minR_ = r;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    delta_ = 0.5*thickness_;
    pitch_ = l_/n_;
    h_ = 0.5*l_/(constants::pi * n_);
    offset_ = 0.0;
    
    logger(DEBUG, "Helicoid04(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destroyer of worlds
Helicoid04::~Helicoid04()
{
    logger(DEBUG, "~Helicoid04() finished, destroyed the Helicoid.");
}

// A pointer to the copy of the helicoid
Helicoid04* Helicoid04::copy() const
{
    return new Helicoid04(*this);
}

// Function to set the Helicoid parameters
void Helicoid04::set(Vec3D start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfRevelations, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = length;
    maxR_ = bladeRadius;
    minR_ = shaftRadius;
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
    std::cout << "Shaft radius : " << minR_ << "\n";
    std::cout << "Number of turns : " << n_ << "\n";
    std::cout << "Angular velocity : " << omega_ << "\n";
    std::cout << "Thickness : " << thickness_ << "\n";
    std::cout << "Half thickness : " << delta_ << "\n";
    std::cout << "Pitch length : " << pitch_ << "\n";
    std::cout << "Rescaled length : " << h_ << "\n";
}

// Function to change the helicoid radius
void Helicoid04::setRadius(Mdouble bladeRadius)
{
    maxR_ = bladeRadius;
}

// Function to change the thickness (and the half-thickness as well)
void Helicoid04::setThickness(Mdouble thickness)
{
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
}

// Function to change the helicoid rotation velocity
void Helicoid04::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// Function to make the screw rotate by adding \omega*t to the angular offset
void Helicoid04::move_time(Mdouble dt)
{
    offset_ += omega_ * dt;
}

// Helicoid Input stream
void Helicoid04::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> start_
    >> dummy >> l_
    >> dummy >> maxR_
    >> dummy >> minR_
    >> dummy >> n_
    >> dummy >> omega_
    >> dummy >> thickness_
    >> dummy >> offset_;
}

// I think I can kill this since it refers to the old implementation I have no idea about
void Helicoid04::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> start_
    >> dummy >> l_
    >> dummy >> maxR_
    >> dummy >> minR_
    >> dummy >> n_
    >> dummy >> omega_
    >> dummy >> offset_;
}

// Helicoid writer
void Helicoid04::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << start_
    << " Length " << l_
    << " Blade Radius " << maxR_
    << " Shaft Radius " << minR_
    << " Revolutions " << n_
    << " Omega " << omega_
    << " Thickness " << thickness_
    << " Offset " << offset_;
}

// Returns the string "Helicoid"
std::string Helicoid04::getName() const
{
    return "Helicoid";
}

// The contact detection algorithm
// Takes as input a pointer to a particle and evaluates if there is a collision.
// If so updates the distance between the particle centre and the collision point, the collision normal and returns true.
// Returns false if there is no collision.
bool Helicoid04::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial distance of the particle from the helicoid axis
    Mdouble rho2 = pow(p.getPosition().X - start_.X, 2) + pow(p.getPosition().Y - start_.Y, 2);
    
    // if the particle is outside of the cylinder that contains the helicoid returns false
    if (rho2 > pow(maxR_ + p.getWallInteractionRadius(this), 2)) return false;
    if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius(this)) return false;
    if (p.getPosition().Z < start_.Z - p.getWallInteractionRadius(this)) return false;
    
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
    
    // the modulation factor shaping the blade
    Mdouble modFactor = 1.0 - (rho - minR_)/(maxR_ - minR_);
    
    // check for the collision threshold
    if (fabs(deltaZ)*cosEta > p.getWallInteractionRadius(this) + delta_*modFactor) return false;
    
    // trigonometric functions relative to the particle position
    Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
    Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
    
    // the component normal to the centre of the blade
    Vec3D normalVector;
    
    // the normal to the center of the blade at the particle position
    normalVector.X = h_*sinXi;
    normalVector.Y = -h_*cosXi;
    normalVector.Z = rho;
    
    // radial component of the normal to the blade surface
    Vec3D radialVector;
    
    // radial vector in cartesian coordinates
    radialVector.X = cosXi;
    radialVector.Y = sinXi;
    radialVector.Z = 0.0;
    
    // the "oriented" distance between the particle centre and the centre on the blade
    Mdouble normalDistance = deltaZ*cosEta;
    
    // trigonometric functions relative to the normal_to_the_blade vs. radial_vector angle
    Mdouble cosGamma;
    Mdouble sinGamma;
    
    // collision computation
    if (rho >= maxR_) // COLLISION WITH THE EDGE
    {
        // the distance between the particle centre and the contact point on the edge
        Mdouble radialDistance = rho - maxR_;
        
        // updates the angle
        cosGamma = 1.0/sqrt(1.0+pow(radialDistance/normalDistance,2));
        sinGamma = 1.0/sqrt(1.0+pow(normalDistance/radialDistance,2));
        
        // checks for collision with the edge
        if (fabs(normalDistance) > p.getWallInteractionRadius(this)*cosGamma) return false;
        
        // distance between the contact point and the particle's centre
        distance = sqrt(pow(normalDistance,2) + pow(radialDistance,2));
        
        // takes the right direction of the vector and normalizes
        normalVector *= - deltaZ;
        normalVector /= normalVector.getLength();
        
        // the normal is the weighted sum between the normalized radial and normal components
        normal_return = normalVector*cosGamma - radialVector*sinGamma;
        
        return true;
    }
    else    // COLLISION WITH THE SIDE OF THE BLADE
    {
        // the half thickness to blade radius ratio
        Mdouble lambda = delta_/(maxR_ - minR_);
        
        // updates the angle
        // Note: it won't depend on the position of the particle since the shape of the blade is not changing
        cosGamma = 1.0/sqrt(1.0+pow(lambda,2));
        sinGamma = 1.0/sqrt(1.0+pow(1.0/lambda,2));
        
        // distance between the contact point and the particle's centre
        distance = (fabs(normalDistance) - delta_*modFactor)*cosGamma;
        
        // takes the right direction of the vector normal to the center of the blade and normalizes
        normalVector *= - deltaZ;
        normalVector /= normalVector.getLength();
        
        // the normal is the weighted sum between the normalized radial and normal components
        normal_return = normalVector*cosGamma - radialVector*sinGamma;
        normal_return /= normal_return.getLength();
        
        return true;
    }
}

// Checks for the interaction between a particle p at a time timeStamp.
// In case of interaction returns a pointer to the BaseInteraction happened between the Helicoid and the
// BaseParticle at time timeStamp
BaseInteraction* Helicoid04::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
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








