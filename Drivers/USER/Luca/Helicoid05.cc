
#include "Helicoid05.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"

// Last update 22.02.17
/*
 -- BLADE SECTION --
 / / / / / / / / / / / / / / / / / / / / / / / / / /
 / / / / / ............................. / / / / / /
 / / / / / :       :            :       :/ / / / / /
 / / / / / :   2   :      1     :   2   :/ / / / / /
 / / / / / :.......:____________:.......:/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :   3   |      |     |   3   :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 / / / / / :       |      |     |       :/ / / / / /
 ..................|   E  |  E  |...................
           :       |      |     |       :
      4    :   5   |      |     |   5   :     4
 __________________/      |     \___________________
           :              |             :           
      E    :       E      |     E       :     E
           :              |             :
 
 REGION 1: BLADE SIDE COLLISION
 REGION 2: EDGE COLLISION
 REGION 3: BLADE SURFACE COLLISION
 REGION 4: SHAFT COLLISION
 REGION 5: CORNER COLLISION (ROUNDED EDGE)
 REGION E: UNACCESSIBLE REGION (COLLISION ERROR)
 ELSE: NO COLLISION
 */

// ToDo:
// - add the screw paraemter check and exit(0) if they are not met with an error message
// - implement the finite screw
// - implement the angular velocity update inside this driver
// - implement the right/left handed screw

// Default constructor
Helicoid05::Helicoid05()
{
    // user defined variables
    start_.setZero();
    l_ = 1.0;
    rMax_ = 1.0;
    rMin_ = 0.5*rMax_;
    n_ = 1.0;
    omega_ = 1.0;
    thickness_ = 0.0;
    offset_ = 0.0;
    rightHandeness_ = 1;
    // derived variables
    delta_ = 0.0;
    pitch_ = 1.0;
    h_ = 0.5/constants::pi;
    
    logger(DEBUG, "Helicoid05() constructor finished.");
}

// Copy constructor
Helicoid05::Helicoid05(const Helicoid05& other)
    : BaseWall(other)
{
    // user defined variables
    start_ = other.start_;
    l_ = other.l_;
    rMax_ = other.rMax_;
    rMin_ = other.rMin_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    offset_ = other.offset_;
    rightHandeness_ = other.rightHandeness_;
    // derived variables
    delta_ = 0.5*other.thickness_;
    pitch_ = other.l_/other.n_;
    h_ = 0.5*other.l_/(constants::pi * other.n_);
    
    logger(DEBUG, "Helicoid05(const Helicoid05&) copy constructor finished.");
}

// Initialized constructor
// Parameters in: initial starting point of the screw, length, blade radius, number of turns, angular velocity, thickness of the blade.
Helicoid05::Helicoid05(Vec3D start, Mdouble l, Mdouble R, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness, bool rightHandeness)
{
    // user defined variables
    start_ = start;
    l_ = l;
    rMax_ = R;
    rMin_ = r;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    offset_ = 0.0;
    rightHandeness_ = rightHandeness;
    // derived variables
    delta_ = 0.5*thickness_;
    pitch_ = l_/n_;
    h_ = 0.5*l_/(constants::pi * n_);
    
    logger(DEBUG, "Helicoid05(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble, bool) constructor finished.");
}

// Destructor
Helicoid05::~Helicoid05()
{
    logger(DEBUG, "~Helicoid05() finished, destroyed the screw.");
}

// A pointer to the copy of the screw
Helicoid05* Helicoid05::copy() const
{
    return new Helicoid05(*this);
}

// ----------   SET FUNCTIONS   ----------
// Function to set the screw parameters
void Helicoid05::set(Vec3D start, Mdouble length, Mdouble bladeRadius, Mdouble shaftRadius, Mdouble numberOfTurns, Mdouble omega, Mdouble thickness, bool rightHandeness)
{
    // check if the quantities are correct
    if (length <= 0.0) {std::cout << "\nERROR: INVALID LENGTH."; exit(0);}
    if (bladeRadius <= 0.0) {std::cout << "\nERROR: INVALID BLADE RADIUS."; exit(0);}
    if (shaftRadius >= bladeRadius) {std::cout << "\nERROR: INVALID SHAFT RADIUS."; exit(0);}
    if (numberOfTurns <= 0.0) {std::cout << "\nERROR: INVALID NUMBER OF SCREW TURNS."; exit(0);}
    if (thickness >= length) {std::cout << "\nERROR: INVALID THICKNESS."; exit(0);}
    
    // user defined variables
    start_ = start;
    l_ = length;
    rMax_ = bladeRadius;
    rMin_ = shaftRadius;
    n_ = numberOfTurns;
    omega_ = omega;
    thickness_ = thickness;
    offset_ = 0.0;
    rightHandeness_ = rightHandeness;
    // derived variables
    delta_ = 0.5 * thickness_;
    pitch_ = l_/n_;
    h_ = 0.5 * l_/(constants::pi * n_);
    
    std::cout << "\n\nScrew parameters set.\n";
    std::cout << "User defined variables:\n";
    std::cout << "Start : " << start_ << "\n";
    std::cout << "Length : " << l_ << "\n";
    std::cout << "Blade radius : " << rMax_ << "\n";
    std::cout << "Shaft radius : " << rMin_ << "\n";
    std::cout << "Number of turns : " << n_ << "\n";
    std::cout << "Angular velocity : " << omega_ << "\n";
    std::cout << "Thickness : " << thickness_ << "\n";
    std::cout << "Initial angular offset : " << offset_ << "\n";
    std::cout << "Right handeness : " << rightHandeness_ << "\n";
    std::cout << "Derived variables:\n";
    std::cout << "Half thickness : " << delta_ << "\n";
    std::cout << "Pitch length : " << pitch_ << "\n";
    std::cout << "Rescaled length : " << h_ << "\n\n";
}

// Function to change the screw radius
void Helicoid05::setRadius(Mdouble bladeRadius)
{
    rMax_ = bladeRadius;
}

// Function to change the thickness (and the half-thickness as well)
void Helicoid05::setThickness(Mdouble thickness)
{
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
}

// Function to change the screw rotation velocity
void Helicoid05::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// Function to change the screw angular offset
void Helicoid05::setOffset(Mdouble off)
{
    offset_ = off;
}

// ----------   GET FUNCTIONS   ------------------------------------------------------
// The contact detection algorithm
// All the quantities in cylindrical coordinates are referred to the screw axis
bool Helicoid05::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
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
    Mdouble cosEta = 1.0/sqrt(1.0 + pow(h_/rho, 2));
    
    // angular coordinate of the particle
    // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
    Mdouble xi;
    
    // trigonometric functions relative to the particle angle
    Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
    Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
    
    // oriented axial distance between the particle's centre and the blade centre
    Mdouble deltaZ;
    
    // normal to the blade at the particle position
    Vec3D normalVector;
    if (rightHandeness_) // right-handed thread
    {
        xi = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
        if (xi < 0.0) xi += 2.0*constants::pi;
        cosXi = (p.getPosition().X - start_.X)/rho;
        sinXi = (p.getPosition().Y - start_.Y)/rho;
        
        deltaZ = fmod((p.getPosition().Z - start_.Z) - h_*(xi + offset_) -(int)((p.getPosition().Z - start_.Z)/pitch_),pitch_);
        if (deltaZ > 0.5*pitch_) {deltaZ -= pitch_;}
        else if (deltaZ < -0.5*pitch_) {deltaZ += pitch_;}
        
        normalVector.X = h_*sinXi;
        normalVector.Y = -h_*cosXi;
        normalVector.Z = rho;
        // takes the right direction of the vector and normalizes
        normalVector *= - deltaZ;
        normalVector /= normalVector.getLength();
    }
    else // left-handed thread
    {
        xi = atan2(-(p.getPosition().Y - start_.Y), p.getPosition().X - start_.X);
        if (xi < 0.0) xi += 2.0*constants::pi;
        xi += 0.5*constants::pi;
        xi = fmod(xi, 2.0*constants::pi);
        
        cosXi = (p.getPosition().X - start_.X)/rho;
        sinXi = (p.getPosition().Y - start_.Y)/rho;
        
        deltaZ = fmod((p.getPosition().Z - start_.Z) - h_*(xi +0.5*constants::pi - offset_) -(int)((p.getPosition().Z - start_.Z)/pitch_),pitch_);
        if (deltaZ > 0.5*pitch_) {deltaZ -= pitch_;}
        else if (deltaZ < -0.5*pitch_) {deltaZ += pitch_;}
        
        normalVector.X = -h_*cosXi;
        normalVector.Y = h_*sinXi;
        normalVector.Z = rho;
        // takes the right direction of the vector and normalizes
        normalVector *= - deltaZ;
        normalVector /= normalVector.getLength();
    }
    
//    // takes the right direction of the vector and normalizes
//    normalVector *= - deltaZ;
//    normalVector /= normalVector.getLength();
    
    // radial vector at the particle position
    Vec3D radialVector;
    radialVector.X = cosXi;
    radialVector.Y = sinXi;
    radialVector.Z = 0.0;
    
    // component of the particle-blade_center distance normal to the blade surface
    Mdouble deltaN = fabs(deltaZ)*cosEta - delta_;
    
    // collision check
    if (rho >= rMax_) // check for collision with the blade edge
    {
        // radial component of the distance between the particle centre and the blade edge
        Mdouble deltaR = rho - rMax_;
        
        if (deltaN <= 0.0)
        {
            // region 1
            // the distance between the collision point and the particle centre
            distance = deltaR;
            
            // if the distance is negative prints an error message
            if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SIDE: OVERLAP > PARTICLE RADIUS\n";
            
            // the collision has no other component other than the radial vector
            normal_return = -radialVector;
            
//            std::cout << "Region 1 " << xi << "\t" << deltaZ << "\n";
            return true;
        }
        else
        {
            // distance between the particle centre and the and the edge
            distance = sqrt(pow(deltaN, 2) + pow(deltaR, 2));
            
            // if the particle-blade_edge distance is higher than the particle radius there is no collision
            if (distance > p.getWallInteractionRadius(this)) return false;
            
            // if the distance is negative prints an error message
            if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW EDGE: OVERLAP > PARTICLE RADIUS\n";
            
            // region 2
            // the normal to the edge through the particle centre
            // IMPORTANT: it needs to be normalized!
            normal_return = deltaN*normalVector - deltaR*radialVector;
            normal_return /= normal_return.getLength();
            
//            std::cout << "Region 2 " << xi << "\t" << deltaZ << "\n";
            return true;
        }
        
        // error message if some case has not been accounted for
        std::cout << "\nERROR: THE SCREW COLLISION FUNCTION DID NOT RETURN PROPERLY FOR rho > rMax\n";
        
        return false;
    }
    else if (rho >= rMin_ + p.getWallInteractionRadius(this)) // check for collision with the blade side only
    {
        // if the particle-blade_surface distance is higher than the collision threshold there is no collision
        if (deltaN > p.getWallInteractionRadius(this)) return false;
        
        // region 3
        // the distance between the collision point and the particle centre
        distance = deltaN;
        
        // if the distance is negative prints an error message
        if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SURFACE: OVERLAP > PARTICLE RADIUS\n";
        
        // the collision has no other component other than the normal to the blade
        normal_return = normalVector;
        
//        std::cout << "Region 3 " << xi << "\t" << deltaZ << "\n";
        return true;
    }
    else if (rho >= rMin_) // collision with the shaft
    {
        // radial distance between the particle centre and the shaft
        Mdouble deltaR = rho - rMin_;
        
        if (deltaN > p.getWallInteractionRadius(this)) // region 4
        {
            // the distance between the collision point and the particle centre
            distance = deltaR;
            
            // if the distance is negative prints an error message
            if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW SHAFT: OVERLAP > PARTICLE RADIUS\n";
            
            // the collision has no other component other than the radial vector
            normal_return = -radialVector;
            
            return true;
        }
        else // region 5
        {
            // distance between the particle centre and the and the point (rMin + rp, delta + rp)
            distance = sqrt(pow(deltaN - p.getWallInteractionRadius(this), 2) + pow(deltaR - p.getWallInteractionRadius(this), 2));
            
            // transform the former to be the actual distance of the particle from the contact point on the rounded corner
            distance = p.getWallInteractionRadius(this) - distance;
            
            // if the distance is negative prints an error message
            if (distance < 0.0) std::cout << "\nCOLLISION ERROR WITH THE SCREW CORNER: OVERLAP > PARTICLE RADIUS\n";
            
            // the normal to the edge through the particle centre
            // IMPORTANT: it needs to be normalized!
            normal_return = fabs(deltaN - p.getWallInteractionRadius(this))*normalVector - fabs(deltaR - p.getWallInteractionRadius(this))*radialVector;
            normal_return /= normal_return.getLength();
            
            return true;
        }
        
        // error message if some case has not been accounted for
        std::cout << "\nERROR: THE SCREW COLLISION FUNCTION DID NOT RETURN PROPERLY FOR rho >= rMin + p.getWallInteractionRadius(this)\n";
        
        return false;
    }
    else
    {
        
        // if the distance from the centre is smaller than the shaft radius prints an error
        std::cout << "\nERROR: rho < rMin\n";
        
        return false;
    }
    
    // error message if some case has not been accounted for
    std::cout << "\nERROR: THE SCREW COLLISION FUNCTION DID NOT RETURN PROPERLY\n";
    
    return false;
}

// Function that returns the current angular offset
Mdouble Helicoid05::getOffset()
{
    return offset_;
}

// Function that returns the screw handeness
bool Helicoid05::getRightHandeness()
{
    return rightHandeness_;
}

// ----------   ROTATION FUNCTIONS   ------------------------------------------------------
// Function to make the screw rotate by adding \omega*t to the angular offset
void Helicoid05::rotate(Mdouble dt)
{
    offset_ += omega_ * dt;
}

// Function to make the screw rotate by incrementing the angular offset
void Helicoid05::incrementOffset(Mdouble off)
{
    offset_ += off;
}

// ----------------------------------------------------------------------------------------
// Input stream
void Helicoid05::read(std::istream& is)
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
    >> dummy >> offset_
    >> dummy >> rightHandeness_;
}

// Old input stream (no idea what this is)
void Helicoid05::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> start_
    >> dummy >> l_
    >> dummy >> rMax_
    >> dummy >> rMin_
    >> dummy >> n_
    >> dummy >> omega_
    >> dummy >> offset_
    >> dummy >> rightHandeness_;
}

// Writer
void Helicoid05::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << start_
    << " Length " << l_
    << " Blade radius " << rMax_
    << " Shaft radius " << rMin_
    << " Turns " << n_
    << " Omega " << omega_
    << " Thickness " << thickness_
    << " Offset " << offset_
    << " Right handeness " << rightHandeness_;
}

// Returns the string "Helicoid05"
std::string Helicoid05::getName() const
{
    return "Helicoid05";
}

// Checks for the interaction between a particle p at a time timeStamp
// In case of interaction returns a pointer to the BaseInteraction happened between the Screw and the BaseParticle at time timeStamp
BaseInteraction* Helicoid05::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
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









