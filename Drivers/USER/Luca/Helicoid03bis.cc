#include "Helicoid03bis.h"
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
Helicoid03bis::Helicoid03bis()
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
    
    logger(DEBUG, "Helicoid03bis() constructor finished.");
}

// Copy constructor
Helicoid03bis::Helicoid03bis(const Helicoid03bis& other)
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
    
    logger(DEBUG, "Helicoid03bis(const Helicoid03bis&) copy constructor finished.");
}

// The good constructor
// Parameters in: the initial starting point of the screw axis, the length, the radius, the number of revelations,
// the angular velocity, the thickness of the blade.
Helicoid03bis::Helicoid03bis(Vec3D start, Mdouble l, Mdouble R, Mdouble n, Mdouble omega, Mdouble thickness)
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
    
    logger(DEBUG, "Helicoid03bis(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

// Destructor
Helicoid03bis::~Helicoid03bis()
{
    logger(DEBUG, "~Helicoid03bis() finished, destroyed the Helicoid.");
}

// A pointer to the copy of the helicoid
Helicoid03bis* Helicoid03bis::copy() const
{
    return new Helicoid03bis(*this);
}

// Function to set the Helicoid parameters
void Helicoid03bis::set(Vec3D start, Mdouble length, Mdouble bladeRadius, Mdouble numberOfRevelations, Mdouble omega, Mdouble thickness)
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
void Helicoid03bis::setRadius(Mdouble bladeRadius)
{
    maxR_ = bladeRadius;
}

// Function to change the thickness (and the half-thickness as well)
void Helicoid03bis::setThickness(Mdouble thickness)
{
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
}

// Function to change the helicoid rotation velocity
void Helicoid03bis::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// Function to make the screw rotate by adding \omega*t to the angular offset
void Helicoid03bis::move_time(Mdouble dt)
{
    offset_ += omega_ * dt;
}

Mdouble Helicoid03bis::getOffset()
{
    return offset_;
}

void Helicoid03bis::setOffset(Mdouble off)
{
    offset_ = off;
}

// Function to make the screw rotate by incrementing the angular offset (should substitute the move_time function)
void Helicoid03bis::incrementOffset(Mdouble off)
{
    offset_ += off;
}

// Helicoid Input stream
void Helicoid03bis::read(std::istream& is)
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
void Helicoid03bis::oldRead(std::istream& is)
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
void Helicoid03bis::write(std::ostream& os) const
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
std::string Helicoid03bis::getName() const
{
    return "Helicoid";
}

// The contact detection algorithm
// Takes as input a pointer to a particle and evaluates if there is a collision.
// If so updates the distance between the particle centre and the collision point, the collision normal and returns true.
// Returns false if there is no collision.
bool Helicoid03bis::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // components of particle distance from the blade
    Mdouble dZ, dR, dN;
    
    // z coordinate check
    if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius() || p.getPosition().Z < start_.Z - p.getWallInteractionRadius()) return false;
    
    // radial distance of the particle from the helicoid axis
    Mdouble rho = sqrt(pow(p.getPosition().X - start_.X, 2) + pow(p.getPosition().Y - start_.Y, 2));
    
    // r coordinate check
    if (rho > maxR_ + p.getWallInteractionRadius()) return false;
    
    // cosine of the helix angle
    Mdouble cosEta = 1.0/sqrt(1.0+pow(h_/rho,2));
    
    // angular coordinate of the particle in the cylindrical frame of reference centered on the helicoid axis
    // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
    Mdouble xi = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
    if (xi < 0.0) xi += 2.0*constants::pi;
    
    // the "oriented" axial distance between the particle's centre and the blade centre
    Mdouble deltaZ = fmod((p.getPosition().Z - start_.Z) - h_*(xi + offset_) - (int)((p.getPosition().Z - start_.Z)/pitch_),pitch_);
    if (deltaZ > 0.5*pitch_) {deltaZ -= pitch_;}
    else if (deltaZ < -0.5*pitch_) {deltaZ += pitch_;}
    
    // n coordinate check
    if (fabs(deltaZ)*cosEta > p.getWallInteractionRadius() + delta_) return false;
    
    // sine of the helix angle
    Mdouble sinEta = 1.0/sqrt(1.0+pow(rho/h_,2));
    
    // the distance from the blade in r, n, z direction are assigned
    dR = rho - maxR_;
    dZ = fabs(p.getPosition().Z - 0.5*(l_ + start_.Z)) - 0.5*(l_ - start_.Z);
    dN = fabs(deltaZ)*cosEta - delta_;
    
    // trigonometric functions relative to the particle position
    Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
    Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
    
    // collision cases
    if (dZ <= -p.getWallInteractionRadius()*cosEta) // (-,-,0) collision
    {
        if (dR <= 0.0) // (0,n,0) collision
        {
            // n vector
            Vec3D nVector;
            nVector.X = sinXi*sinEta;
            nVector.Y = -cosXi*sinEta;
            nVector.Z = cosEta;
            
            // the collision distance
            distance = dN;
            
            // gets the right n vector direction
            nVector *= -deltaZ;
            nVector /= nVector.getLength();
            
            // the collision vector
            normal_return = nVector;
            
            return true;
        }
        else // (r,-,0) collision
        {
            // r vector
            Vec3D rVector;
            rVector.X = -cosXi;
            rVector.Y = -sinXi;
            rVector.Z = 0.0;
            
            if (dN >= 0.0) // (r,n,0) collision
            {
                // n vector
                Vec3D nVector;
                nVector.X = sinXi*sinEta;
                nVector.Y = -cosXi*sinEta;
                nVector.Z = cosEta;
                
                // the collision distance
                distance = sqrt(pow(dN,2.0) + pow(dR,2.0));
                
                if (distance > p.getWallInteractionRadius()) return false;
                
                // gets the right n vector direction
                nVector *= -deltaZ;
                nVector /= nVector.getLength();
                
                // the collision vector
                normal_return = (dN*nVector + dR*rVector)/distance;
                
                return true;
            }
            else // (r,0,0) collision
            {
                // the collision distance
                distance = dR;
                
                // the collision vector
                normal_return = rVector;
                
                return true;
            }
        }
        
        std::cout << "\nERROR: THE COLLISION FUNCTION DID NOT RETURN PROPERLY!\n";
    }
    else
    {
        // z vector
        Vec3D zVector;
        zVector.X = 0.0;
        zVector.Y = 0.0;
        zVector.Z = -1.0;
        
        if (dN <= 0.0) // (-,0,-) collision
        {
            if (dR >= 0.0) // (r,0,-) collision
            {
                // r vector
                Vec3D rVector;
                rVector.X = -cosXi;
                rVector.Y = -sinXi;
                rVector.Z = 0.0;
                
                if (dZ <= 0.0) // (r,0,0) collision
                {
                    // the collision distance
                    distance = dR;
                    
                    // the collision vector
                    normal_return = rVector;
                    
                    return true;
                }
                else // (r,0,z) collision
                {
                    // the collision distance
                    distance = sqrt(pow(dZ,2.0) + pow(dR,2.0));
                    
                    if (distance > p.getWallInteractionRadius()) return false;
                    
                    // gets the right z vector direction
                    zVector *= p.getPosition().Z - 0.5*(l_ + start_.Z);
                    zVector /= zVector.getLength();
                    
                    // the collision vector
                    normal_return = (dR*rVector + dZ*zVector)/distance;
                    
                    return true;
                }
            }
            else // (0,0,z) collision
            {
                // the collision distance
                distance = dZ;
                
                // gets the right z vector direction
                zVector *= p.getPosition().Z - 0.5*(l_ + start_.Z);
                zVector /= zVector.getLength();
                
                // the collision vector
                normal_return = zVector;
                
                return true;
            }
        }
        else // (-,n,z) collision
        {
            // phi vector
            Vec3D phiVector;
            phiVector.X = sinXi;
            phiVector.Y = -cosXi;
            phiVector.Z = 0.0;
            
            // polar angle related variable
            double dPhi = dN/sinXi;
            double zPhiDistance = sqrt(pow(dPhi,2.0) + pow(dZ,2.0));
            
            if (zPhiDistance > p.getWallInteractionRadius()) return false;
            
            // gets the right phi vector direction
            phiVector *= -deltaZ;
            phiVector /= phiVector.getLength();
            
            // gets the right z vector direction
            zVector *= p.getPosition().Z - 0.5*(l_ + start_.Z);
            zVector /= zVector.getLength();
            
            if (deltaZ*(p.getPosition().Z - 0.5*(l_ + start_.Z)) > 0.0) // edgy collision
            {
                if (dR <= 0.0) // (0,n,z) collision
                {
                    // the collision vector
                    normal_return = (dPhi*phiVector + dZ*zVector)/zPhiDistance;
                    
                    return true;
                }
                else // (r,n,z) collision
                {
                    // r vector
                    Vec3D rVector;
                    rVector.X = -cosXi;
                    rVector.Y = -sinXi;
                    rVector.Z = 0.0;
                    
                    // the collision distance
                    distance = sqrt(pow(dR,2.0) + pow(zPhiDistance,2.0));
                    
                    if (distance > p.getWallInteractionRadius()) return false;
                    
                    // the collision vector
                    normal_return = (dPhi*phiVector + dZ*zVector + dR*rVector)/distance;
                    
                    return true;
                }
            }
            else // flatty collision
            {
                if (dR <= 0.0) // (0,n,z) collision
                {
                    if (dZ > p.getWallInteractionRadius()*cosEta) // actual collision with the edge
                    {
                        // the collision vector
                        normal_return = (dPhi*phiVector + dZ*zVector)/zPhiDistance;
                        
                        return true;
                    }
                    else // still collision with the blade only
                    {
                        // the collision vector
                        normal_return = sinEta*phiVector + cosEta*zVector;
                        
                        return true;
                    }
                }
                else // (r,n,z) collision
                {
                    // r vector
                    Vec3D rVector;
                    rVector.X = -cosXi;
                    rVector.Y = -sinXi;
                    rVector.Z = 0.0;
                    
                    if (dZ > p.getWallInteractionRadius()*cosEta) // actual collision with the edge
                    {
                        // the collision distance
                        distance = sqrt(pow(dR,2.0) + pow(zPhiDistance,2.0));
                        
                        if (distance > p.getWallInteractionRadius()) return false;
                        
                        // the collision vector
                        normal_return = (dPhi*phiVector + dZ*zVector + dR*rVector)/distance;
                        
                        return true;
                    }
                    else // still collision with the blade only
                    {
                        // n vector
                        Vec3D nVector;
                        nVector.X = sinXi*sinEta;
                        nVector.Y = -cosXi*sinEta;
                        nVector.Z = cosEta;
                        
                        // the collision distance
                        distance = sqrt(pow(dN,2.0) + pow(dR,2.0));
                        
                        if (distance > p.getWallInteractionRadius()) return false;
                        
                        // gets the right n vector direction
                        nVector *= -deltaZ;
                        nVector /= nVector.getLength();
                        
                        // the collision vector
                        normal_return = (dN*nVector + dR*rVector)/distance;

                        return true;
                    }
                }
            }
        }
        
        std::cout << "\nERROR: THE COLLISION FUNCTION DID NOT RETURN PROPERLY!\n";
    }
    
    
//    
//    
//    // axial distance (non zero only if there is an axial component in the collision)
//    if (p.getPosition().Z > l_ + start_.Z) {dZ = p.getPosition().Z - l_ - start_.Z;}
//    else if (p.getPosition().Z < start_.Z) {dZ = p.getPosition().Z - start_.Z;}
//    else {dZ = 0.0;}
//    
//    // radial distance (non zero only if there is a radial component in the collision)
//    if (rho > maxR_) {dR = rho - maxR_;}
//    else {dR = 0.0;}
//    
//    // normal distance (non zero only if there is an axial component in the collision)
//    if (fabs(deltaZ)*cosEta > delta_) {dN = fabs(deltaZ)*cosEta - delta_;}
////    else if (deltaZ*cosEta < -delta_) {dN = deltaZ*cosEta + delta_;}
//    else {dN = 0.0;}
//    
//
//    
//    if (dZ) // (-,-,z) collision
//    {
//        // axial vector
//        Vec3D zVector;
//        zVector.X = 0.0;
//        zVector.Y = 0.0;
//        zVector.Z = -1.0;
//        
//        if (dR) // (r,-,z) collision
//        {
//            // trigonometric functions relative to the particle position
//            Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
//            Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
//            
//            // radial vector
//            Vec3D rVector;
//            rVector.X = -cosXi;
//            rVector.Y = -sinXi;
//            rVector.Z = 0.0;
//            
//            if (dN) // (r,n,z) collision
//            {
//                // normal vector
//                Vec3D nVector;
//                nVector.X = h_*sinXi;
//                nVector.Y = -h_*cosXi;
//                nVector.Z = rho;
//                
//                // distance from the contact point
//                distance = sqrt(pow(dZ,2) + pow(dR,2) + pow(dN,2));
//                if (distance > p.getWallInteractionRadius()) return false;
//                
//                // takes the correct orientation of zVector and normalizes
//                zVector *= -dZ;
//                zVector /= zVector.getLength();
//                
//                // takes the correct orientation of nVector and normalizes
//                nVector *= -deltaZ;
//                nVector /= nVector.getLength();
//                
//                // the collision vector
//                normal_return = dZ*zVector + dR*rVector + dN*nVector;
//                normal_return /= normal_return.getLength();
//                
////                std::cout << "(r,n,z)\n";
//                
//                return true;
//            }
//            else // (r,0,z) collision
//            {
//                // distance from the contact point
//                distance = sqrt(pow(dZ,2) + pow(dR,2));
//                if (distance > p.getWallInteractionRadius()) return false;
//                
//                // takes the correct orientation of zVector and normalizes
//                zVector *= -dZ;
//                zVector /= zVector.getLength();
//                
//                // the collision vector
//                normal_return = dZ*zVector + dR*rVector;
//                normal_return /= normal_return.getLength();
//                
////                std::cout << "(r,0,z)\n";
//                
//                return true;
//            }
//        }
//        else // (0,-,z) collision
//        {
//            if (dN) // (0,n,z) collision
//            {
//                // trigonometric functions relative to the particle position
//                Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
//                Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
//                
//                // normal vector
//                Vec3D nVector;
//                nVector.X = h_*sinXi;
//                nVector.Y = -h_*cosXi;
//                nVector.Z = rho;
//                
//                // distance from the contact point
//                distance = sqrt(pow(dZ,2) + pow(dN,2));
//                if (distance > p.getWallInteractionRadius()) return false;
//                
//                // takes the correct orientation of zVector and normalizes
//                zVector *= -dZ;
//                zVector /= zVector.getLength();
//                
//                // takes the correct orientation of nVector and normalizes
//                nVector *= -deltaZ;
//                nVector /= nVector.getLength();
//                
//                // the collision vector
//                normal_return = dZ*zVector + dN*nVector;
//                normal_return /= normal_return.getLength();
//                
////                std::cout << "(0,n,z)\n";
//                
//                return true;
//            }
//            else // (0,0,z) collision
//            {
//                // distance from the contact point
//                distance = fabs(dZ);
//                                
//                // takes the correct orientation of zVector and normalizes
//                zVector *= -dZ;
//                zVector /= zVector.getLength();
//                
//                // the collision vector
//                normal_return = zVector;
//                
////                std::cout << "(0,0,z)\n";
////                std::cout << distance << "\n";
//                
//                return true;
//            }
//        }
//    }
//    else // (-,-,0) collision
//    {
//        if (dR) // (r,-,0) collision
//        {
//            // trigonometric functions relative to the particle position
//            Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
//            Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
//            
//            // radial vector
//            Vec3D rVector;
//            rVector.X = -cosXi;
//            rVector.Y = -sinXi;
//            rVector.Z = 0.0;
//            
//            if (dN) // (r,n,0) collision
//            {
//                // normal vector
//                Vec3D nVector;
//                nVector.X = h_*sinXi;
//                nVector.Y = -h_*cosXi;
//                nVector.Z = rho;
//                
//                // distance from the contact point
//                distance = sqrt(pow(dR,2) + pow(dN,2));
//                if (distance > p.getWallInteractionRadius()) return false;
//                
//                // takes the correct orientation of nVector and normalizes
//                nVector *= -deltaZ;
//                nVector /= nVector.getLength();
//                
//                // the collision vector
//                normal_return = dR*rVector + dN*nVector;
//                normal_return /= normal_return.getLength();
//                
//                return true;
//            }
//            else // (r,0,0) collision
//            {
//                // distance from the contact point
//                distance = dR;
//                
//                // the collision vector
//                normal_return = rVector;
//                
//                return true;
//            }
//        }
//        else // (0,n,0) collision
//        {
//            // trigonometric functions relative to the particle position
//            Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
//            Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
//            
//            // normal vector
//            Vec3D nVector;
//            nVector.X = h_*sinXi;
//            nVector.Y = -h_*cosXi;
//            nVector.Z = rho;
//            
//            // distance from the contact point
//            distance = dN;
//            
//            // takes the correct orientation of nVector and normalizes
//            nVector *= -deltaZ;
//            nVector /= nVector.getLength();
//            
//            // the collision vector
//            normal_return = nVector;
//            
//            return true;
//        }
//    }
}

// Checks for the interaction between a particle p at a time timeStamp.
// In case of interaction returns a pointer to the BaseInteraction happened between the Helicoid and the
// BaseParticle at time timeStamp
std::vector<BaseInteraction *> Helicoid03bis::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
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









