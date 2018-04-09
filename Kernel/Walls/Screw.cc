//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Screw.h"
#include "array"
#include "InteractionHandler.h"
#include "WallHandler.h"
#include "DPMBase.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
using mathsFunc::square;
/*!
 * \details Make a Screw which is centered in the origin, has a length of 1, one
 * revelation, a radius of 1, turns with 1 revelation per second, is infinitely thin
 * and starts at its normal initial point.
 */
Screw::Screw()
{
    start_.setZero();
    l_ = 1.0;
    maxR_ = 1.0;
    n_ = 1.0;
    omega_ = 1.0;
    offset_ = 0.0;
    thickness_ = 0.0;
    setOrientationViaNormal({0,0,1});//default screw is in z-direction
    logger(DEBUG, "Screw() constructor finished.");              
}

/*!
 * \param[in] other The Screw that has to be copied.
 */
Screw::Screw(const Screw& other)
    : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    maxR_ = other.maxR_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    offset_ = other.offset_;
    logger(DEBUG, "Screw(const Screw&) copy constructor finished.");
}

/*!
 * \param[in] start A Vec3D which denotes the centre of the lower end of the Screw.
 * \param[in] l The length of the Screw, must be positive.
 * \param[in] r The radius of the Screw, must be positive.
 * \param[in] n The number of revelations of the Screw, must be positive.
 * \param[in] omega The rotation speed of the Screw in rev/s.
 * \param[in] thickness The thickness of the Screw, must be non-negative.
 * \details Make a Screw by assigning all input parameters to the data-members of
 * this class, and setting the offset_ to 0.
 */
Screw::Screw(Vec3D start, Mdouble l, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = l;
    maxR_ = r;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    offset_ = 0.0;
    logger(DEBUG, "Screw(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");           
}

Screw::~Screw()
{
    logger(DEBUG, "~Screw() finished, destroyed the Screw.");
}

/*!
 * \return A pointer to a copy of this Screw.
 */
Screw* Screw::copy() const
{
    return new Screw(*this);
}

bool Screw::getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const
{
    //transform coordinates into position-orientation frame
    Vec3D position = P.getPosition() - getPosition();
    getOrientation().rotateBack(position);
    ///\todo do this for all walls
    BaseSpecies* s = getHandler()->getDPMBase()->speciesHandler.getMixedObject(P.getSpecies(),getSpecies());
    if (getDistanceAndNormalLabCoordinates(position, P.getRadius()+s->getInteractionDistance(), distance, normal_return)) {
        getOrientation().rotate(normal_return);
        return true;
    } else {
        return false;
    }
}

/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between 
 * a given BaseParticle and this Screw. If there is a collision, this
 * function also computes the distance between the BaseParticle and Screw
 * and the normal of the IntersectionOfWalls at the intersection point. 
 * \todo Make this function readable and explain the steps in the details.
 */
bool Screw::getDistanceAndNormalLabCoordinates(Vec3D position, Mdouble wallInteractionRadius, Mdouble& distance, Vec3D& normal_return) const
{
    Mdouble RSquared = square(position.Y - start_.Y) + square(position.Z - start_.Z);
    Mdouble X = position.X - start_.X;
    //first do a simple check if particle is within the cylindrical hull of the screw
    if (RSquared > square(maxR_ + wallInteractionRadius + thickness_)
        || X > l_ + wallInteractionRadius + thickness_
        || X < - wallInteractionRadius - thickness_)
    {
        //std::cout << "failed " << position << std::endl;
        return false;
    }

    //else:
    Mdouble R = sqrt(RSquared);
    Mdouble A = atan2(position.Z - start_.Z, position.Y - start_.Y);
    
    //after subtracting the start position and transforming the particle position from (XYZ) into (RAZ) 
    //coordinates, we compute the distance to the wall at, located at (r,a,z)=(r,2*pi*(offset+N*q+k/2),q*L), 0<q<1.

    //To find the contact point we have to minimize (with respect to r and q)
    //distance^2=(x-x0-r*cos(2*pi*(offset+N*q)))^2+(y-y0-r*sin(2*pi*(offset+N*q)))^2+(z-z0-q*L)^2
    //Using polar coordinates (i.e. x-x0=R*cos(A), y-y0=R*sin(A) and Z=z-z0)
    //distance^2=R^2+r^2-2*R*r*cos(A-2*pi*(offset+N*q))+(Z-q*L)^2
    
    //Assumption: d(distance)/dr=0 at minDistance (should there be also a q-derivative?)
    //Differentiate with respect to r and solve for zero:
    //0=2*r-2*R*cos(A-2*pi*(offset+N*q)
    //r=R*cos(A-2*pi*(offset+N*q))
    
    //Substitue back
    //distance^2=R^2+R^2*cos^2(A-2*pi*(offset+N*q))-2*R^2*cos^2(A-2*pi*(offset+N*q))+(Z-q*L)^2
    //distance^2=R^2*sin^2(A-2*pi*(offset+N*q))+(Z-q*L)^2
    
    //So we have to minimize:
    //distance^2=R^2*sin^2(A-2*pi*(offset+N*q))^2 + (Z-q*L)^2 = f(q)
    //f'(q)=(-2*pi*N)*R^2*sin(2*A-4*pi*(offset+N*q)) + 2*L*(Z-q*L)    (D[Sin[x]^2,x]=Sin[2x])
    //f''(q)=(4*pi^2*N^2)*R^2*cos(2*A-4*pi*(offset+N*q)) - 2*L*L
    //For this we use the Euler algoritm
    
    Mdouble q; //Current guess
    Mdouble dd; //Derivative at current guess
    Mdouble ddd; //Second derivative at current guess
    Mdouble q0 = X / l_; //assume closest point q0 is at same z-location as the particle
            
    //Set initial q to the closest position on the screw with the same angle as A
    //The initial guess will be in the minimum of the sin closest to q0
    //Minima of the sin are at
    //A-2*pi*(offset+N*q)=k*pi (k=integer)
    //q=A/(2*pi*N)-k/(2*N)-offset/N (k=integer)
    
    Mdouble k = round(A / constants::pi - 2.0 * (offset_ + n_ * q0)); // k: |A-a(q0)-k*pi|=min
    q = A / (2.0 * constants::pi * n_) - k / (2.0 * n_) - offset_ / n_; // q: a(q)=A
    
    //Now apply Newton's method
    do
    {
        Mdouble arg = 2.0 * A - 4.0 * constants::pi * (n_ * q + offset_);
        dd = -2.0 * constants::pi * n_ * RSquared * sin(arg) - 2.0 * l_ * (X - q * l_);
        ddd = 8.0 * constants::sqr_pi * n_ * n_ * RSquared * cos(arg) + 2.0 * l_ * l_;
        q -= dd / ddd;
    } while (fabs(dd / ddd) > 1e-14);
    
    //Calculate r
    Mdouble r = R * cos(2.0 * constants::pi * (offset_ + n_ * q) - A);
    
    //Check if the location is actually on the screw:
    //First posibility is that the radius is too large:
    if (fabs(r) > maxR_) //Left boundary of the coil
    {
        r = mathsFunc::sign(r) * maxR_;
        unsigned int steps = 0;
        //This case reduces to the coil problem
        do
        {
            dd = -4.0 * R * r * constants::pi * n_ * sin(A - 2.0 * constants::pi * (n_ * q + offset_)) - 2.0 * l_ * (X - q * l_);
            ddd = 8.0 * R * r * constants::sqr_pi * n_ * n_ * cos(A - 2.0 * constants::pi * (n_ * q + offset_)) + 2.0 * l_ * l_;
            q -= dd / ddd;
            steps++;
        } while (fabs(dd / ddd) > 1e-14);
    }
    //Second possibility is that it occured before the start of after the end
    if (q < 0)
    {
        q = 0;
        r = R * cos(A - 2.0 * constants::pi * (offset_ + q * n_));
        if (fabs(r) > maxR_)
        {
            r = mathsFunc::sign(r) * maxR_;
        }
    }
    else if (q > 1)
    {
        q = 1;
        r = R * cos(A - 2.0 * constants::pi * (offset_ + q * n_));
        if (fabs(r) > maxR_)
        {
            r = mathsFunc::sign(r) * maxR_;
        }
    }
    
    Mdouble distanceSquared = R * R * pow(sin(A - 2 * constants::pi * (offset_ + n_ * q)), 2) + pow(X - q * l_, 2);
    //If distance is too large there is no contact
    if (distanceSquared >= square(wallInteractionRadius + thickness_))
    {
        return false;
    }
    
    Vec3D ContactPoint;
    distance = sqrt(distanceSquared) - thickness_;
    ContactPoint.Y = start_.Y + r * cos(2.0 * constants::pi * (offset_ + n_ * q));
    ContactPoint.Z = start_.Z + r * sin(2.0 * constants::pi * (offset_ + n_ * q));
    ContactPoint.X = start_.X + q * l_;
    normal_return = (ContactPoint - position);
    normal_return.normalize();
    return true;
}

/*!
 * \param[in] dt The time for which the Screw has to be turned.
 */
void Screw::move_time(Mdouble dt)
{
    //offset_ += omega_ * dt;
}

/**
 * \todo the move and rotate functions should only pass the time step, as teh velocity can be accessed directly by the object
 * \param angularVelocityDt
 */
void Screw::rotate(const Vec3D &angularVelocityDt)
{
    BaseInteractable::rotate(angularVelocityDt);
    offset_ += omega_ * getHandler()->getDPMBase()->getTimeStep();
}


/*!
 * \param[in,out] is Input stream from which the Screw must be read.
 */
void Screw::read(std::istream& is)
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

/*!
 * \param[in,out] is Input stream from which the Screw must be read.
 * \details Read the Screw in old style, please note that the thickness is not 
 * read in this function, so it has either to be set manually or it is 0.0 from
 * the default constructor.
 */
void Screw::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> start_
            >> dummy >> l_
            >> dummy >> maxR_
            >> dummy >> n_
            >> dummy >> omega_
            >> dummy >> offset_;
}

/*!
 * \param[in,out] os Output stream to which the Screw must be written.
 */
void Screw::write(std::ostream& os) const
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

/*!
 * \return The string "Screw".
 */
std::string Screw::getName() const
{
    return "Screw";
}

void Screw::writeVTK (VTKContainer &vtk) const {
    unsigned nr = 10;
    unsigned nz = 60*abs(n_);

    unsigned nPoints = vtk.points.size();
    vtk.points.reserve(nPoints+nr*nz);
    Vec3D contactPoint;
    for (unsigned iz=0; iz<nz; iz++) {
        for (unsigned ir=0; ir<nr; ir++) {
            double q = (double)iz/nz;
            double r = (double)ir/nr*maxR_;
            contactPoint.Y = start_.Y - r * cos(2.0 * constants::pi * (offset_ + n_ * q));
            contactPoint.Z = start_.Z - r * sin(2.0 * constants::pi * (offset_ + n_ * q));
            contactPoint.X = start_.X + q * l_;
            getOrientation().rotate(contactPoint);
            contactPoint += getPosition();
            vtk.points.push_back(contactPoint);
        }
    }

    unsigned nCells = vtk.triangleStrips.size();
    vtk.triangleStrips.reserve(nCells+(nz-1));
    for (unsigned iz=0; iz<nz-1; iz++) {
        std::vector<double> cell;
        cell.reserve(2*nr);
        for (unsigned ir=0; ir<nr; ir++) {
            cell.push_back(nPoints+ir+iz*nr);
            cell.push_back(nPoints+ir+(iz+1)*nr);
        }
        vtk.triangleStrips.push_back(cell);
    }
}

void Screw::writeVTK (std::string filename) const
{
    VTKContainer vtk;
    writeVTK(vtk);

    std::stringstream file;
    file << "# vtk DataFile Version 2.0\n"
     << getName() << "\n"
     "ASCII\n"
     "DATASET UNSTRUCTURED_GRID\n"
     "POINTS " << vtk.points.size() << " double\n";
    for (const auto& vertex : vtk.points)
        file << vertex << '\n';
    file << "\nCELLS " << vtk.triangleStrips.size() << ' ' << 4*vtk.triangleStrips.size() << "\n";
    for (const auto& face : vtk.triangleStrips)
        file << "3 " << face[0] << ' ' << face[1] << ' ' << face[2] << '\n';
    file << "\nCELL_TYPES " << vtk.triangleStrips.size() << "\n";
    for (const auto& face : vtk.triangleStrips)
        file << "5\n";
    helpers::writeToFile(filename,file.str());
}
