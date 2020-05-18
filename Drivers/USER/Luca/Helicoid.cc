//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include "Helicoid.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// Last update 24.08.15
// TODO - I'm not using the offset_ variable. Check to remove all teh calls and so on.
// TODO - The actual computation of teh normal is kinda retarded since I compute it the otehr way around in teh notes.

/*!
 * \details Make a Helicoid which is centered in the origin, has a length of 1, one
 * revelation, a radius of 1, turns with 1 revelation per second, is infinitely thin
 * and starts at its normal initial point.
 */
Helicoid::Helicoid()
{
    start_.setZero();
    l_ = 1.0;
    maxR_ = 1.0;
    n_ = 1.0;
    omega_ = 1.0;
    offset_ = 0.0;
    thickness_ = 0.0;
    logger(DEBUG, "Helicoid() constructor finished.");              
}

/*!
 * \param[in] other The Helicoid that has to be copied.
 */
Helicoid::Helicoid(const Helicoid& other)
    : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    maxR_ = other.maxR_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    offset_ = other.offset_;
    logger(DEBUG, "Helicoid(const Helicoid&) copy constructor finished.");
}


/*!
 * \param[in] start A Vec3D which denotes the centre of the lower end of the Helicoid.
 * \param[in] l The length of the Helicoid, must be positive.
 * \param[in] r The radius of the Helicoid, must be positive.
 * \param[in] n The number of revelations of the Helicoid, must be positive.
 * \param[in] omega The rotation speed of the Helicoid in rev/s.
 * \param[in] thickness The thickness of the Helicoid, must be non-negative.
 * \details Make a Helicoid by assigning all input parameters to the data-members of
 * this class, and setting the offset_ to 0.
 */
Helicoid::Helicoid(Vec3D start, Mdouble l, Mdouble r, Mdouble n, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = l;
    maxR_ = r;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    offset_ = 0.0;
    logger(DEBUG, "Helicoid(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

Helicoid::~Helicoid()
{
    logger(DEBUG, "~Helicoid() finished, destroyed the Helicoid.");
}

/*!
 * \return A pointer to a copy of this Helicoid.
 */
Helicoid* Helicoid::copy() const
{
    return new Helicoid(*this);
}

// HACK - BEGIN
// need the proper constructor, like in the coil
void Helicoid::set(Vec3D start, Mdouble length, Mdouble radius, Mdouble numberOfRevelations, Mdouble omega, Mdouble thickness)
{
    start_ = start;
    l_ = length;
    maxR_ = radius;
    n_ = numberOfRevelations;
    omega_ = omega;
    thickness_ = thickness;
    offset_ = 0.0;
    
    delta_ = 0.5 * thickness_;
    h_ = l_/(2.0 * constants::pi * n_);
    cosEta_ = cos(atan2(h_,2.0 * constants::pi * maxR_));
    
    std::cout << "\n\n Helicoid parameters set.\n";
    std::cout << "Start : " << start_ << "\n";
    std::cout << "Length : " << l_ << "\n";
    std::cout << "Radius : " << maxR_ << "\n";
    std::cout << "Number of turns : " << n_ << "\n";
    std::cout << "Angular velocity : " << omega_ << "\n";
    std::cout << "Thickness : " << thickness_ << "\n";
    std::cout << "Half thickness : " << delta_ << "\n";
    std::cout << "Rescaled length : " << h_ << "\n";
    std::cout << "Cosine of helicoid angle : " << cosEta_ << "\n\n";
}

// so I can change the radius on the run
// the quantities depending on the radius should be re-evaluated as well
void Helicoid::setRadius(Mdouble radius)
{
    maxR_ = radius;
    cosEta_ = cos(atan2(h_,2.0 * constants::pi * maxR_));
}

// so I can change the thickness on the run
// the quantities depending on the thickness should be re-evaluated as well
void Helicoid::setThickness(Mdouble thickness)
{
    thickness_ = thickness;
    delta_ = 0.5 * thickness_;
}

// so I can change the omega on the run
void Helicoid::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// this stuff is for checking the collision on the run
// std::cout << p.getId() << "\t" << fabs(phi)/tau << "\t" << signPhi << "\t" << epsilon/(p.getWallInteractionRadius(this)+delta_) << "\t" << (distance+delta_)/(p.getWallInteractionRadius(this)+delta_) << "\t" << p.getPosition().X << "\t" << p.getPosition().Y << "\t" << p.getPosition().Z << "\t" << ContactPoint.X << "\t" << ContactPoint.Y << "\t" << ContactPoint.Z << "\t" << normal_return.X << "\t" << normal_return.Y << "\t" << normal_return.Z << "\n";

// the contact detection algorithm
bool Helicoid::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial distance of the particle from the helicoid axis
    Mdouble rho2 = pow(p.getPosition().X - start_.X, 2) + pow(p.getPosition().Y - start_.Y, 2);
    
    // if the particle is outside of the cylinder that contains the helicoid returns false
    // check for the radius
    if (rho2 > pow(maxR_ + p.getWallInteractionRadius(this), 2)) return false;
    // check for the height
    if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius(this) + delta_) return false;
    if (p.getPosition().Z < start_.Z - p.getWallInteractionRadius(this) - delta_) return false;
    
    // radial distance of the particle from the helicoid axis
    Mdouble rho = sqrt(rho2);
    
    // angular coordinate of the particle in the cylindrical frame of reference centered on the helicoid axis
    // IMPORTANT: this angle needs to be defined in the interval [-p1, +pi] radians!
    Mdouble xi = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
    
    // angle of the helicoid at the particle's centre
    Mdouble theta = fmod((p.getPosition().Z - start_.Z)/h_ - offset_, 2.0 * constants::pi); // <- HERE
    
    // angular distance between the helicoid and the particle's centre
    Mdouble phi = 0.0;
    if (fabs(theta - xi) <= constants::pi) phi = theta - xi;
    else phi = theta - xi - 2.0 * constants::pi;
    
    // angular distance threshold for the collision
    Mdouble tau = (p.getWallInteractionRadius(this) + delta_)/(h_ * cosEta_);
    
    // if the absolute value of phi is bigger than tau there is no collision
    if (fabs(phi) > tau) return false;
    
    // if not we need to evaluate the following stuff
    // overlap between particle and helicoid
    Mdouble epsilon = h_ * cosEta_ * (tau - fabs(phi));
    
    // distance between the contact point and the particle's centre
    distance = (p.getWallInteractionRadius(this) - epsilon)/cosEta_;
    
    // sign of phi
    Mdouble signPhi = 0.0;
    if (std::signbit(phi)) signPhi = -1.0;
    else signPhi = 1.0;
    
    // normalization factor for the normal
    Mdouble normFactor = 1/sqrt(pow(h_,2) + rho2);
    
    // contact point computation
    Vec3D ContactPoint;
    ContactPoint.X = p.getPosition().X - signPhi * distance * normFactor * h_ * sin(xi);
    ContactPoint.Y = p.getPosition().Y + signPhi * distance * normFactor * h_ * cos(xi);
    ContactPoint.Z = p.getPosition().Z - signPhi * distance * normFactor * rho;
    
    // computation of the normal and normalization
    normal_return = ContactPoint - p.getPosition();
    normal_return /= normal_return.getLength();
    
    return true;
}

// HACK - END

/*!
 * \param[in] dt The time for which the Helicoid has to be turned.
 */
void Helicoid::move_time(Mdouble dt)
{
    offset_ += omega_ * dt;
}

/*!
 * \param[in,out] is Input stream from which the Helicoid must be read.
 */
void Helicoid::read(std::istream& is)
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
 * \param[in,out] os Output stream to which the Helicoid must be written.
 */
void Helicoid::write(std::ostream& os) const
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
 * \return The string "Helicoid".
 */
std::string Helicoid::getName() const
{
    return "Helicoid";
}
