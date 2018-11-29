//Copyright (c) 2013-2017, The MercuryDPM Developers Team. All rights reserved.
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

#include "Coil2.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"

Coil2::Coil2()
        : BaseWall()
{
    tangent0_.setZero();
    tangent1_.setZero();
    length_ = 1;
    radius_ = 1;
    windings_ = 1;
    thickness_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "Coil2() finished" << std::endl;
#endif              
}

Coil2::Coil2(Vec3D orientation, Vec3D tangent, Mdouble L, Mdouble r, Mdouble N, Mdouble thickness)
        : BaseWall()
{
    set(orientation, tangent, L, r, N, thickness);
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "Coil2() finished" << std::endl;
#endif              
}

void Coil2::set(Vec3D orientation, Vec3D tangent, Mdouble L, Mdouble r, Mdouble N, Mdouble thickness)
{
    orientation.normalise();
    setOrientation(orientation);
    tangent0_ = tangent - Vec3D::dot(tangent, orientation) * orientation;
    tangent0_.normalise();
    tangent1_ = Vec3D::cross(orientation, tangent0_);
    length_ = L;
    radius_ = r;
    windings_ = N;
    thickness_ = thickness;
}

Coil2* Coil2::copy() const
{
    return new Coil2(*this);
}

bool Coil2::getDistanceAndNormal(const BaseParticle &P, Mdouble &distance, Vec3D &normal_return) const
{
    //we use the coordinate system (normal,tangent0,tangent1), with axes (orientation_, tangent0_, tangent1_) and origin position_
    //first, calculate (normal,tangent0,tangent1), which are the coordinates of the position of P
    //in the coordinate system with axes directions (orientation_, tangent0_, tangent1_) and origin (position_)
    Vec3D PO = P.getPosition() -getPosition();
    ///\todo IFCD: now put .getAxis after getOrientation to make it compile again. Please check if it is correct!
    Mdouble normal = Vec3D::dot(PO, getOrientation().getAxis());
    Vec3D tangentialVector = PO - normal * getOrientation().getAxis();
    Mdouble tangentSquared = tangentialVector.getLengthSquared();
    //check if the position of P is far enough from the coil axis to prohibit any contact
    if (   tangentSquared > mathsFunc::square(radius_ + P.getWallInteractionRadius() + thickness_)
        || tangentSquared < mathsFunc::square(radius_ - P.getWallInteractionRadius() - thickness_)
        || normal > length_ + P.getWallInteractionRadius() + thickness_
        || normal < - P.getWallInteractionRadius() - thickness_)
        return false;
    Mdouble tangent0 = Vec3D::dot(PO, tangent0_);
    Mdouble tangent1 = Vec3D::dot(PO, tangent1_);

    //get the angle and and distance from the coil axis
    Mdouble angle = atan2(tangent1, tangent0);
    Mdouble r = sqrt(tangentSquared);

    ///To find the contact point we have to minimize (with respect to q)
    ///distance^2=(tangent0-radius*cos(2*pi*windings*q))^2+(tangent1-radius*sin(2*pi*windings*q))^2+(normal-q*length)^2
    ///Using polar coordinates (i.e. tangent0=r*cos(angle), y-y0=r*sin(angle)
    ///distance^2=r^2+radius^2-2*r*radius*cos(angle-2*Pi*(windings*q))+(normal-q*L)^2
    ///This, we have to minimize. We use the Euler algoritm
    
    Mdouble q; //Current guess
    Mdouble dd; //Derivative at current guess
    Mdouble ddd; //Second derivative at current guess
    Mdouble q0 = normal / length_; //Minimum of the parabolic part
            
    ///The initial guess will be in the maximum of the cos closest to the minimum of the parabolic part
    ///Minima of the cos are at
    ///angle-2*Pi*windings*q=2*k*Pi (k=integer)
    ///q=angle/(2*Pi*windings)-k/windings (k=integer)
    
    Mdouble k = round(angle / 2.0 / constants::pi - windings_ * q0);
    q = angle / (2 * constants::pi * windings_) - k / windings_;
    
    //Now apply Newton's method
    do
    {
        dd = -4.0 * r * radius_ * constants::pi * windings_ * sin(angle - 2.0 * constants::pi * (windings_ * q)) - 2.0 * length_ * (normal - q * length_);
        ddd = 8.0 * r * radius_ * constants::sqr_pi * windings_ * windings_ * cos(angle - 2.0 * constants::pi * (windings_ * q)) + 2.0 * length_ * length_;
        q -= dd / ddd;
    } while (fabs(dd / ddd) > 1e-14);
    
    //Check if the location is actually on the coil, otherwise a point collision with the end of the coil calculated
    if (q < 0) //Left boundary of the coil
        q = 0;
    else if (q > 1) //right boundary of the coil
        q = 1;

    Mdouble distanceSquared = r * r + radius_ * radius_ - 2 * r * radius_ * std::cos(angle - 2 * constants::pi * windings_ * q) + mathsFunc::square(normal - q * length_);
    //If distance is too large there is no contact
    if (distanceSquared >= (P.getWallInteractionRadius() + thickness_) * (P.getWallInteractionRadius() + thickness_))
    {
        //std::cout<<"Particle is out of second bound checking, distance^2="<<distance<<" max="<<(P.getRadius()+thickness_)*(P.getRadius()+thickness_)<<std::endl;
        return false;
    }

    Vec3D contactPoint;
    distance = sqrt(distanceSquared) - thickness_;
    //here, we ignore the thickness
    ///\todo IFCD: now put .getAxis after getOrientation to make it compile again. Please check if it is correct!
    contactPoint = getPosition() + q*length_*getOrientation().getAxis()
        + radius_ * std::cos(2.0 * constants::pi * windings_ * q) * tangent0_
        + radius_ * std::sin(2.0 * constants::pi * windings_ * q) * tangent1_;
    normal_return = contactPoint - P.getPosition();
    normal_return.normalise();
    return true;
}

///Allows the wall to be moved to a new position (also orthogonal to the normal), and setting the velocity
void Coil2::move_time(Mdouble dt)
{
//    offset_ += omega_ * dt;
}

///reads wall
void Coil2::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> tangent0_
        >> dummy >> length_
        >> dummy >> radius_
        >> dummy >> windings_
        >> dummy >> thickness_;
    
    ///\todo IFCD: now put .getAxis after getOrientation to make it compile again. Please check if it is correct!
    set(getOrientation().getAxis(), tangent0_, length_, radius_, windings_, thickness_);
}

void Coil2::oldRead(std::istream& is)
{
//    std::string dummy;
//    is >> dummy >> tangent0_
//        >> dummy >> length_
//        >> dummy >> radius_ >> dummy >> windings_
//        >> dummy >> omega_
//        >> dummy >> offset_;
}

///outputs wall
void Coil2::write(std::ostream& os) const
{
    BaseWall::write(os);
    os  << " tangent " << tangent0_
        << " length " << length_
        << " radius " << radius_
        << " windings " << windings_
        << " thickness " << thickness_;
}

std::string Coil2::getName() const
{
    return "Coil2";
}
