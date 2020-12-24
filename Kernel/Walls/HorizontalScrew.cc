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

#include "HorizontalScrew.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
using mathsFunc::square;

/*!
 * \details Make a HorizontalScrew which is centered in the origin, has a length of 1, one
 * revelation, a radius of 1, turns with 1 revelation per second, is infinitely thin
 * and starts at its normal initial point.
 */
HorizontalScrew::HorizontalScrew()
{
    start_.setZero();
    l_ = 1.0;
    minR_ = 0.0;
    lowerR_ = 1.0;
    diffR_ = 0.0;
    n_ = 1.0;
    omega_ = 1.0;
    offset_ = 0.0;
    thickness_ = 0.0;
    bladeLength_=0;
    bladeWidth_=0;
    bladeMounts_={};
    logger(DEBUG, "HorizontalScrew() constructor finished.");
}

/*!
 * \param[in] other The HorizontalScrew that has to be copied.
 */
HorizontalScrew::HorizontalScrew(const HorizontalScrew& other)
    : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    minR_ = other.minR_;
    lowerR_ = other.lowerR_;
    diffR_ = other.diffR_;
    n_ = other.n_;
    omega_ = other.omega_;
    thickness_ = other.thickness_;
    offset_ = other.offset_;
    bladeLength_=other.bladeLength_;
    bladeWidth_=other.bladeWidth_;
    bladeMounts_=other.bladeMounts_;
    logger(DEBUG, "HorizontalScrew(const HorizontalScrew&) copy constructor finished.");
}

/*!
 * \param[in] start A Vec3D which denotes the centre of the lower end of the HorizontalScrew.
 * \param[in] l The length of the HorizontalScrew, must be positive.
 * \param[in] r The radius of the HorizontalScrew, must be positive.
 * \param[in] n The number of revelations of the HorizontalScrew, must be positive.
 * \param[in] omega The rotation speed of the HorizontalScrew in rev/s.
 * \param[in] thickness The thickness of the HorizontalScrew, must be non-negative.
 * \details Make a HorizontalScrew by assigning all input parameters to the data-members of
 * this class, and setting the offset_ to 0.
 */
HorizontalScrew::HorizontalScrew(Vec3D start, Mdouble l, Mdouble minR, Mdouble lowerR, Mdouble diffR, Mdouble n, Mdouble omega, Mdouble thickness, const ParticleSpecies* s)
{
    start_ = start;
    l_ = l;
    minR_ = minR;
    lowerR_ = lowerR;
    diffR_ = diffR;
    n_ = n;
    omega_ = omega;
    thickness_ = thickness;
    offset_ = 0.0;
    bladeLength_=0;
    bladeWidth_=0;
    bladeMounts_={};
    setSpecies(s);
    logger(DEBUG, "HorizontalScrew(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

HorizontalScrew::~HorizontalScrew()
{
    logger(DEBUG, "~HorizontalScrew() finished, destroyed the HorizontalScrew.");
}

/*!
 * \return A pointer to a copy of this HorizontalScrew.
 */
HorizontalScrew* HorizontalScrew::copy() const
{
    return new HorizontalScrew(*this);
}

/*!
 * \param[in] p BaseParticle we want to calculate the distance and whether it collided of.
 * \param[out] distance The distance of the BaseParticle to this wall.
 * \param[out] normal_return If there was a collision, the normal vector to this wall will be placed here.
 * \return A boolean which says whether or not there was a collision.
 * \details This function computes whether or not there is a collision between
 * a given BaseParticle and this HorizontalScrew. If there is a collision, this
 * function also computes the distance between the BaseParticle and HorizontalScrew
 * and the normal of the IntersectionOfWalls at the intersection point.
 * \todo Make this function readable and explain the steps in the details.
 */
bool HorizontalScrew::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    Mdouble RSquared = square(p.getPosition().X - start_.X) + square(p.getPosition().Y - start_.Y);
    Mdouble Z = p.getPosition().Z - start_.Z;
    //first do a simple check if particle is within the cylindrical hull of the screw
    Mdouble maxR = std::max(minR_,lowerR_+diffR_*Z/l_)+bladeWidth_;//todo: fix
    Mdouble interactionRadius = p.getWallInteractionRadius(this) + thickness_;
    if (RSquared > square(maxR + interactionRadius)
        || Z > l_ + interactionRadius
        || Z < - interactionRadius)
    {
        return false;
    }

    //else:
    Mdouble R = sqrt(RSquared);
    Mdouble A = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);

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
    Mdouble q0 = Z / l_; //assume closest point q0 is at same z-location as the particle

    //Set initial q to the closest position on the screw with the same angle as A
    //The initial guess will be in the minimum of the sin closest to q0
    //Minima of the sin are at
    //A-2*pi*(offset+N*q)=k*pi (k=integer)
    //q=A/(2*pi*N)-k/(2*N)-offset/N (k=integer)

    Mdouble k = round(A / constants::pi - 2.0 * (offset_ + n_ * q0)); // k: |A-a(q0)-k*pi|=min
    q = A / (2.0 * constants::pi * n_) - k / (2.0 * n_) - offset_ / n_; // q: a(q)=A

    //this makes the trioliet screw unique: only one turn
    if (((int)k)%2==0)
        return false;

    //Now apply Newton's method to find the rel height q of the contact point
    do
    {
        Mdouble arg = 2.0 * A - 4.0 * constants::pi * (n_ * q + offset_);
        dd = -2.0 * constants::pi * n_ * RSquared * sin(arg) - 2.0 * l_ * (Z - q * l_);
        ddd = 8.0 * constants::sqr_pi * n_ * n_ * RSquared * cos(arg) + 2.0 * l_ * l_;
        q -= dd / ddd;
    } while (fabs(dd / ddd) > 1e-14);

    //Calculate r
    Mdouble r = R * cos(2.0 * constants::pi * (offset_ + n_ * q) - A);
    maxR = std::max(minR_,lowerR_+diffR_*q);//todo: fix
    for (auto b : bladeMounts_)
    {
        Mdouble dq = q-b;
        if (dq>0 && dq<bladeLength_)
        {
            maxR += bladeWidth_*(dq/bladeLength_);
        }
    }

    //Check if the location is actually on the screw:
    //First possibility is that the radius is too large:
    if (fabs(r) > maxR) //Left boundary of the coil
    {
        //either the contact point is too far from the screw ...
        if (fabs(r) > maxR+thickness_)
            return false;
        //... or we have to compute the contact around the edge
        r = mathsFunc::sign(r) * maxR;
        unsigned int steps = 0;
        //This case reduces to the coil problem
        do
        {
            dd = -4.0 * R * r * constants::pi * n_ * sin(A - 2.0 * constants::pi * (n_ * q + offset_)) - 2.0 * l_ * (Z - q * l_);
            ddd = 8.0 * R * r * constants::sqr_pi * n_ * n_ * cos(A - 2.0 * constants::pi * (n_ * q + offset_)) + 2.0 * l_ * l_;
            q -= dd / ddd;
            steps++;
        } while (fabs(dd / ddd) > 1e-14);
    }
    //Second possibility is that it occurred before the start of after the end
    if (q < 0)
    {
        q = 0;
        r = R * cos(A - 2.0 * constants::pi * (offset_ + q * n_));
        if (fabs(r) > maxR)
        {
            r = mathsFunc::sign(r) * maxR;
        }
    }
    else if (q > 1)
    {
        q = 1;
        r = R * cos(A - 2.0 * constants::pi * (offset_ + q * n_));
        if (fabs(r) > maxR)
        {
            r = mathsFunc::sign(r) * maxR;
        }
    }

    Mdouble distanceSquared = square(R*sin(A - 2 * constants::pi * (offset_ + n_ * q))) + square(Z - q * l_);
    //If distance is too large there is no contact
    if (distanceSquared >= square(interactionRadius))
    {
        return false;
    }

    Vec3D ContactPoint;
    distance = sqrt(distanceSquared) - thickness_;
    ContactPoint.X = start_.X + r * cos(2.0 * constants::pi * (offset_ + n_ * q));
    ContactPoint.Y = start_.Y + r * sin(2.0 * constants::pi * (offset_ + n_ * q));
    ContactPoint.Z = start_.Z + q * l_;
    normal_return = (ContactPoint - p.getPosition());
    normal_return.normalise();
    return true;
}

Mdouble HorizontalScrew::getLength() const
{
    return l_;
}

Vec3D HorizontalScrew::getStart() const
{
    return start_;
}

/*!
 * \param[in] dt The time for which the HorizontalScrew has to be turned.
 */
void HorizontalScrew::move_time(Mdouble dt)
{
    offset_ += omega_ * dt;
}

/*!
 * \param[in,out] is Input stream from which the HorizontalScrew must be read.
 */
void HorizontalScrew::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    unsigned n=0;
    is >> dummy >> start_
            >> dummy >> l_
            >> dummy >> minR_
            >> dummy >> lowerR_
            >> dummy >> diffR_
            >> dummy >> n_
            >> dummy >> omega_
            >> dummy >> thickness_
            >> dummy >> offset_
            >> dummy >> bladeLength_
            >> dummy >> bladeWidth_
            >> dummy >> n;
    for (unsigned i=0; i<n; i++) {
        Mdouble val;
        is >> val;
        bladeMounts_.push_back(val);
    }
}

/*!
 * \param[in,out] os Output stream to which the HorizontalScrew must be written.
 */
void HorizontalScrew::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << start_
            << " length " << l_
            << " minRadius " << minR_
            << " lowerRadius " << lowerR_
            << " diffRadius " << diffR_
            << " revolutions " << n_
            << " omega " << omega_
            << " thickness " << thickness_
            << " offset " << offset_
            << " bladeLength " << bladeLength_
            << " bladeWidth " << bladeWidth_
            << " bladeMounts " << bladeMounts_.size();
    for (auto n : bladeMounts_)
        os << " " << n;
}

/*!
 * \return The string "HorizontalScrew".
 */
std::string HorizontalScrew::getName() const
{
    return "HorizontalScrew";
}

void HorizontalScrew::setBlades(const Mdouble bladeWidth, const Mdouble bladeLength,const std::vector<Mdouble> bladeMounts)
{
    bladeWidth_ = bladeWidth;
    bladeLength_ = bladeLength;
    bladeMounts_ = bladeMounts;
}

void HorizontalScrew::writeVTK (VTKContainer& vtk) const
{
    unsigned nr = 10;
    unsigned nz = 180;

    unsigned nPoints = vtk.points.size();
    Vec3D contactPoint;
    for (unsigned iz=0; iz<nz; iz++) {
        double maxR = std::max(minR_,lowerR_+(double)iz/nz*diffR_);//todo: fix
        for (auto b : bladeMounts_)
        {
            Mdouble dq = (double)iz/nz-b;
            if (dq>0 && dq<bladeLength_)
            {
                maxR += bladeWidth_*(dq/bladeLength_);
            }
        }
        for (unsigned ir=0; ir<nr; ir++) {
            double q = (double)iz/nz;
            double r = (double)ir/nr*maxR;
            contactPoint.X = start_.X - r * cos(2.0 * constants::pi * (offset_ + n_ * q));
            contactPoint.Y = start_.Y - r * sin(2.0 * constants::pi * (offset_ + n_ * q));
            contactPoint.Z = start_.Z + q * l_;
            vtk.points.push_back(contactPoint);
        }
    }

    unsigned nCells = vtk.triangleStrips.size();
    //vtk.triangleStrips.reserve(nCells+(nz-1));
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
