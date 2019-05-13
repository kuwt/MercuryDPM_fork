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

#include "Shaft.h"
#include "InteractionHandler.h"
#include "Math/ExtendedMath.h"
#include "Particles/BaseParticle.h"
#include <math.h>

// Last update 23.06.16
// ToDo: check for simplifications in the cos of arc whatever quantities
// ToDo: check for different boundaries depending on maximum possible overlap and particle radius
// ToDo: re-define cosEta_
// ToDo: i believe it can be re-written to have a unique loop just defining if (...) cosBeta = whatever else cosBeta = 1.0;
// ToDo: substitute all the sin(arctan) and similar with the algebraic values

/*!
 * \details Make a Shaft which is centered in the origin, has a length of 1, one
 * revelation, a radius of 1, turns with 1 revelation per second, is infinitely thin
 * and starts at its normal initial point.
 */
Shaft::Shaft()
{
    start_.setZero();
    l_ = 1.0;
    radius_ = 1.0;
    omega_ = 1.0;
    offset_ = 0.0;
    logger(DEBUG, "Shaft() constructor finished.");              
}

/*!
 * \param[in] other The Shaft that has to be copied.
 */
Shaft::Shaft(const Shaft& other)
    : BaseWall(other)
{
    start_ = other.start_;
    l_ = other.l_;
    radius_ = other.radius_;
    omega_ = other.omega_;
    offset_ = other.offset_;
    logger(DEBUG, "Shaft(const Shaft&) copy constructor finished.");
}


/*!
 * \param[in] start A Vec3D which denotes the centre of the lower end of the Shaft.
 * \param[in] l The length of the Shaft, must be positive.
 * \param[in] r The radius of the Shaft, must be positive.
 * \param[in] n The number of revelations of the Shaft, must be positive.
 * \param[in] omega The rotation speed of the Shaft in rev/s.
 * \param[in] thickness The thickness of the Shaft, must be non-negative.
 * \details Make a Shaft by assigning all input parameters to the data-members of
 * this class, and setting the offset_ to 0.
 */
Shaft::Shaft(Vec3D start, Mdouble l, Mdouble r, Mdouble omega)
{
    start_ = start;
    l_ = l;
    radius_ = r;
    omega_ = omega;
    offset_ = 0.0;
    logger(DEBUG, "Shaft(Vec3D, Mdouble, Mdouble, Mdouble, Mdouble, Mdouble) constructor finished.");
}

Shaft::~Shaft()
{
    logger(DEBUG, "~Shaft() finished, destroyed the Shaft.");
}

/*!
 * \return A pointer to a copy of this Shaft.
 */
Shaft* Shaft::copy() const
{
    return new Shaft(*this);
}

// HACK - BEGIN
// need the proper constructor, like in the coil
void Shaft::set(Vec3D start, Mdouble length, Mdouble r, Mdouble omega)
{
    start_ = start;
    l_ = length;
    radius_ = r;
    omega_ = omega;
    offset_ = 0.0;
    
    std::cout << "\n\n Shaft parameters set.\n";
    std::cout << "Start : " << start_ << "\n";
    std::cout << "Length : " << l_ << "\n";
    std::cout << "Radius : " << radius_ << "\n";
    std::cout << "Angular velocity : " << omega_ << "\n";
}

// so I can change the radii on the run
void Shaft::setRadius(Mdouble r)
{
    radius_ = r;
}

// so I can change the omega on the run
void Shaft::setOmega(Mdouble omega)
{
    omega_ = omega;
}

// Function to make the screw rotate by adding \omega*t to the angular offset
void Shaft::rotate(Mdouble dt)
{
    offset_ += omega_ * dt;
}

// the contact detection algorithm
bool Shaft::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    // squared radial distance of the particle from the Shaft axis
    Mdouble rho2 = pow(p.getPosition().X - start_.X, 2) + pow(p.getPosition().Y - start_.Y, 2);
    
    // if the particle is outside of the Shaft returns false
    // check for the radius
    if (rho2 > pow(radius_ + p.getWallInteractionRadius(this), 2)) return false;
    // check for the height
    if (p.getPosition().Z > l_ + start_.Z + p.getWallInteractionRadius(this)) return false;
    if (p.getPosition().Z < start_.Z - p.getWallInteractionRadius(this)) return false;
    
    // radial distance of the particle from the Shaft axis
    Mdouble rho = sqrt(rho2);
    
    // angular coordinate of the particle
    // IMPORTANT: this angle needs to be defined in the interval [0, +2*pi[ radians!
    Mdouble xi = atan2(p.getPosition().Y - start_.Y, p.getPosition().X - start_.X);
    if (xi < 0.0) xi += 2.0*constants::pi;
    
    // trigonometric functions relative to the particle angle
    Mdouble cosXi = (p.getPosition().X - start_.X)/rho;
    Mdouble sinXi = (p.getPosition().Y - start_.Y)/rho;
    
    normal_return.X = -(cosXi*cos(offset_) - sinXi*sin(offset_));
    normal_return.Y = -(sinXi*cos(offset_) + cosXi*sin(offset_));
    normal_return.Z = 0.0;
    
    // distance between the contact point and the particle's centre
    distance = rho - radius_;
    
    return true;
}


/*!
 * \param[in] dt The time for which the Shaft has to be turned.
 */
void Shaft::move_time(Mdouble dt)
{
    offset_ += omega_ * dt;
}

/*!
 * \param[in,out] is Input stream from which the Shaft must be read.
 */
void Shaft::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> start_
            >> dummy >> l_
            >> dummy >> radius_
            >> dummy >> omega_
            >> dummy >> offset_;
}

/*!
 * \param[in,out] is Input stream from which the Shaft must be read.
 * \details Read the Shaft in old style, please note that the thickness is not 
 * read in this function, so it has either to be set manually or it is 0.0 from
 * the default constructor.
 */
void Shaft::oldRead(std::istream& is)
{
    std::string dummy;
    is >> dummy >> start_
            >> dummy >> l_
            >> dummy >> radius_
            >> dummy >> omega_
            >> dummy >> offset_;
}

/*!
 * \param[in,out] os Output stream to which the Shaft must be written.
 */
void Shaft::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " start " << start_
            << " Length " << l_
            << " Radius " << radius_
            << " Omega " << omega_
            << " Offset " << offset_;
}

/*!
 * \return The string "Shaft".
 */
std::string Shaft::getName() const
{
    return "Shaft";
}

/*!
 * \param[in] p Pointer to the BaseParticle which we want to check the interaction for.
 * \param[in] timeStamp The time at which we want to look at the interaction.
 * \param[in] interactionHandler A pointer to the InteractionHandler in which the interaction can be found.
 * \return A pointer to the BaseInteraction that happened between this Shaft
 * and the BaseParticle at the timeStamp.
 */
BaseInteraction*
Shaft::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    if (getDistanceAndNormal(*p,distance,normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        c->setContactPoint(p->getPosition()-(p->getRadius() - 0.5 * c->getOverlap()) * c->getNormal());
        return c;
    }
    return nullptr;
}








