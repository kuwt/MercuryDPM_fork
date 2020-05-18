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

#include "Walls/SineWall.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"

SineWall::SineWall()
{
    l_ = 1.0;
    sw_wavn_ = 0.0;
    sw_phshift_ = 0.0;
    sw_amp_ = 0.0;
    logger(DEBUG, "SineWall() constructor finished");
}

SineWall::SineWall(const SineWall& other) : BaseWall(other)
{
    l_ = other.l_;
    sw_wavn_ = other.sw_wavn_;
    sw_phshift_ = other.sw_phshift_;
    sw_amp_ = other.sw_amp_;
    logger(DEBUG, "SineWall copy constructor finished");
}

SineWall::SineWall(Mdouble length, Mdouble sw_wavn, Mdouble sw_phshift, Mdouble sw_amp)
{
    l_ = length;
    sw_wavn_ = sw_wavn;
    sw_phshift_ = sw_phshift;
    sw_amp_ = sw_amp;
    logger(DEBUG, "SineWall(params) constructor finished");
}

SineWall::~SineWall()
{
    logger(DEBUG, "SineWall destructor finished");
}

void SineWall::set(Mdouble length, Mdouble sw_wavn, Mdouble sw_phshift, Mdouble sw_amp)
{
    l_ = length;
    sw_wavn_ = sw_wavn;
    sw_phshift_ = sw_phshift;
    sw_amp_ = sw_amp;
}

SineWall* SineWall::copy() const
{
    return new SineWall(*this);
}

bool SineWall::getDistanceAndNormal(const BaseParticle& p, Mdouble& distance, Vec3D& normal_return) const
{
    /* Define some shortcuts */
    double x0 = p.getPosition().X;
    double y0 = p.getPosition().Y;
    double z0 = p.getPosition().Z;
    double a = p.getRadius();
    
    // more shortcuts
    double A = sw_amp_;
    double k = sw_wavn_;
    double ph = sw_phshift_;
    
    /* First, check whether the particle is definitely out of contact with the
     * chute. This will be so if the particle's interaction radius is small and
     * the particle's z-position is high. */
    // TODO do this properly
    if (z0 > A + a)
        return false;
        
        /* Or, if the particle lies under the surface, then it is expelled
         * perpendicularly. This shouldn't happen very often provided that the step
         * size is small.*/
    else if (z0 < A * sin(k * x0 + ph))
    {
        fprintf(stderr, "Particle at x0=%f, z0=%f is under SineWall\n", x0, z0);
        // exit(-1);
        
        distance = A * sin(k * x0 + ph) - z0;
        normal_return = Vec3D(0, 0, 1);
        return true;
    }
        
        /* Otherwise, find the contact point by using Newton's method to minimise
         * (half of the squared) distance */
    else
    {
        // Iterate according to Newton-Raphson
        double q = x0;
        double correction = 0;
        // fprintf(stderr, "x0 = %lf, z0 = %lf, q = %lf\n", x0, z0, q);
        do
        {
            double z = A * sin(k * q + ph);
            double dzdq = A * k * cos(k * q + ph);
            double d2zdq2 = -A * pow(k, 2) * sin(k * q + ph);
            
            // (Half of the square of the) distance. Not used.
            double R = 0.5 * (pow(q - x0, 2) + pow(z - z0, 2));
            // derivative of half the squared distance -- to be zeroed
            double dRdq = q - x0 + (z - z0) * dzdq;
            // second derivative of half the squared distance
            double d2Rdq2 = 1 + (z - z0) * d2zdq2 + pow(dzdq, 2);
            // The required correction to our current estimate. (Note the sign
            // convention that this is to be subtracted.)
            correction = dRdq / d2Rdq2;
            //fprintf(stderr, "x0 = %.12lf, z0 = %.12lf, q = %.12lf, dRdq = %e, d2Rdq2 = %e, correction = %.e\n", 
            //       x0, z0, q, dRdq, d2Rdq2, correction);
            // N-R update
            q -= correction;
            // Perhaps this form is useful for reducing rounding errors?
            // q = (q*d2Rdq2 - dRdq)/d2Rdq2;         
        } while (fabs(correction) > 1e-14);
        //fprintf(stderr, "correction = %e\n", fabs(correction));
        
        // The standard Newton-Raphson iteration won't work, because the problem is
        // too oscillatory (and N-R will converge on a near-root minimum, not the
        // root). Neither does the following, a modified N-R which takes second and
        // third derivatives into account. 
        /*
           double q = x0;
           double correction = 0;
           do 
           {
        // derivative of half the squared distance
        double dRdq = q - x0 - A*pow(k,2)*z0*cos(k*q+ph) + 0.5*pow(A,2)*pow(k,3)*sin(2*k*q+2*ph);
        // second derivative
        double d2Rdq2 = 1 + pow(A,2)*pow(k,4)*cos(2*k*q+2*ph) + A*pow(k,3)*z0*sin(k*q+ph);
        // third derivative
        double d2Rdq2d = A*pow(k,4)*cos(k*q+ph) * (z0 - 4*A*k*sin(k*q+ph));

        double discriminant = d2Rdq2*d2Rdq2 - 2*dRdq*d2Rdq2d;
        double newq = discriminant > 0 ? 
        q - (d2Rdq2 - sqrt(discriminant)) / d2Rdq2d
        : q - dRdq / d2Rdq2;
        double correction = newq/q - 1;
        fprintf(stderr, "x0 = %.12lf, z0 = %.12lf, q = %.12lf, newq = %.12lf, correction = %e\n",
        x0, z0, q, newq, correction);
        q = newq;
        } while (fabs(correction) > 1e-7);
        fprintf(stderr, "converged, fabs(correction) = %e\n", fabs(correction));
        */
        
        Mdouble distanceSquared
                = pow(q - x0, 2) + pow(sw_amp_ * sin(sw_wavn_ * q + sw_phshift_) - z0, 2);
        
        if (distanceSquared >= mathsFunc::square(p.getWallInteractionRadius(this))
            && z0 > sw_amp_ * sin(sw_wavn_ * q + sw_phshift_))
        {
            return false;
        }
        else
        {
            Vec3D ContactPoint;
            distance = sqrt(distanceSquared);
            ContactPoint.X = q;
            ContactPoint.Y = y0;
            ContactPoint.Z = sw_amp_ * sin(sw_wavn_ * q + sw_phshift_);
            // normal_return = ContactPoint - p.getPosition();
            normal_return = ContactPoint - p.getPosition();
            normal_return /= normal_return.getLength();
            return true;
        }
    }
}

BaseInteraction*
SineWall::getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler)
{
    Mdouble distance;
    Vec3D normal;
    if (getDistanceAndNormal(*p, distance, normal))
    {
        BaseInteraction* c = interactionHandler->getInteraction(p, this, timeStamp);
        c->setNormal(-normal);
        c->setDistance(distance);
        c->setOverlap(p->getRadius() - distance);
        // c->setContactPoint( p->getPosition() - (p->getRadius() - 0.5*c->getOverlap()) * c->getNormal());
        c->setContactPoint(p->getPosition() + distance * normal);
        /// \todo Hacked please fix
        return c;
    }
    else
        return nullptr;
}

/*!
 * \param[in,out] is The input stream from which the SineWall is read.
 */
void SineWall::read(std::istream& is)
{
    BaseWall::read(is);
    std::string dummy;
    is >> dummy >> l_
       >> dummy >> sw_wavn_
       >> dummy >> sw_phshift_
       >> dummy >> sw_amp_;
}

/*!
 * \param[in,out] os The outpus stream to which the SineWall is written.
 */
void SineWall::write(std::ostream& os) const
{
    BaseWall::write(os);
    os << " Length " << l_
       << " sw_wavn " << sw_wavn_
       << " sw_phshift " << sw_phshift_
       << " sw_amp " << sw_amp_;
}

/*!
 * \return The string "SineWall".
 */
std::string SineWall::getName() const
{
    return "SineWall";
}
