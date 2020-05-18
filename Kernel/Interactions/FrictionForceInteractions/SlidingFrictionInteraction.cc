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


#include "SlidingFrictionInteraction.h"
#include "Species/FrictionForceSpecies/SlidingFrictionSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>
#include <DPMBase.h>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
SlidingFrictionInteraction::SlidingFrictionInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    slidingSpring_.setZero();
    isSuperQuadricInteraction_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SlidingFrictionInteraction::SlidingFrictionInteraction() finished"<<std::endl;
#endif
}

//used for mpi
SlidingFrictionInteraction::SlidingFrictionInteraction()
        : BaseInteraction()
{
    slidingSpring_.setZero();
    isSuperQuadricInteraction_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SlidingFrictionInteraction::SlidingFrictionInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
SlidingFrictionInteraction::SlidingFrictionInteraction(const SlidingFrictionInteraction& p)
        : BaseInteraction(p)
{
    slidingSpring_ = p.slidingSpring_;
    isSuperQuadricInteraction_ = p.isSuperQuadricInteraction_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SlidingFrictionInteraction::SlidingFrictionInteraction(const SlidingFrictionInteraction& p) finished"<<std::endl;
#endif
}

/*!
 *
 */
SlidingFrictionInteraction::~SlidingFrictionInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SlidingFrictionInteraction::~SlidingFrictionInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void SlidingFrictionInteraction::write(std::ostream& os) const
{
    //BaseInteraction::write(os);
    os << " slidingSpring " << slidingSpring_;
}

/*!
 * \param[in,out] is
 */
void SlidingFrictionInteraction::read(std::istream& is)
{
    //BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> slidingSpring_;
}

/*!
 *
 */
void SlidingFrictionInteraction::computeFrictionForce()
{
    //If tangential forces are absent
    if (getAbsoluteNormalForce() == 0.0) return;
    
    const SlidingFrictionSpecies* species = getSpecies();//dynamic_cast
    isSuperQuadricInteraction_ = species->getIsSuperquadricSpecies();
    
    if (species->getSlidingFrictionCoefficient() != 0.0)
    {
        //Compute the tangential component of relativeVelocity_
        // relativeVelocity = [v_p + (r_p-c) x w_p] - [v_i + (r_i-c) x w_i]
        // tangentialRelativeVelocity = relativeVelocity - (relativeVelocity . n) n
        Vec3D tangentialRelativeVelocity = getRelativeVelocity() - getNormal() * getNormalRelativeVelocity();
        
        if (species->getSlidingStiffness() != 0.0)
        {
            if (!isSuperQuadricInteraction_)
            {
                computeSlidingSpring(tangentialRelativeVelocity);
            }
            else
            {
                computeSlidingSpringSuperQuadric(tangentialRelativeVelocity);
            }
            
            //Calculate test force acting on P including viscous force
            tangentialForce_ = -species->getSlidingStiffness() * slidingSpring_ -
                               species->getSlidingDissipation() * tangentialRelativeVelocity;
            if (getBaseSpecies()->getNormalForce()->getConstantRestitution()) tangentialForce_ *=2.0*getEffectiveMass();
            
            //tangential forces are modelled by a spring-damper of elasticity kt and viscosity dispt (sticking),
            //but the force is limited by Coulomb friction (sliding):
            Mdouble tangentialForceSquared = tangentialForce_.getLengthSquared();
            if (tangentialForceSquared <=
                mathsFunc::square(species->getSlidingFrictionCoefficientStatic() * getAbsoluteNormalForce()))
            {
                //if sticking (|ft|<=mu*|fn|), apply the force
                addForce(tangentialForce_);
            }
            else
            {
                //if sliding, resize the tangential force such that |ft|=mu*|fn|
                tangentialForce_ *= species->getSlidingFrictionCoefficient() * getAbsoluteNormalForce() /
                                    std::sqrt(tangentialForceSquared);
                addForce(tangentialForce_);
                //resize the tangential spring accordingly such ft=-k*delta-nu*relVel
                slidingSpring_ = -(tangentialForce_ + species->getSlidingDissipation() * tangentialRelativeVelocity) /
                                 species->getSlidingStiffness();
            }
        }
        else //if no spring stiffness is set
        {
//            if (species->getSlidingDissipation()==0.0)
//            {
//                std::cerr << "SlidingFrictionInteraction::getForce(): warning:  both sliding stiffness and dissipation are zero" << std::endl;
//            }
            Mdouble tangentialRelativeVelocitySquared = tangentialRelativeVelocity.getLengthSquared();
            if (tangentialRelativeVelocitySquared * mathsFunc::square(species->getSlidingDissipation()) <=
                mathsFunc::square(species->getSlidingFrictionCoefficientStatic() * getAbsoluteNormalForce()))
                tangentialForce_ = -species->getSlidingDissipation() * tangentialRelativeVelocity;
            else //if sliding, set force to Coulomb limit
                tangentialForce_ = -(species->getSlidingFrictionCoefficient() * getAbsoluteNormalForce() /
                                     std::sqrt(tangentialRelativeVelocitySquared)) * tangentialRelativeVelocity;
            
            addForce(tangentialForce_);
        }
    }

    /* And if species->getSlidingFrictionCoefficient() == 0.0, then there is no
     * friction force at all. */
}

void SlidingFrictionInteraction::computeSlidingSpring(const Vec3D& tangentialRelativeVelocity)
{
    //used to Integrate the spring
    if (dynamic_cast<BaseParticle*>(getI()) == nullptr)  //if particle-wall
    {
        slidingSpringVelocity_ = tangentialRelativeVelocity;
    }
    else //if particle-particle
    {
        slidingSpringVelocity_ = (tangentialRelativeVelocity -
                                  Vec3D::dot(slidingSpring_, getP()->getVelocity() - getI()->getVelocity()) *
                                  getNormal() / getDistance());
    }
    // v_s = v_t - [xi . (v_p-v_i)/|r_pi|] n
    
    
    //integrate(getHandler()->timeStep_);
    // xi = xi' + dt v_s
    slidingSpring_ += slidingSpringVelocity_ * getHandler()->getDPMBase()->getTimeStep();
    // Stefan does [EJECE-12/2008] sth. like xi = xi' - (xi . n) n + dt*v_t, see SlidingFrictionSuperQuadricInteraction
}

void SlidingFrictionInteraction::computeSlidingSpringSuperQuadric(const Vec3D& tangentialRelativeVelocity)
{
    //used to Integrate the spring
    if (dynamic_cast<BaseParticle*>(getI()) == nullptr)  //if particle-wall
    {
        slidingSpringVelocity_ = tangentialRelativeVelocity;
    }
    else //if particle-particle
    {
        slidingSpringVelocity_ = tangentialRelativeVelocity;
    }
    
    // From Stefan [EJECE-12/2008]:
    // xi = xi' - (xi . n) n and renormalising, then add v_t * dt
    logger(DEBUG, "sliding spring, normal before rotation-step: %, %", slidingSpring_, getNormal());
    if (slidingSpring_.getLength() > 1e-10)
    {
        const Mdouble springLength = slidingSpring_.getLength();
        slidingSpring_ -= Vec3D::dot(slidingSpring_, getNormal()) * getNormal();
        slidingSpring_.normalise();
        slidingSpring_ *= springLength;
        logger.assert(std::abs(slidingSpring_.getLength() - springLength) < 1e-10, "Spring length not the same after rotation");
    }
    logger(DEBUG, "sliding spring after rotation-step: %\n", slidingSpring_);
    slidingSpring_ += slidingSpringVelocity_ * getHandler()->getDPMBase()->getTimeStep();
    logger.assert(std::abs(Vec3D::dot(slidingSpring_, getNormal())) < 1e-10, "sliding spring not perpendicular to normal");
}


/*!
 * \param[in] timeStep the dt
 */
void SlidingFrictionInteraction::integrate(Mdouble timeStep)
{
    slidingSpring_ += slidingSpringVelocity_ * timeStep;
}

/*!
 * \return Mdouble
 */
Mdouble SlidingFrictionInteraction::getElasticEnergy() const
{
    if (getBaseSpecies()->getNormalForce()->getConstantRestitution()) {
        return getEffectiveMass() * getSpecies()->getSlidingStiffness() * slidingSpring_.getLengthSquared();
    } else {
        return 0.5 * getSpecies()->getSlidingStiffness() * slidingSpring_.getLengthSquared();
    }
}

/*!
 * \return Mdouble
 */
Mdouble SlidingFrictionInteraction::getTangentialOverlap() const
{
    ///\todo TWnow this should be positive
    return -slidingSpring_.getLength();
}

/*!
 * \return const Vec3D
 */
const Vec3D SlidingFrictionInteraction::getTangentialForce() const
{
    return tangentialForce_;
}

/*!
 * \return const SlidingFrictionSpecies*
 */
const SlidingFrictionSpecies* SlidingFrictionInteraction::getSpecies() const
{
    return static_cast<const SlidingFrictionSpecies*>(getBaseSpecies()->getFrictionForce());
;
}

/*!
 * \return std::string
 */
std::string SlidingFrictionInteraction::getBaseName() const
{
    return "SlidingFriction";
}

///\bug may need to remove the consts
void SlidingFrictionInteraction::setSlidingSpring(const Vec3D slidingSpring)
{
    slidingSpring_ = slidingSpring;
}

/*!
 *
 */
void SlidingFrictionInteraction::reverseHistory()
{
    slidingSpring_ = -slidingSpring_;
    slidingSpringVelocity_ = -slidingSpringVelocity_;
    tangentialForce_ = -tangentialForce_;
}

void SlidingFrictionInteraction::rotateHistory(Matrix3D& rotationMatrix)
{
    slidingSpring_ = rotationMatrix * slidingSpring_;
    slidingSpringVelocity_ = rotationMatrix * slidingSpringVelocity_;
    tangentialForce_ = rotationMatrix * tangentialForce_;
}

Vec3D SlidingFrictionInteraction::getSlidingSpring() const
{
    return slidingSpring_;
}

void SlidingFrictionInteraction::moveSlidingSpring(const Vec3D displacement)
{
    slidingSpring_ += displacement;
}

void SlidingFrictionInteraction::setIsSuperQuadricInteraction(bool isSuperQuadricInteraction)
{
    isSuperQuadricInteraction_ = isSuperQuadricInteraction;
}
