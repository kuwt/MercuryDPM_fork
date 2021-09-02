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


#include "MindlinInteraction.h"
#include "Species/FrictionForceSpecies/MindlinSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <DPMBase.h>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
MindlinInteraction::MindlinInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    slidingSpring_.setZero();
    //k_edit
    //setting the "previousSlidingSpring_" initially to zero, as before a contact
    //this parameter should (obviously) not carry a value.
    slidingSpringPrevious_.setZero();
    //similarly, setting the turning point tangential forces to zero, as a "new"
    //interaction will clearly have no such history
    tangentialForceTurningPointLU_.setZero();
    tangentialForceTurningPointUL_.setZero();
    //k_edit
    //setting the K_t0 parameter by default to zero...
    tangentialStiffnessZero_ = 0;
    //...and its previous value...
    tangentialStiffnessZeroPrevious_ = 0;
    
    absoluteNormalForcePrevious_ = 0;
    //k_edit
    //...and, similarly, the "normal" tangential stiffness, K_t
    tangentialStiffness_ = 0;
    //k_edit
    //setting the tangential force direction by default to -1 (i.e. acting against the motion!)
    tangentialForceDirection_ = -1;
    //k_edit
    //setting to zero the temporary 'holders' which store intermediate values of tangential force
    //and displacement during multiple-step force calculations
    tangentialForceTemp_.setZero();
    tangentialForceTemp2_.setZero();
    tangentialDisplacementTemp_.setZero();
    tangentialDisplacementTemp2_.setZero();
    //Similarly for temporary turning point force values...
    tangentialForceTurningPointLUTemp_.setZero();
    tangentialForceTurningPointULTemp_.setZero();
    //...and the turning point displacement values
    tangentialDisplacementTurningPointUL_.setZero();
    tangentialDisplacementTurningPointLU_.setZero();
    //k_edit
    //setting the variable which stores the minimum displacement for simple loading to zero
    tangentialDisplacementSL_ = 0;
    //k_edit
    //setting the priorLoadingFlag to zero as, when an interaction is created for the first time,
    //there will, by definition, have been no prior loading
    priorLoadingFlag_ = false;
    //ensuring that the initial tangential velocity of a new interaction is zero until updated
    initialTangentialVelocity_.setZero();

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinInteraction::MindlinInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
MindlinInteraction::MindlinInteraction(const MindlinInteraction& p)
        : BaseInteraction(p)
{
    slidingSpring_ = p.slidingSpring_;
    //k_edit
    //setting the K_t0 parameter to that possessed by the interaction
    //to be copied....
    tangentialStiffnessZero_ = p.tangentialStiffnessZero_;
    //...and its previous value...
    tangentialStiffnessZeroPrevious_ = p.tangentialStiffnessZeroPrevious_;
    //...and the same for the "normal" tangential stiffness parameter
    tangentialStiffness_ = p.tangentialStiffness_;
    //k_edit
    //allowing turning point forces to be copied when an interaction is copied
    tangentialForceTurningPointLU_ = p.tangentialForceTurningPointLU_;
    tangentialForceTurningPointUL_ = p.tangentialForceTurningPointUL_;
    //k_edit
    //allowing the direction of the tangential force to be "remembered" also
    //by the copied interaction
    tangentialForceDirection_ = p.tangentialForceDirection_;
    //k_edit
    //copying the temporary 'holders' which store intermediate values of tangential force
    //and displacement during multiple-step force calculations
    tangentialForceTemp_ = p.tangentialForceTemp_;
    tangentialForceTemp2_ = p.tangentialForceTemp2_;
    tangentialDisplacementTemp_ = p.tangentialDisplacementTemp_;
    tangentialDisplacementTemp2_ = p.tangentialDisplacementTemp2_;
    //and similarly for the temporary turning point variables
    tangentialForceTurningPointLUTemp_ = p.tangentialForceTurningPointLUTemp_;
    tangentialForceTurningPointULTemp_ = p.tangentialForceTurningPointULTemp_;
    //...and the displacement equivalents
    tangentialDisplacementTurningPointUL_ = p.tangentialDisplacementTurningPointUL_;
    tangentialDisplacementTurningPointLU_ = p.tangentialDisplacementTurningPointLU_;
    //k_edit
    //copying the variable which stores the minimum displacement for simple loading
    tangentialDisplacementSL_ = p.tangentialDisplacementSL_;
    //k_edit
    priorLoadingFlag_ = p.priorLoadingFlag_;
    //copying the initial tangential velocity
    initialTangentialVelocity_ = p.initialTangentialVelocity_;
    
    absoluteNormalForcePrevious_ = p.absoluteNormalForcePrevious_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinInteraction::MindlinInteraction(const MindlinInteraction& p) finished"<<std::endl;
#endif
}

/*!
 *
 */
MindlinInteraction::MindlinInteraction()
{
    //I don't see why we are restricted in the use of interactions
//#ifdef MERCURY_USE_MPI
//    logger(FATAL,"MindlinInteractions are currently not implemented in parallel MercuryDPM");
//#endif
}

/*!
 *
 */
MindlinInteraction::~MindlinInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"MindlinInteraction::~MindlinInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void MindlinInteraction::write(std::ostream& os) const
{
    //BaseInteraction::write(os);
    os << " slidingSpring " << slidingSpring_;
}

/*!
 * \param[in,out] is
 */
void MindlinInteraction::read(std::istream& is)
{
    //BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> slidingSpring_;
}

//k_new


//k_edit
//Allows the K_t0 parameter (for Mindlin model) to be set
void MindlinInteraction::setTangentialStiffnessZero(Mdouble newKt0)
{
    tangentialStiffnessZero_ = newKt0;
}

//k_edit
//Allows the K_t0 parameter (for Mindlin model) to be accessed
Mdouble MindlinInteraction::getTangentialStiffnessZero()
{
    return tangentialStiffnessZero_;
}

//k_edit
//Allows the "normal" K_t parameter (for Mindlin model) to be accessed
Mdouble MindlinInteraction::getTangentialStiffness()
{
    return tangentialStiffness_;
}

//k_edit
//A quick way to update the K_t0  parameter when necessary
//**************k_TODO at present, only valid for particle-wall interactions and/or interactions*********************
//between similar particles; need to include the reduced parameters for dissimilar particles!
//takes as arguments the relevant particle radius and the relevant shear modulus
//2.0/3.0 as a factor from Di Renzo and Di Maio (2005) eqn (15)
void MindlinInteraction::updateTangentialStiffnessZero(Mdouble rad, Mdouble shearMod)
{
    //K_t0 = 8G_eq * sqrt(R_eq * delta_n) (Di Renzo and Di Maio)
    tangentialStiffnessZero_ = 8 * shearMod * sqrt(rad * getOverlap());
}

//k_new
//A general function to update the tangential stiffness corresponding to an interaction for *any* loading path
//Takes as arguments:
// i) the relevant getSlidingFrictionCoefficientStatic() for the current species,
// ii) The vector direction corresponding to the particles current tangential motion relative to its initial loading
// iii) A boolean to determine whether the turning point force should be considered (false if a 'virgin loading' OR
// if the current tangential force exceeds the magnitude of the relevant turning point force
// iv) a boolean to indicate whether loading or unloading
void MindlinInteraction::updateK_t(const Mdouble fric, const Vec3D direction, const bool useTurningPoint,
                                   const bool isLoading)
{
    //creating a simple positive or negative integer 'multiplier', m, whose values are either +1 or -1
    //that can be used to alter the relevant equations depending on if we are loading or unloading
    Mdouble m;
    //A parameter to temporarily hold the relevant turning point so that this function can work both for
    //loading and loading.
    Vec3D turningPoint;
    //Changing the directions and turning points used in calculation depending on if we are loading or unloading
    if (isLoading)
    {
        m = 1;
        turningPoint = tangentialForceTurningPointUL_;
    }
    else
    {
        m = -1;
        turningPoint = tangentialForceTurningPointLU_;
    }
    //if the current part of the interaction can be treated as a virgin loading
    //(i.e. the turning point can be ignored), we simply update K_t in the 'standard' manner
    if (!useTurningPoint)
    {
        //determining the direction of the current tangential force relative to the
        //current direction of motion
        Mdouble forceDirection = Vec3D::dot(Vec3D::getUnitVector(tangentialForceTemp_), direction);
        //finding the scalar ratio f_t1 / (mu f_n1)
        //(note for a first loading f_t1 is simply equal to f_t, and similarly f_n1 f_n = const
        Mdouble forceRatio = tangentialForceTemp_.getLength() / (fric * getAbsoluteNormalForce());
        //changing the direction of the force ratio, if necessary, based on its relative orientation
        //with respect to the current tangential motion.
        if (forceDirection < 0)
            forceRatio = -forceRatio;
        //using the scalar ratio to update tangential stiffness
        tangentialStiffness_ = tangentialStiffnessZero_ * cbrt(1 - forceRatio);
    }
        //alternatively, if we are reloading and within the limits where the turning point force must be considered
    else
    {
        Vec3D resultantForce = m * (tangentialForceTemp_ - turningPoint);
        //QUESTION: Does it make sense to use the total force here? Or should I just use the tangentialForceTemp?
        Mdouble forceDirection = Vec3D::dot(Vec3D::getUnitVector(resultantForce), direction);
        Mdouble forceRatio = resultantForce.getLength() / (2 * fric * getAbsoluteNormalForce());
        //changing the direction of the force ratio, if necessary, based on its relative orientation
        //with respect to the current tangential motion.
        if (forceDirection < 0)
            forceRatio = -forceRatio;
        //using the scalar ratio to update tangential stiffness
        tangentialStiffness_ = tangentialStiffnessZero_ * cbrt(1 - m * forceRatio);
    }
}

void MindlinInteraction::computeFrictionForce()
{
    //If tangential forces are absent
    if (getAbsoluteNormalForce() == 0.0) return;
    
    const MindlinSpecies* species = getSpecies();//dynamic_cast
    
    //Checking if the relevant particle species has a non-zero friction coefficient
    if (species->getSlidingFrictionCoefficient() != 0.0)
    {
        //Compute the tangential component of relativeVelocity_
        Vec3D tangentialRelativeVelocity = getRelativeVelocity() - getNormal() * getNormalRelativeVelocity();
        {
            // *************************************************************************************************************************
            // **********************************************GENERAL CALCULATIONS*******************************************************
            // *************************************************************************************************************************
            //used to Integrate the spring
            //if particle-wall
            if (dynamic_cast<BaseParticle*>(getI()) == nullptr)
            {
                slidingSpringVelocity_ = tangentialRelativeVelocity;
            }
                //if particle-particle
            else
            {
                slidingSpringVelocity_ = (tangentialRelativeVelocity -
                                          Vec3D::dot(slidingSpring_, getP()->getVelocity() - getI()->getVelocity()) *
                                          getNormal() / getDistance());
            }
            //integrate(getHandler()->timeStep_);
            slidingSpring_ += slidingSpringVelocity_ * getHandler()->getDPMBase()->getTimeStep();
            
            //1) Calculating the current value of K_t0
            //This is identical for both unloading and loading (under constant normal force) and hence can
            //potentially later be moved for neatness/efficiency
            Mdouble r1 = 1.0 * getEffectiveRadius();
            Mdouble shearModulus = species->getEffectiveShearModulus();
            updateTangentialStiffnessZero(r1, shearModulus);
            
            //2) Calculating the relevant tangential dissipation constant based on the current stiffness (and other relevant parameters)
            Mdouble slidingDissipationCoefficient =
                    species->getSlidingDissipation() * sqrt(getEffectiveMass() * tangentialStiffnessZero_);
            
            //3) Using the parameters determined above, calculating the tangential force as per Di Maio and Di Renzo eqn. (37)
            tangentialForce_ = - 2.0/3.0 * tangentialStiffnessZero_ * slidingSpring_;
            
            //applying the Coulomb criterion (Di Maio and Di Renzo eqn. (36)) to ensure falsely large tangential forces are avoided
            //if all is OK, i.e. if tangential force is smaller than mu*f_n, we simply add the tangential force calculated.
            if (tangentialForce_.getLength() < (species->getSlidingFrictionCoefficient() * getAbsoluteNormalForce()))
            {
                tangentialForce_ -= slidingDissipationCoefficient * tangentialRelativeVelocity;
                addForce(tangentialForce_);
            }
                //otherwise, we take the force as mu*f_n and recalculate the spring length accordingly
            else
            {
                tangentialForce_ *= species->getSlidingFrictionCoefficient() * getAbsoluteNormalForce() /
                                    tangentialForce_.getLength();
                //adding the force as mu*f_n
                addForce(tangentialForce_);
                //updating the spring length accordingly
                slidingSpring_ = -(tangentialForce_ / (2.0/3.0 * tangentialStiffnessZero_));
            }
        }
    }
    slidingSpringPrevious_ = slidingSpring_;
}

/*!
 * \param[in] timeStep the dt
 */
void MindlinInteraction::integrate(Mdouble timeStep)
{
    slidingSpring_ += slidingSpringVelocity_ * timeStep;
}

/*!
 * \return Mdouble
 */
Mdouble MindlinInteraction::getElasticEnergy() const
{
    
    //new_edit replacing outdated 'getSlidingStiffness()'
    //return 0.5 * getSpecies()->getSlidingStiffness() * slidingSpring_.getLengthSquared();
    return 0.5 * tangentialStiffness_ * slidingSpring_.getLengthSquared();
}

/*!
 * \return Mdouble
 */
Mdouble MindlinInteraction::getTangentialOverlap() const
{
    ///\todo TWnow this should be positive
    //k_edit Indeed it should - todo done :)
    //k_edit
    //Updating sliding spring to give positive / negative values dependent on the direction of the extension
    return slidingSpring_.getLength() * Vec3D::dot(Vec3D::getUnitVector(slidingSpring_), Vec3D(1.0, 1.0, 1.0));
}

/*!
 * \return const Vec3D
 */
const Vec3D MindlinInteraction::getTangentialForce() const
{
    return tangentialForce_;
}

//k_edit
//Overwrites the function of the same name in the BaseInteraction class
//so that the tangential force direction can be known in the base interaction
//class and hence used to correctly output a **signed** tangential force value.
//(the value was previously always positive!)
const Mdouble MindlinInteraction::getTangentialForceDirection() const
{
    return tangentialForceDirection_;
}

/*!
 * \return const MindlinSpecies*
 */
const MindlinSpecies* MindlinInteraction::getSpecies() const
{
    return static_cast<const MindlinSpecies*>(getBaseSpecies()->getFrictionForce());
;
}

/*!
 * \return std::string
 */
std::string MindlinInteraction::getBaseName() const
{
    return "Mindlin";
}

/*!
 *
 */
void MindlinInteraction::reverseHistory()
{
    slidingSpring_ = -slidingSpring_;
    slidingSpringVelocity_ = -slidingSpringVelocity_;
    tangentialForce_ = -tangentialForce_;
}

void MindlinInteraction::rotateHistory(Matrix3D& rotationMatrix)
{
    slidingSpring_ = rotationMatrix * slidingSpring_;
    slidingSpringVelocity_ = rotationMatrix * slidingSpringVelocity_;
    tangentialForce_ = rotationMatrix * tangentialForce_;
}

//k_edit
/*!
 * \details Return the absolute normal force. This is the magnitude of the normal
 *          force.
 *          \todo Ant: Check this comment.
 * \return  Mdouble that contains the absolute norm (length) of the normal
 *          force.
 */
Mdouble MindlinInteraction::getAbsoluteNormalForcePrevious() const
{
    return absoluteNormalForcePrevious_;
}

//k_edit
/*!
 * \details set absolute normal force.
 * \param[in] absoluteNormalForce   Mdouble contain the value of the absolute
 *                                  normal force.
 */
void MindlinInteraction::setAbsoluteNormalForcePrevious(Mdouble absoluteNormalForcePrevious)
{
    absoluteNormalForcePrevious_ = absoluteNormalForcePrevious;
}
