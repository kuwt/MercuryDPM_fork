//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


#include "MindlinRollingTorsionInteraction.h"
#include "Species/FrictionForceSpecies/MindlinRollingTorsionSpecies.h"
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
MindlinRollingTorsionInteraction::MindlinRollingTorsionInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                   unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp), MindlinInteraction(P, I, timeStamp)
{
    rollingSpring_.setZero();
    torsionSpring_.setZero();
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinRollingTorsionInteraction::MindlinRollingTorsionInteraction() finished"<<std::endl;
#endif
}

//used for mpi
MindlinRollingTorsionInteraction::MindlinRollingTorsionInteraction()
        : BaseInteraction(), MindlinInteraction()
{
    rollingSpring_.setZero();
    torsionSpring_.setZero();
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinRollingTorsionInteraction::MindlinRollingTorsionInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
MindlinRollingTorsionInteraction::MindlinRollingTorsionInteraction(const MindlinRollingTorsionInteraction& p)
        : BaseInteraction(p), MindlinInteraction(p)
{
    rollingSpring_ = p.rollingSpring_;
    torsionSpring_ = p.torsionSpring_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MindlinRollingTorsionInteraction::MindlinRollingTorsionInteraction(const MindlinRollingTorsionInteraction& p) finished"<<std::endl;
#endif
}

/*!
 *
 */
MindlinRollingTorsionInteraction::~MindlinRollingTorsionInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"MindlinRollingTorsionInteraction::~MindlinRollingTorsionInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in/out] os output file stream
 */
void MindlinRollingTorsionInteraction::write(std::ostream& os) const
{
    MindlinInteraction::write(os);
    os << " rollingSpring " << rollingSpring_;
    os << " torsionSpring " << torsionSpring_;
}

/*!
 * \param[in/out] is input file stream
 */
void MindlinRollingTorsionInteraction::read(std::istream& is)
{
    MindlinInteraction::read(is);
    std::string dummy;
    is >> dummy >> rollingSpring_;
    is >> dummy >> torsionSpring_;
}

/*!
 * \details Calls the MindlinInteraction::computeFrictionForce() as well, see MindlinInteraction.cc.
 */
void MindlinRollingTorsionInteraction::computeFrictionForce()
{
    MindlinInteraction::computeFrictionForce();
    
    const MindlinRollingTorsionSpecies* species = getSpecies();
    //If tangential forces are present
    if (getAbsoluteNormalForce() == 0.0) return;
    
    if (species->getRollingFrictionCoefficient() != 0.0)
    {
        Mdouble rollingStiffness = tangentialStiffnessZero_;
        double effectiveRadius = getEffectiveRadius();

        //From Luding 2008, objective rolling velocity (eq 15) w/o 2.0!
        Vec3D rollingRelativeVelocity = - effectiveRadius *
                                        Vec3D::cross(getNormal(),
                                                     getP()->getAngularVelocity() - getI()->getAngularVelocity());

        const Mdouble springLength = rollingSpring_.getLength();
        if (springLength > 1e-10)
        {
            rollingSpring_ -= Vec3D::dot(rollingSpring_, getNormal()) * getNormal();
            rollingSpring_ *= springLength / rollingSpring_.getLength();
            // logger.assert(std::abs(slidingSpring_.getLength() - springLength) < 1e-10, "Spring length not the same after rotation");
        }

        //integrate(getHandler()->timeStep_);
        rollingSpringVelocity_ = rollingRelativeVelocity;
        rollingSpring_ += rollingSpringVelocity_ * getHandler()->getDPMBase()->getTimeStep();
        //logger(INFO,"rollingSpring.normalDirection %",Vec3D::dot(rollingSpring_/rollingSpring_.getLength(),getNormal()));


        //Calculate test force acting on P including viscous force
        Vec3D rollingForce = -rollingStiffness * rollingSpring_ -
                             species->getRollingDissipation() * rollingRelativeVelocity;
        
        //tangential forces are modelled by a spring-damper of elasticity kt and viscosity dispt (sticking),
        //but the force is limited by Coulomb friction (rolling):
        Mdouble rollingForceSquared = rollingForce.getLengthSquared();
        if (rollingForceSquared <=
            mathsFunc::square(species->getRollingFrictionCoefficientStatic() * getAbsoluteNormalForce()))
        {
            //if sticking (|ft|<=mu*|fn|), apply the force
        }
        else
        {
            //if rolling, resize the tangential force such that |ft|=mu*|fn|
            rollingForce *= species->getRollingFrictionCoefficient() * getAbsoluteNormalForce() /
                            std::sqrt(rollingForceSquared);
            //resize the tangential spring accordingly such ft=-k*delta-nu*relVel
            rollingSpring_ = -(rollingForce + species->getRollingDissipation() * rollingRelativeVelocity) /
                             rollingStiffness;
        }
        //Add (virtual) rolling force to torque
        addTorque(effectiveRadius * Vec3D::cross(getNormal(), rollingForce));
    } //end if rolling force
    
    if (species->getTorsionFrictionCoefficient() != 0.0)
    {
        Mdouble torsionStiffness = tangentialStiffnessZero_;
        ///\todo TW: Why do we not use the corrected diameter here, as in the rolling case? And check if Stefan uses radius or diameter
        Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
        
        //From Luding2008, spin velocity (eq 16) w/o 2.0!
        Vec3D torsionRelativeVelocity = effectiveDiameter * Vec3D::dot(getNormal(), getP()->getAngularVelocity() -
                                                                                    getI()->getAngularVelocity()) *
                                        getNormal();
        
        //Integrate the spring
        torsionSpringVelocity_ = torsionRelativeVelocity;
        //integrate(getHandler()->timeStep_);
        torsionSpring_ +=
                Vec3D::dot(torsionSpring_ + torsionSpringVelocity_ * getHandler()->getDPMBase()->getTimeStep(),
                           getNormal()) * getNormal();
        
        //Calculate test force acting on P including viscous force
        Vec3D torsionForce = -torsionStiffness * torsionSpring_ -
                             species->getTorsionDissipation() * torsionRelativeVelocity;
        
        //tangential forces are modelled by a spring-damper of elasticity kt and viscosity dispt (sticking),
        //but the force is limited by Coulomb friction (torsion):
        Mdouble torsionForceSquared = torsionForce.getLengthSquared();
        if (torsionForceSquared <=
            mathsFunc::square(species->getTorsionFrictionCoefficientStatic() * getAbsoluteNormalForce()))
        {
            //if sticking (|ft|<=mu*|fn|), apply the force
        }
        else
        {
            //if torsion, resize the tangential force such that |ft|=mu*|fn|
            torsionForce *= species->getTorsionFrictionCoefficient() * getAbsoluteNormalForce() /
                            std::sqrt(torsionForceSquared);
            //resize the tangential spring accordingly such ft=-k*delta-nu*relVel
            torsionSpring_ = -(torsionForce + species->getTorsionDissipation() * torsionRelativeVelocity) /
                             torsionStiffness;
        }
        //Add (virtual) rolling force to torque
        addTorque(effectiveDiameter * torsionForce);
    } //end if torsion force
}

/*!
 * \param[in] timeStep the amount of time by which the solution is evolved
 */
void MindlinRollingTorsionInteraction::integrate(Mdouble timeStep)
{
    MindlinInteraction::integrate(timeStep);
    rollingSpring_ += rollingSpringVelocity_ * timeStep;
    torsionSpring_ += Vec3D::dot(torsionSpring_ + torsionSpringVelocity_ * timeStep, getNormal()) * getNormal();
}

/*!
 * \return Mdouble the total elastic energy stored in the frictional springs
 */
Mdouble MindlinRollingTorsionInteraction::getElasticEnergy() const
{
    return MindlinInteraction::getElasticEnergy()
           + 0.5 * tangentialStiffnessZero_ * rollingSpring_.getLengthSquared()
           + 0.5 * tangentialStiffnessZero_ * torsionSpring_.getLengthSquared();
}

/*!
 * \return const MindlinRollingTorsionSpecies*
 */
const MindlinRollingTorsionSpecies* MindlinRollingTorsionInteraction::getSpecies() const
{
    return static_cast<const MindlinRollingTorsionSpecies*>(getBaseSpecies()->getFrictionForce());
;
}

/*!
 * \return std::string
 */
std::string MindlinRollingTorsionInteraction::getBaseName() const
{
    return "Friction";
}

/*!
 *
 */
void MindlinRollingTorsionInteraction::reverseHistory()
{
    MindlinInteraction::reverseHistory();
    //rollingSpring_=-rollingSpring_;
    //rollingSpringVelocity_=-rollingSpringVelocity_;
    torsionSpring_ = -torsionSpring_;
    torsionSpringVelocity_ = -torsionSpringVelocity_;
}

/*!
 *
 */
void MindlinRollingTorsionInteraction::rotateHistory(Matrix3D& rotationMatrix)
{
    MindlinInteraction::rotateHistory(rotationMatrix);
    rollingSpring_ = rotationMatrix * rollingSpring_;
    rollingSpringVelocity_ = rotationMatrix * rollingSpringVelocity_;
    torsionSpring_ = rotationMatrix * torsionSpring_;
    torsionSpringVelocity_ = rotationMatrix * torsionSpringVelocity_;
}

Vec3D MindlinRollingTorsionInteraction::getRollingSpring() const
{
    return rollingSpring_;
}

Vec3D MindlinRollingTorsionInteraction::getTorsionSpring() const
{
    return torsionSpring_;
}

void MindlinRollingTorsionInteraction::setRollingSpring(Vec3D rollingSpring)
{
    rollingSpring_ = rollingSpring;
}

void MindlinRollingTorsionInteraction::setTorsionSpring(Vec3D torsionSpring)
{
    torsionSpring_ = torsionSpring;
}





