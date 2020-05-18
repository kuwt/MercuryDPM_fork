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

#ifndef MINDLININTERACTION_H
#define MINDLININTERACTION_H

#include "Interactions/BaseInteraction.h"
#include "Math/Vector.h"

class BaseParticle;

class MindlinSpecies;

class BaseInteractable;

/*!
 * \class MindlinInteraction
 * \brief Computes the forces corresponding to sliding friction.
 */
class MindlinInteraction : public virtual BaseInteraction
{
public:
    /*!
     * \brief An alias name for MindlinSpecies data type.
     */
    typedef MindlinSpecies SpeciesType;
    
    /*!
     * \brief Constructor.
     */
    MindlinInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    /*!
     * \brief Copy constructor.
     */
    MindlinInteraction(const MindlinInteraction& p);
    
    /*!
     * \brief Empty constructor
     */
    MindlinInteraction();
    
    /*!
     * \brief Destructor.
     */
    ~MindlinInteraction() override;
    
    /*!
     * \brief Computes the tangential force generated due to compression in the sliding spring.
     *        Does take into account if the interaction is between particle-particle or particle-wall.
     */
    void computeFrictionForce();
    
    /*!
     * \brief Interaction read function, which accepts an std::istream as input.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Interaction write function, which accepts an std::ostream as input.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Increments the amount of compression in sliding spring.
     */
    void integrate(Mdouble timeStep) override;
    
    /*!
     * \brief Returns the amount of elastic energy stored in sliding spring.
     */
    Mdouble getElasticEnergy() const override;
    
    /*!
     * \brief Returns the amount of tangential overlap which is needed by BaseInteraction::writeToFstat().
     */
    Mdouble getTangentialOverlap() const override;
    
    /*!
     * \brief Returns the type/name of interaction (sliding friction interaction)
     */
    std::string getBaseName() const;
    
    //setters and getters
    /*!
     * \brief Returns the sliding friction force vector.
     */
    const Vec3D getTangentialForce() const override;
    
    //k_edit
    //Returns a scalar value letting us know the direction of the tangential force
    //(i.e. positive or negative) so that accurate values can be output to the fstat file
    const Mdouble getTangentialForceDirection() const;
    
    //k_edit
    /*!
     * \brief Returns the absolute value of the norm (length) of the previous Normal force vector.
     */
    Mdouble getAbsoluteNormalForcePrevious() const;
    
    //k_edit
    /*!
     * \brief allows the previous normal force to be (re)set from external classes
     */
    void setAbsoluteNormalForcePrevious(Mdouble absoluteNormalForcePrevious);
    
    /*!
     * \brief Returns a const pointer of type MindlinSpecies*
     */
    const MindlinSpecies* getSpecies() const;
    
    /*!
     * \brief A useful feature if one wants to return to the initial state of the spring. However, 
     *        reverse history decrements the current state to the state corresponding to previous time step. 
     *        Decrements the value of slidingSpring_.
     */
    void reverseHistory() override;
    
    void rotateHistory(Matrix3D& rotationMatrix) override;
    
    //k_edit
    //Allows the K_t0 parameter (for Mindlin model) to be set
    void setTangentialStiffnessZero(Mdouble newKt0);
    
    //k_edit
    //Allows the K_t0 parameter (for Mindlin model) to be accessed
    Mdouble getTangentialStiffnessZero();
    
    //k_edit
    //...and a similar function for the "normal" K_t
    Mdouble getTangentialStiffness();
    
    //k_edit
    //A quick way to update the K_t0 parameter when necessary...
    void updateTangentialStiffnessZero(Mdouble rad, double shearMod);
    //k_edit
    //...and a second, similar method by which to update K_t0 for loading curves in which
    //multiple steps are required to increment the tangential force
    //void updateTangentialStiffnessZeroSecondStep(Mdouble rad, double shearMod);
    
    //k_edit
    //A way to easily update K_t for all possible cases,
    //i.e. initial loading (under constant normal force)...
    void updateTangentialStiffnessInitial(Mdouble fric);
    
    void updateTangentialStiffnessInitial2(Mdouble fric, Vec3D direction);
    
    //...unloading (under constant normal force)...
    void updateTangentialStiffnessUnloading(Mdouble fric, Vec3D direction);
    
    //...reloading (under constant normal force)...
    void updateTangentialStiffnessReloading(Mdouble fric, Vec3D direction);
    
    //...reloading (i.e. tangential force increasing) and normal force varying (moving from state 1 to 2)...
    void updateTangentialStiffnessReloadingTanUp(Mdouble fric, Vec3D direction);
    
    //...and with increasing normal force and decreasing tangential...
    void updateTangentialStiffnessUnloadingTanDown(Mdouble fric, Vec3D direction);
    
    
    //k_new
    void updateK_t(Mdouble fric, Vec3D direction, bool useTurningPoint, bool isLoading);

protected:
    /*!
     * \brief Stores the amount of sliding spring (\f$\delta\f$) compression from the expression \f$f_t=-k*\delta-\nu*relVel\f$.
     *        Set in the member function integrate(), used in computeFrictionForce().
     */
    Vec3D slidingSpring_;
    //k_edit
    //Introducing a parameter to store the value of "slidingSpring_" (i.e. delta_t) for a previous time step
    Vec3D slidingSpringPrevious_;
    /*!
     * \brief Stores the rate at which the sliding spring compressed or relaxed. Set in the member function
     *        computeFrictionForce() and used in integrate().
     */
    Vec3D slidingSpringVelocity_;
    /*!
     * \brief Computes the tangential force such that \f$|f_t|=\mu*|f_n|\f$. Set and computed in computeFrictionForce().
     */
    Vec3D tangentialForce_;
    //k_edit
    //adding a parameter to store the tangential force from a previous time step
    Vec3D tangentialForcePrevious_;
    //a parameter to give the (scalar) DIRECTION of the force, i.e. "with" or "against" the tangential motion!
    Mdouble tangentialForceDirection_;
    //k_edit
    //A pair of parameters to store the force at turning point, for use in calculating unloading/reloading
    //curves using the Mindlin model (see Di Renzo and Di Maio 2004)
    //"UL" and "LU" are used to distinguish between turning points between unloading and loading and
    //those between loading and unloading
    Vec3D tangentialForceTurningPointLU_;
    Vec3D tangentialForceTurningPointUL_;
    //A pair of variables to temporarily store updated versions of these values for use in multi-step
    //Mindlin calculations
    Vec3D tangentialForceTurningPointLUTemp_;
    Vec3D tangentialForceTurningPointULTemp_;
    //a pair of vectors to store the current tangential displacements at each turning point
    Vec3D tangentialDisplacementTurningPointUL_;
    Vec3D tangentialDisplacementTurningPointLU_;
    
    
    //k_edit
    // the constant K_t0 used to calculate the tangential stiffness (K_t) for the
    //Mindlin model (see, for example, Di Renzo and Di Maio, 2004)
    Mdouble tangentialStiffnessZero_;
    //...and a variable to save its value from the previous time step
    Mdouble tangentialStiffnessZeroPrevious_;
    // the constant K_t as used in Di Maio and Di Renzo
    Mdouble tangentialStiffness_;
    //k_edit
    //a simple flag to let functions know whether loading has been performed before or not
    //(0 if no prior loading, 1 if prior loading)
    bool priorLoadingFlag_;
    
    //k_edit
    //a pair of variables to temporarily store the intermediate values of f_t and delta_t
    //during the multiple-point force calculations necessary when both normal and tangential
    //forces are varying
    Vec3D tangentialForceTemp_;
    Vec3D tangentialDisplacementTemp_;
    
    //k_edit
    //A second set of temporary vectors for when >2 steps are required.
    //although not technically necessary, helps the clarity of the code GREATLY!
    Vec3D tangentialForceTemp2_;
    Vec3D tangentialDisplacementTemp2_;
    
    //k_edit
    //the minimal tangential displacement for which the simple loading assumption can be expected to hold
    //in Herz-Mindlin theory
    Mdouble tangentialDisplacementSL_;
    
    //records the initial tangential velocity such that all subsequent motion can be directly compared to this!
    Vec3D initialTangentialVelocity_;
    
    //k_edit
    //parameter allowing previous values of the normal force to be stored and used in calculations
    Mdouble absoluteNormalForcePrevious_;
    
};

#endif
