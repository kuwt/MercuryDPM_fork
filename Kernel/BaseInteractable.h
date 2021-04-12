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

#ifndef BASEINTERACTABLE_H
#define BASEINTERACTABLE_H

#include <functional>
#include <vector>
#include "BaseObject.h"
#include "Math/Vector.h"
#include "Math/Quaternion.h"

class BaseParticle;

class ParticleSpecies;

class BaseInteraction;

class InteractionHandler;

/*!
 * \class BaseInteractable
 * \brief   Defines the basic properties that a interactable object can have.
 * \details Inherits from class BaseObject (public)
 *          Also it includes a lot of code to deal with interactable objects
 *          that have a prescibed motion.
 *          Most of the code in here is MercuryDPM internal. The only place an
 *          user will interface with this code is for setting the lambda
 *          functions that prescribe the motion of infinite mass particles.
 *          \todo Check prescribed objects have infinite mass.
 */
class BaseInteractable : public BaseObject
{
public:
    /*!
     * \brief Default BaseInteractable constructor. 
     */
    BaseInteractable();
    
    /*!
     * \brief Copy constructor. 
     */
    BaseInteractable(const BaseInteractable& p);
    
    /*!
     * \brief Destructor, it simply destructs the BaseInteractable and all the 
     *        objects it contains.
     */
    ~BaseInteractable() override;
    
    /*!
     * \brief Reads a BaseInteractable from an input stream.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Write a BaseInteractable to an output stream.
     * \param[in] os The output stream to which the BaseInteractable is written.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the index of the species associated with the interactable object.
     * \return  Unsigned int which is the unique index of the species
     */
    unsigned int getIndSpecies() const
    { return indSpecies_; }
    
    /*!
     * \brief Sets the index of the Species of this BaseInteractable.
     * \details This set the species associated with this interactable.
     * This function should not be used and BaseInteractable::setSpecies
     * should be used instead.
     * See also BaseInteractable::setSpecies
     */
    virtual void setIndSpecies(unsigned int indSpecies)
    { indSpecies_ = indSpecies; }
    
    /*!
     * \brief Returns a pointer to the species of this BaseInteractable.
     * \details This function return a ParticleSpecies* for the current interacable.
     * Please note, this is a ParticleSpecies; not, a BaseSpecies as interactables must have physically properties as well.
     * \return  constant ParticleSpecies* pointer to the species storing the physical
     * properties of this interactable.
     */
    const ParticleSpecies* getSpecies() const
    {
        return species_;
    }
    
    /*!
     * \brief Sets the species of this BaseInteractable.
     */
    void setSpecies(const ParticleSpecies* species);
    
    /*!
     * \brief Returns the force on this BaseInteractable.
     * \details Return the current force being to the BaseInteractable.
     * Note, the code works by first computing the forces of each
     * interaction and then it loops over all BaseInteracables applying
     * forces to them from the interactions they are involved in.
     * \return  const Vec3D reference that is the total force applied to this         interactable.
     */
    const Vec3D& getForce() const
    { return force_; }
    
    /*!
     * \brief Returns the torque on this BaseInteractable.
     * \details Return the current torque being to the BaseInteractable.
     * Note, the code works by first computing the forces of each
     * interaction and then it loops over all BaseInteracables applying
     * forces to them from the interactions they are involved in.
     * \return  const Vec3D reference that is the total force applied to this
     * interactable.
     */
    const Vec3D& getTorque() const
    { return torque_; }
    
    /*!
     * \brief Sets the force on this BaseInteractable.
     * \details This sets the force being applied to this interactable.
     * Note, first the code computes all forces in the interactions and
     * then loops over all interactable objects applying the forces from
     * the interactions to the interactables involved in the interaction.
     * \param[in]   force   Vec3D which is the force to be applied.
     */
    void setForce(const Vec3D& force)
    { force_ = force; }
    
    /*!
     * \brief Sets the torque on this BaseInteractable.
     * \details This sets the torque being applied to this interactable.
     * Note, first the code computes all force/torques in the interactions
     * and then loops over all interactable objects applying the torques
     * from the interactions to the interactables involved in the
     * interaction.
     * \param[in]   torque   Vec3D which is the force to be applied.
     */
    void setTorque(const Vec3D& torque)
    { torque_ = torque; }
    
    /*!
     * \brief Adds an amount to the force on this BaseInteractable.
     * \details Incremental version of BaseInteractable::setForce.
     * Also see BaseInteraction::setForce for were this is used.
     * \param[in]   addForce    Vec3D incremental force which is added to the total
     * force of the interactable.
     */
    void addForce(const Vec3D& addForce);
    
    /*!
     * \brief Adds an amount to the torque on this BaseInteractable.
     * \details Incremental version of BaseInteractable::setTorque.
     * Also see BaseInteraction::setTorque for were this is used.
     * \param[in]   addTorque    Vec3D incremental force which is added to the total
     * torque of the interactable.
     */
    void addTorque(const Vec3D& addTorque);
    
    /*
     * This function sets the forces and torques acting on an interactable (particle/wall) to zero at teh beginning of computeOneTimeStep.
     *
     * For omp simulations, it further sets does the first step of summing up all forces (and torques) acting on an interactable.
     * It is step 1 of the add-by-reduction approach:
     *  1. Declare a vector forceOMP_, one element per thread, and set it to zero.
     *  2. When a new force needs to be added to the interactable, don't add it immediately.
     *     Add it to forceOMP_[threadNum] instead
     *  3. Add all forces in forceOMP_ together to force_.
     * This reduction approach allows step 2 to be run in parallel, while step 3 needs to be done in sequence.
     * It is efficient, since step 3 is much quicker than step 2.
     */
    void resetForceTorque(int numberOfOMPthreads);
    
    /*
     * This function does the last step of summing up all forces (and torques) acting on an interactable.
     * It is step 3 of the add-by-reduction approach:
     *  1. Declare a vector forceOMP_, one element per thread, and set it to zero.
     *  2. When a new force needs to be added to the interactable, don't add it immediately.
     *     Add it to forceOMP_[threadNum] instead
     *  3. Add all forces in forceOMP_ together to force_.
     * This reduction approach allows step 2 to be run in parallel, while step 3 needs to be done in sequence.
     * It is efficient, since step 3 is much quicker than step 2.
     */
    void sumForceTorqueOMP();
    
    /*!
     * \brief Returns the position of this BaseInteractable.
     * \details Returns the reference to a Vec3D which contains the position of the
     *          interactionable.
     *          Please note the interpretation of this depends on which
     *          interactable. For particles this is the centre of the particle;
     *          where for walls it is one point of the wall given \f$r.n=p\f$
     * \return  Returns a reference to a Vec3D returns the position of the
     *          interactable.
     */
    const Vec3D& getPosition() const
    { return position_; }
    
    /*!
     * \brief Returns the orientation of this BaseInteractable.
     * \details Returns the reference to a Vec3D which contains the orientation of the
     * interactionable.
     * Please note the interpretation of this depends on which
     * interactable. Please see derived objects for details.
     * \return  Returns a reference to a Vec3D returns the position of the
     * interactable.
     */
    const Quaternion& getOrientation() const
    { return orientation_; }
    
    /*!
     * \brief Sets the position of this BaseInteractable.
     * \details Interpretation depends on which interactable is being considered
     * See also BaseInteractable::getPosistion.
     * \param[in] position  Reference to Vec3D storing the position of the particle.
    */
    void setPosition(const Vec3D& position)
    { position_ = position; }
    
    /*!
     * \brief Sets the orientation of this BaseInteractable by defining the vector that results from the rotation of the (1,0,0) vector.
     */
    void setOrientationViaNormal(Vec3D normal);
    
    /*!
     * \brief Sets the orientation of this BaseInteractable by defining the euler angles.
     */
    void setOrientationViaEuler(Vec3D eulerAngle);
    
    /*!
     * \brief Sets the orientation of this BaseInteractable.
     * \details Interpretation depends on which interactable is being considered
     * See also BaseInteractable::getOrientation.
     *
     * \param[in] orientation  Reference to Vec3D storing the orientation
     * of the particle.
     */
    void setOrientation(const Quaternion& orientation)
    { orientation_ = orientation; }
    
    /*!
     * \brief Moves this BaseInteractable by adding an amount to the position.
     */
    virtual void move(const Vec3D& move);
    
    /*!
     * \brief Rotates this BaseInteractable.
     */
    virtual void rotate(const Vec3D& angularVelocityDt);
    
    /*!
     * \brief Returns a list of interactions which belong to this interactable.
     * \return  An list of pointers to all the interactions which this interacable is involved in.
     */
    const std::vector<BaseInteraction*>& getInteractions() const
    { return interactions_; }
    
    /*!
     * \brief Adds an interaction to this BaseInteractable.
     */
    void addInteraction(BaseInteraction* I);
    
    /*!
     * \brief Removes an interaction from this BaseInteractable.
     */
    bool removeInteraction(BaseInteraction* I);
    
    /*!
     * \brief Copies interactions to this BaseInteractable whenever a periodic 
     *        copy made.
     */
    void copyInteractionsForPeriodicParticles(const BaseInteractable& p);
    
    /*!
     * \brief set the velocity of the BaseInteractable.
     */
    void setVelocity(const Vec3D& velocity);
    
    /*!
     * \brief set the angular velocity of the BaseInteractble.
     */
    void setAngularVelocity(const Vec3D& angularVelocity);
    
    /*!
     * \brief adds an increment to the velocity.
     * \details See also BaseInteractable::setVelocity
     * \param[in] velocity Vec3D containing the velocity increment which to increase
     * the velocity by.
     */
    void addVelocity(const Vec3D& velocity)
    { velocity_ += velocity; }
    
    /*!
     * \brief add an increment to the angular velocity.
     */
    void addAngularVelocity(const Vec3D& angularVelocity);
    
    /*!
     * \brief Returns the velocity of this interactable.
     */
    virtual const Vec3D& getVelocity() const;
    
    /*!
     * \brief Returns the angular velocity of this interactable.
     */
    virtual const Vec3D& getAngularVelocity() const;
    
    /*!
     * \brief Allows the position of an infinite mass interactable to be 
     *        prescribed.
     */
    void setPrescribedPosition(const std::function<Vec3D(double)>& prescribedPosition);
    
    /*!
     * \brief Computes the position from the user defined prescribed position 
     *        function.
     */
    void applyPrescribedPosition(double time);
    
    /*!
     * \brief Allows the velocity of an infinite mass interactable to be 
     *        prescribed.
     */
    void setPrescribedVelocity(const std::function<Vec3D(double)>& prescribedVelocity);
    
    /*!
     * \brief Computes the velocity from the user defined prescribed velocity
     *        function.
     */
    void applyPrescribedVelocity(double time);
    
    /*!
     * \brief Allows the orientation of the infinite mass interactbale to be
     *        prescribed. 
     */
    void setPrescribedOrientation(const std::function<Quaternion(double)>& prescribedOrientation);
    
    /*!
     * \brief Computes the orientation from the user defined prescribed 
     *        orientation function.
     */
    void applyPrescribedOrientation(double time);
    
    /*!
     * \brief Allows the angular velocity of the infinite mass interactable to 
     *         be prescribed.
     */
    void setPrescribedAngularVelocity(const std::function<Vec3D(double)>& prescribedAngularVelocity);
    
    /*!
     * \brief Computes the angular velocity from the user defined prescribed 
     *        angular velocity.
     */
    void applyPrescribedAngularVelocity(double time);
    
    /*!
     * \brief Returns the interaction between this object and a given 
     *        BaseParticle
     * \todo TW make sure this function sets normal, distance, overlap, contact point
     * \todo AT why is this a BaseParticle and not a BaseInteratable.
     */
    virtual BaseInteraction* getInteractionWith(BaseParticle* P, unsigned timeStamp, InteractionHandler* interactionHandler)=0;
    
    /*!
     * \brief Returns the velocity at the contact point, use by many force laws.
     */
    virtual const Vec3D getVelocityAtContact(const Vec3D& contact) const;
    
    /*!
     * \brief This is part of integrate routine for objects with infinite mass
     */
    void integrateBeforeForceComputation(double time, double timeStep);
    
    /*!
     * \brief This is part of the integration routine for objects with infinite
     *        mass.
     */
    void integrateAfterForceComputation(double time, double timeStep);
    
    /*!
     * used to distinguish particles which belong to the flow and fixed particles/walls
     */
    virtual bool isFixed() const =0;

    /*! returns the inverse mass. This value is zero for walls and gets overridden for particles that have finite mass*/
    virtual Mdouble getInvMass() const
    { return 0.0; }

    /*! returns the inverse radius, or curvature, of the surface. This value is zero for walls and gets overridden for particles that have finite radius
     * \todo should be wall-type dependent
     * */
    virtual Mdouble getCurvature(const Vec3D& labFixedCoordinates) const
    { return 0.0; }

private:
    /*!
     * User defined function which if set describes the position of the 
     * interactable
     */
    std::function<Vec3D(double)> prescribedPosition_;
    
    /*!
     * User defined function which if set describes the velocity of the
     * interactable.
     */
    std::function<Vec3D(double)> prescribedVelocity_;
    
    /*!
     * User defined function which if set describes the orientation of the 
     * interactable
     */
    std::function<Quaternion(double)> prescribedOrientation_;
    
    /*!
     * User defined functions which if set describes the angular velocity of the
     * interactable.
     */
    std::function<Vec3D(double)> prescribedAngularVelocity_;
    
    /*!
     * Stores the position of the interactable. Exactly what is stored depends
     * on the type of interactable.
     */
    Vec3D position_;
    
    /*!
     * Stores the orientation of the interactable.  Exactly what is stored 
     * depends on the type of interatable
     */
    Quaternion orientation_;
    
    /*!
     * Store the angular velocity of the interactable.
     */
    Vec3D angularVelocity_;
    
    /*!
     * Stores the force applied to the interactable. 
     */
    Vec3D force_;
    
    /*!
     * Stores the torque applied to the interactable.
     */
    Vec3D torque_;
    
    /*
     * This variable is needed in omp simulations to sum up all forces (and torques) acting on an interactable.
     * It is used in an add-by-reduction approach:
     *  1. Declare a vector forceOMP_, one element per thread, and set it to zero.
     *  2. When a new force needs to be added to the interactable, don't add it immediately.
     *     Add it to forceOMP_[threadNum] instead
     *  3. Add all forces in forceOMP_ together to force_.
     * This reduction approach allows step 2 to be run in parallel, while step 3 needs to be done in sequence.
     * It is efficient, since step 3 is much quicker than step 2.
     */
    std::vector<Vec3D> forceOMP_;
    
    /*!
     * See #BaseInteractable::forceOMP_.
     */
    std::vector<Vec3D> torqueOMP_;

    /*!
     * Point to the ParticlesSpecies which stores density and other material
     * properties of the interactable. 
     */
    const ParticleSpecies* species_;
    
    /*!
     * Stores the index on the species associated with this interactable.
     */
    unsigned int indSpecies_;
    
    /*!
     * Stores the velocity of this interactable.
     */
    Vec3D velocity_;
    
    /*!
     * List of interactions this interactable is involved with. 
     */
    std::vector<BaseInteraction*> interactions_;
};

#endif

