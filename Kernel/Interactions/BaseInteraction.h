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

#ifndef BASEINTERACTION_H
#define BASEINTERACTION_H

#include <vector>
#include "BaseObject.h"
#include "Math/Vector.h"
#include "Math/Matrix.h"
#include "Logger.h"

class InteractionHandler;

class BaseParticle;

class BaseWall;

class BaseSpecies;

class BaseInteractable;


/*!
 * \class BaseInteraction
 * \brief   Stores information about interactions between two interactable 
 *          objects; often particles but could be walls etc. By info about interactions
 *          one means the overlaps, contact point, forces, torques, relative velocities
 *          etc.
 * \details This class is completely a MercuryDPM internal and should not be visible to
 *          users. It is also an abstract base class. All possible Mercury
 *          interactions are derived from this class. It is basically a policy
 *          class that defines what an interaction must have in Mercury.
 *          Advanced uses that need to add a new interaction must derive it from
 *          here and implement the black functions.
 */
class BaseInteraction : public BaseObject
{
public:
    /*!
     * \brief A constructor takes the BaseInteractable objects which are interacting (come into contact)
     *         and time the interaction starts.
     */
    BaseInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    /// \todo: empty baseInteraction, used in mpi code probably
    BaseInteraction();
    
    ///\brief Copy constructor
    BaseInteraction(const BaseInteraction& p);
    
    ///\brief The destructor. It removes this interactions from the objects that were interacting, and writes the time
    ///to a file when needed.
    ~BaseInteraction() override;
    
    /*!\brief If an interaction needs to do something before it gets erased, add it here.
     * E.g. Liquid bridges rupture at the end of their lifetime, and
     * the liquid bridge volume has to be redistributed.
     * The reason this action is not done in the destructor is
     * that this action should not be taken when erasing ghost interactions
     */
    virtual void actionsOnErase()
    {};
    
    /*!
     * \brief Virtual function that contains the force law between the two objects interacting.
     */
    virtual void computeForce();
    
    /*!
     * \brief Interaction read function, which accepts an std::istream as input.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Interaction print function, which accepts an std::ostream as input.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Writes forces data to the FStat file.
     */
    void writeToFStat(std::ostream& os, Mdouble time) const;
    
    /*!
     * \brief Virtual function which allows interactions to be named.
     */
    std::string getName() const override;
    
    /*!
     * \brief Returns a Mdouble which is the current about of Elastic energy in the interaction.
     */
    virtual Mdouble getElasticEnergy() const;
    
    /*!
     * \brief Sets the interaction distance between the two interacting objects. 
     */
    void setDistance(Mdouble distance);
    
    /*!
     * \brief Sets the normal vector between the two interacting objects.
     */
    void setNormal(Vec3D normal);
    
    /*!
     * \brief Set the overlap between the two interacting object.
     */
    void setOverlap(Mdouble overlap);
    
    /*!
     * \brief Set the location of the contact point between the two interacting objects.
     */
    void setContactPoint(Vec3D contactPoint);
    
    /*!
     * \brief Updates the time step of the interacting. Note, time steps used to find completed interactions.
     */
    void setTimeStamp(unsigned timeStamp);
    
    /*!
     * \brief Set the Species of the interaction; note this can either be a Species or MixedSpecies.
     */
    void setSpecies(const BaseSpecies* species);
    
    /*!
     * \brief Sets the first object involved in the interaction (normally a particle).
     */
    void setP(BaseInteractable* P);
    
    /*!
     * \brief Sets the second object involved in the interaction (often particle or wall).
     */
    void setI(BaseInteractable* I);

    /*!
     * \brief Sets the first object involved in the interaction (normally a particle).
     */
    void importP(BaseInteractable *P);

    /*!
     * \brief Sets the second object involved in the interaction (often particle or wall).
     */
    void importI(BaseInteractable *I);
    
    /*!
     * returns the vector from I to P.
     * Note, this cannot simply computed as P-I because I is possibly a ghost particle and might be shifted.
     */
    Vec3D getIP() const;
    
    /*!
     * returns the vector from I to C.
     * Note, this cannot simply computed as P-I because I is possibly a ghost particle and might be shifted.
     */
    Vec3D getIC() const;
    
    /*!
     * returns the vector from C to P.
     * Note, this cannot simply computed as P-I because I is possibly a ghost particle and might be shifted.
     */
    Vec3D getCP() const;
    
    void setLagrangeMultiplier(Mdouble multiplier)
    {
        lagrangeMultiplier_ = multiplier;
    }
    
    Mdouble getLagrangeMultiplier()
    {
        return lagrangeMultiplier_;
    }
    
    /*!
     * \brief Sets the pointer to the interaction hander which is storing this interaction.
     */
    void setHandler(InteractionHandler* handler);
    
    /*!
     * \brief Gets a point to the interaction handlers to which this interaction belongs.
     */
    InteractionHandler* getHandler() const;
    
    /*!
     * \brief Gets the current force (vector) between the two interacting objects.
     * \details Return a reference to a Vec3D which contains the current value of  the force associated with this interaction.
     * \return  A reference to a Vec3D containing the total force.
     */
    const Vec3D& getForce() const
    { return force_; }
    
    /*!
     * \brief Gets the current torque (vector) between the two interacting objects
     * \details Return a reference to a Vec3D which contains the current value of the torque associated with the interaction.
     * \return  A reference to a Vec3D containing the total torque.
     */
    const Vec3D& getTorque() const
    { return torque_; }
    
    /*!
     * \brief Gets the normal vector between the two interacting objects.
     * \details Returns a reference to a Vec3D which contains the current value of the normal associated with the interactions.
     * \return  A reference to a Vec3D containing the current normal.
     */
    const Vec3D& getNormal() const
    { return normal_; }
    
    /*!
     * \brief Gets constant reference to contact point (vector).
     * \details Returns a reference to a Vec3D which contains the current value of the contact point for the particles in interaction. Normally, does not change about the interaction is created.
     * \return  A reference to a Vec3D containing the contact point.
     */
    const Vec3D& getContactPoint() const
    { return contactPoint_; }
    
    /*!
     * \brief Returns a Mdouble with the current overlap between the two interacting objects.
     */
    Mdouble getOverlap() const
    { return overlap_; }
    
    /*!
     * \brief Returns a Mdouble with the current contact between the two interacting objects.
     */
    Mdouble getContactRadius() const;
    
    /*!
     * \brief Removes this interaction from its interaction hander.
     */
    void removeFromHandler();
    
    /*!
     * \brief This copies the interactions of the original particle and replaces the original with the ghost copy.
     */
    void copySwitchPointer(const BaseInteractable* original, BaseInteractable* ghost) const;
    
    /*!
     * \brief 
     * \todo Thomas please document this; as this is the area you are currently rewriting.
     */
    void gatherContactStatistics();
    
    /*!
     * \brief Returns a pointer to first object involved in the interaction (normally a particle).
     * \details Returns a pointer to the first object in the interactions; normally the particle.
     * \return   Pointer to BaseInteraction but normally will be a BaseParticle.
     */
    BaseInteractable* getP()
    {
        logger.assert(P_ != nullptr, "First particle in interaction % is nullptr", getId());
        return P_;
    }
    
    /*!
     * \brief Returns a pointer to the second object involved in the interaction (often a wall or a particle).
     * \details Returns a pointer to the second object in the interaction; often a particle or a wall i.e. a general interactale hence I.
     * \return  Pointer to BaseInteraction often a particle or a wall.
     */
    BaseInteractable* getI()
    {
        logger.assert(I_ != nullptr, "Second particle in interaction % is nullptr", getId());
        return I_;
    }
    
    /*!
     * \brief Returns a constant pointer to the first object involved in the interaction.
     * \details Returns a constant pointer to the first object in the interactions; normally the particle.
     * Why implement const and non-const setters and getters? See https://www.gamedev.net/forums/topic/550112-c-non-const-getter-in-terms-of-const-getter/
     * \return   Constant pointer to BaseInteraction but normally will be a BaseParticle.
     */
    const BaseInteractable* getP() const
    {
        logger.assert(P_ != nullptr, "First particle in interaction % is nullptr", getId());
        return P_;
    }
    
    /*!
     * \brief Returns a constant pointer to the second object involved in the interaction.
     * \details Returns a constant pointer to the second object in the interaction; often a particle or a wall i.e. a general interactale hence I.
     * \return  Constant pointer to BaseInteraction often a particle or a wall.
     */
    const BaseInteractable* getI() const
    {
        logger.assert(I_ != nullptr, "Second particle in interaction % is nullptr", getId());
        return I_;
    }
    
    /*!
     * \brief Returns an Mdouble which is the time stamp of the interaction.
     * \details Returns the current time stamp of interaction. This should be the last time the interaction was active and used to find the end of the interaction.
     * \return  Mdouble which is the last time the interaction was active.
     */
    Mdouble getTimeStamp() const
    { return timeStamp_; }
    
    /*!
     * \brief integrates variables of the interaction which need to be integrate e.g. the tangential overlap.
     */
    virtual void integrate(Mdouble timeStep);
    
    /*!
     * \brief get the length of the current tangential overlap 
     */
    virtual Mdouble getTangentialOverlap() const;
    
    /*!
     * \brief Returns an Mdouble which is the norm (length) of distance vector 
     */
    Mdouble getDistance() const;
    
    /*!
     * \brief Returns a constant reference to a vector of relative velocity
     */
    const Vec3D& getRelativeVelocity() const;
    
    /*!
     * \brief Returns a double which is the norm  (length) of the relative velocity vector.
     */
    Mdouble getNormalRelativeVelocity() const;
    
    /*!
     * \brief Returns the absolute value of the norm (length) of the Normal force vector.
     */
    Mdouble getAbsoluteNormalForce() const;
    
    /*!
     * \brief Makes a copy of the interaction and returns a pointer to the copy.
     */
    virtual BaseInteraction* copy() const = 0;
    
    void setFStatData(std::fstream& fstat, BaseParticle* P, BaseWall* I);
    
    void setFStatData(std::fstream& fstat, BaseParticle* P, BaseParticle* I);
    
    unsigned int getMultiContactIdentifier() const;
    
    void setMultiContactIdentifier(unsigned int multiContactIdentifier_);
    
    /*!
     * \brief When periodic particles are used, some interactions need certain
     * history properties rotated (e.g. tangential springs).
     * This is the function for that.
     */
    virtual void rotateHistory(Matrix3D& rotationMatrix);
    
    virtual void actionsAfterTimeStep();
    
    virtual unsigned getNumberOfFieldsVTK() const;
    
    virtual std::string getTypeVTK(unsigned i) const;
    
    virtual std::string getNameVTK(unsigned i) const;
    
    virtual std::vector<Mdouble> getFieldVTK(unsigned i) const;
    
    /*!
     * \brief set total force (this is used by the normal force, tangential forces are added use addForce)
     */
    void setForce(Vec3D force);

protected:
    
    /*!
     * Returns a constant vector to the Tangential Force
     */
    virtual const Vec3D getTangentialForce() const;
    
    /*!
     * \brief Returns a Mdouble to the effective radius of the interaction. (Not corrected for the overlap)
     */
    Mdouble getEffectiveRadius() const;
    
    /*!
     * \brief Returns a Mdouble to the effective radius of the interaction. (Not corrected for the overlap)
     */
    Mdouble getEffectiveMass() const;
    
    //functions that only the derived Interactions have to know about:
    /*!
     * \brief add an force increment to the total force.
     */
    void addForce(Vec3D force);
    
    /*
     * \brief add a torque increment to the total torque.
     */
    void addTorque(Vec3D torque);
    
    //these functions are only used for normal forces and should be made private by the normal forces:
    
    /*!
     * \brief set the total force (this is used by the normal force, tangential torques are added use addTorque)
     */
    void setTorque(Vec3D torque);
    
    /*!
     * \brief set the relative velocity of the current of the interactions.
     */
    void setRelativeVelocity(Vec3D relativeVelocity);
    
    /*!
     * \brief set the normal component of the relative velocity. 
     */
    void setNormalRelativeVelocity(Mdouble normalRelativeVelocit);
    
    /*!
     * \brief the absolute values of the norm (length) of the normal force
     */
    void setAbsoluteNormalForce(Mdouble absoluteNormalForce);
    
    virtual Mdouble getElasticEnergyAtEquilibrium(Mdouble adhesiveForce) const
    { return 0; }
    
    /*!
     * \brief When periodic particles some interaction need certain history properties reversing. This is the function for that.
     */
    virtual void reverseHistory();
    
    
    /*!
     * \brief Writes information about a interaction to the interaction file
     */
    void writeInteraction(std::ostream& os, bool created) const;

public:
    
    /*!
     * \brief Return a constant point to BaseSpecies of the interaction.
     */
    const BaseSpecies* getBaseSpecies() const;
    
    // create the mpi data type for the history data class
    virtual void createMPIType();
    
    // creates the array for the history data
    virtual void* createMPIInteractionDataArray(unsigned int numberOfInteractions) const;
    
    virtual void deleteMPIInteractionDataArray(void* dataArray);
    
    /// copies the history interactions into the data array
    virtual void getMPIInteraction(void* historyDataArray, unsigned int index) const;
    
    
    virtual void getInteractionDetails(void* interactionDataArray, unsigned int index, unsigned int& identificationP,
                                       unsigned int& identificationI, bool& isWallInteraction, unsigned& timeStamp);
    
    
    //Sets the interactionData from the MPI class into the current thingy
    virtual void setMPIInteraction(void* interactionDataArray, unsigned int index, bool resetPointers);
    
    void
    setBasicMPIInteractionValues(int P, int I, unsigned timeStamp, Vec3D force, Vec3D torque, bool isWallInteraction,
                                 bool resetPointers);
    
    void setIdentificationP(unsigned int identification);
    
    void setIdentificationI(int identification);
    
    void setWallInteraction(bool flag);
    
    unsigned int getIdentificationP();
    
    int getIdentificationI();
    
    bool isWallInteraction();

    // returns true is the interaction is an 'internal' bond.
    // This is used in BaseInteraction::writeToFstat to detect internal bonds;
    // these internal forces are not written to fstat and thus don't contribute to fstatistics.
    virtual bool isBonded() const {return false;}

private:
    
    
    /*!
     * Pointer to the InteractionHander for this interaction.
     */
    InteractionHandler* handler_;
    
    /*!
     * This first of two interactable objects involved in the interaction; normally a particle hence P.
     */
    BaseInteractable* P_;
    
    /*!
     * This is the second of the two interactable object involved in the interaction; often a particle or a wall, hence I because it is general interactable. 
     */
    BaseInteractable* I_;
    
    
    // used in mpi code, which can't deal with the pointers above
    unsigned int identificationP_;
    int identificationI_;
    bool isWallInteraction_;
    
    /*!
     * This is Vec3D which stores the contact point of the interaction.
     */
    Vec3D contactPoint_;
    
    /*!
     * Variables calculated by the normal force routines used in other force routines (but not changed by them).
     */
    Vec3D relativeVelocity_;
    
    /*!
     * Variables calculated by NormalForceInteraction, used in other force routines (but not changed by them)
     */
    Mdouble normalRelativeVelocity_;
    
    /*!
     * Variables calculated by NormalForceInteraction, used by FrictionForceInteraction and AdhesiveForceInteraction (but not changed by them)
     */
    Mdouble absoluteNormalForce_;
    
    /*!
     * Variables calculated by  NormalForceInteraction, used by FrictionForceInteraction and AdhesiveForceInteraction (but not changed by them)
     */
    Mdouble distance_;
    
    /*!
     * First computed by NormalForceInteraction, but added to by FrictionForceInteraction and AdhesiveForceInteraction.
     */
    Vec3D force_;
    
    /*!
     * First computed by NormalForceInteraction, but added to by FrictionForceInteraction and AdhesiveForceInteraction.
     */
    Vec3D torque_;
    
    /*!
     * Mdouble which store the last timeStamp the interaction was active.
     */
    ///\todo TW it would be safer to use the integer time step here, instead of the double
    unsigned timeStamp_;
    
    /*!
     * Vec3D which stores the normal vector (unit vector from I to P) of the interaction.
     */
    Vec3D normal_;
    
    /*!
     * Mdouble which stores the current overlap.
     */
    Mdouble overlap_;
    
    /*!
     * Pointer to the species of the interaction could be a mixed species or a species.
     */
    const BaseSpecies* species_;
    
    /*!
     * Lagrange multiplier used by contact detection of superquadrics.
     */
    Mdouble lagrangeMultiplier_;

/*!
     * Identifies individual contact such that it can be distinguished. E.g. the face number of a triangulated wall.
     */
    unsigned multiContactIdentifier_;
};

#endif
