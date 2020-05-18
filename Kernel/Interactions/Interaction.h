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

#ifndef INTERACTION_H
#define INTERACTION_H

#include "FrictionForceInteractions/EmptyFrictionInteraction.h"
#include "AdhesiveForceInteractions/EmptyAdhesiveInteraction.h"
#include "InteractionHandler.h"
#include "BaseInteractable.h"

//used for mpi
#include "MpiContainer.h"
#include "MPIInteraction.h"
#include "Logger.h"

class BaseInteractable;

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
class Species;

/*!
 * \class Interaction
 * \brief Contains information about the contact between two interactables, 
 * BaseInteraction::P_ and BaseInteraction::I_;.
 * \details An instance is created every time two BaseInteractables (particles 
 * or walls) get into contact with another, and gets deleted when the contact 
 * ends.
 * 
 * Next, the call process is described for this function for the case of two 
 * particles (it is very similar for a particle-wall contact.
 * 
 * Then, in each time step, every time that a contact is detected, 
 * DPMBase::computeInternalForces is called, which calls 
 * BaseParticle::getInteractionWith to create a new Interaction (setting 
 * the links to the interactables 
 * \link BaseInteraction::P_ P_\endlink and 
 * \link BaseInteraction::I_ I_\endlink 
 * and adds it to the InteractionHandler) or find an existing one. When a new 
 * Interaction is created, the Species determines what kind of Interaction it will be. 
 * DPMBase::computeInternalForces also sets the 
 * \link BaseInteraction::timeStamp_ timeStamp_\endlink, 
 * \link BaseInteraction::normal_ normal_\endlink, 
 * \link BaseInteraction::overlap_ overlap_\endlink, 
 * \link BaseInteraction::distance_ distance_\endlink, and 
 * \link BaseInteraction::contactPoint_ contactPoint_\endlink
 * of the Interaction. 
 * 
 * Then, computeForce is called, which sets the force_, which sets BaseInteraction::force_ and 
 * BaseInteraction::torque_ of the Interaction, and a few temporary values to communicate between 
 * the three different computeForce routines (see diamond inheritance below).
 * 
 * The force_ and torque_ is then used in DPMBase::integrateBeforeForceComputation 
 * and DPMBase::integrateAfterForceComputation to calcutate the new positions 
 * and velocities.
 * \todo MX: I do not think the above is correct. integrateBeforeForceComputation() updates the positions
 *
 * As there are many types of contact forces, the class is templated to allow 
 * for different force models. This is done in a diamond inheritance structure: 
 * First, three kinds of Interactions are created:
 * - NormalForceInteraction: Computes the normal contact force and sets the 
 *   torque to zero (we currently only have spherical particles). Also sets
 *   \link BaseInteraction::relativeVelocity_ relativeVelocity_\endlink, 
 *   \link BaseInteraction::normalRelativeVelocity_ normalRelativeVelocity_\endlink, 
 *   \link BaseInteraction::absoluteNormalForce_ absoluteNormalForce_\endlink.
 * - FrictionForceInteraction: Computes tangential contact force and torques.
 * - AdhesiveForceInteraction: Computes the short-range normal contact force.
 *
 * A full Interaction object is then derived by inheriting from all of the above:
 *  \dot
 *  digraph example {
 *      node [shape=record, fontname=Helvetica, fontsize=10];
 *      a [ label="class BaseInteraction" URL="\ref BaseInteraction"];
 *      f [ label="class ParticleInteraction" URL="\ref ParticleInteraction"];
 *      b [ label="class AdhesiveForceInteraction" URL="\ref AdhesiveForceInteraction"];
 *      c [ label="class FrictionForceInteraction" URL="\ref FrictionForceInteraction"];
 *      d [ label="class NormalForceInteraction" URL="\ref NormalForceInteraction"];
 *      e [ label="class Interaction" URL="\ref Interaction"];
 *      a -> b [ arrowhead="open" ];
 *      a -> c [ arrowhead="open" ];
 *      a -> d [ arrowhead="open" ];
 *      a -> f [ arrowhead="open" ];
 *      f -> e [ arrowhead="open" ];
 *      b -> e [ arrowhead="open" ];
 *      c -> e [ arrowhead="open" ];
 *      d -> e [ arrowhead="open" ];
 *  }
 *  \enddot
 * 
 */
//this class combines normal and tangential force laws
template<class NormalForceInteraction, class FrictionForceInteraction=EmptyFrictionInteraction, class AdhesiveForceInteraction=EmptyAdhesiveInteraction>
class Interaction : public NormalForceInteraction, public FrictionForceInteraction, public AdhesiveForceInteraction
{
public:
    
    ///\brief The default constructor.
    Interaction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    ///\brief Empty constructor
    Interaction();
    
    ///\brief The default copy constructor.
    Interaction(const Interaction& p);
    
    ///\brief The default destructor.
    virtual ~Interaction();
    
    ///\brief Creates a copy of this Interaction.
    Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* copy() const final;
    
    ///\brief Computes the normal, tangential, and adhesive forces.
    void computeForce() final;
    
    ///\brief Computes the interaction parameters based on the forces and torques.
    void computeInteraction();
    
    ///\brief Read Interaction properties from a file.
    void read(std::istream& is) final;
    
    ///\brief Writes Interaction properties to a file.
    void write(std::ostream& os) const final;
    
    ///\brief Returns the name of the Interaction.
    std::string getName() const final;
    
    ///\brief Returns the elastic energy stored in the Interaction.
    Mdouble getElasticEnergy() const final;
    
    ///\brief Integrates the time-dependent parameters of the contact force.
    void integrate(Mdouble timeStep) final;
    
    ///\brief Reverses the parameters of the contact force.
    void reverseHistory() final;
    
    void rotateHistory(Matrix3D& rotationMatrix) final;
    
    void actionsAfterTimeStep();
    
    /*!
     * \brief returns the overlap at which the repulsive elastic force equals a given adhesive force; to be implemented in the normal force
     */
    Mdouble getElasticEnergyAtEquilibrium(Mdouble adhesiveForce) const;
    
    void getMPIInteraction(void* historyDataArray, unsigned int index) const final;
    
    //if resetPointers is true, then the pointers in the interaction will be set to nullptr, indicating it is a dummyInteraction
    //and that it is not actually linked to real particles (which might not be in the domain)
    void setMPIInteraction(void* historyDataArray, unsigned int index, const bool resetPointers) final;
    
    void getInteractionDetails(void* interactionDataArray, unsigned int index, unsigned int& identificationP,
                               unsigned int& identificationI, bool& isWallInteraction, unsigned& timeStamp);
    
    void createMPIType() final;
    
    void* createMPIInteractionDataArray(unsigned int numberOfInteractions) const final;
    
    void deleteMPIInteractionDataArray(void* dataArray) final;
};

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::Interaction(
        BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp), NormalForceInteraction(P, I, timeStamp),
          FrictionForceInteraction(P, I, timeStamp), AdhesiveForceInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"Interaction::Interaction() finished"<<std::endl;
#endif
}


template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::Interaction()
        : BaseInteraction(), NormalForceInteraction(), FrictionForceInteraction(), AdhesiveForceInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"Interaction::Interaction() finished"<<std::endl;
#endif
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::Interaction(
        const Interaction& p)
        : BaseInteraction(p), NormalForceInteraction(p), FrictionForceInteraction(p), AdhesiveForceInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"Interaction::Interaction(const Interaction &p finished"<<std::endl;
#endif
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::~Interaction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"Interaction::~Interaction() finished"<<std::endl;
#endif
}

/*! 
 * \details Useful for polymorphism as it can be called from a BaseInteraction* pointer.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>*
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::copy() const
{
    return new Interaction(*this);
}

/*!
 * \details Writes Interaction properties in human-readable form to a file, typically Files::restartFile.
 * \param [out] os the ostream to which the data is written.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::write(std::ostream& os) const
{
    NormalForceInteraction::write(os);
    FrictionForceInteraction::write(os);
    AdhesiveForceInteraction::write(os);
    
}

/*!
 * \details Reads Interaction properties in human-readable form from a file, typically Files::restartFile.
 * \param [in] is the istream from which the data is read.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::read(std::istream& is)
{
    NormalForceInteraction::read(is);
    FrictionForceInteraction::read(is);
    AdhesiveForceInteraction::read(is);
}

/*!
 * \details Called by BaseInteraction::copySwitchPointer to reverse the  
 * parameters of the contact force in the case that the interactables P_ and I_ 
 * are switched. This is needed to deal with periodic particles.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::reverseHistory()
{
    NormalForceInteraction::reverseHistory();
    FrictionForceInteraction::reverseHistory();
    AdhesiveForceInteraction::reverseHistory();
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::rotateHistory(
        Matrix3D& rotationMatrix)
{
    NormalForceInteraction::rotateHistory(rotationMatrix);
    FrictionForceInteraction::rotateHistory(rotationMatrix);
    AdhesiveForceInteraction::rotateHistory(rotationMatrix);
}

/*!
 * \details Returns the name of the Interaction, which depends on the template parameters.
 * \return the string to which the name is written.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
std::string Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::getName() const
{
    return NormalForceInteraction::getBaseName() + FrictionForceInteraction::getBaseName() +
           AdhesiveForceInteraction::getBaseName() + "Interaction";
}

/*!
 * \details Called by ??? to integrate time-dependent parameters of the 
 * contact force, such as the SlidingFrictionInteraction::slidingSpring_
 * \param[in] the time step.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::integrate(Mdouble timeStep)
{
    FrictionForceInteraction::integrate(timeStep);
}

/*!
 * \details Computes the normal, tangential, and adhesive forces (in that order).
 * The order is important, as the normal force computation also calculates some 
 * additional parameters required by the other force laws:
 * \link BaseInteraction::relativeVelocity_ relativeVelocity_\endlink, 
 * \link BaseInteraction::normalRelativeVelocity_ normalRelativeVelocity_\endlink, 
 * \link BaseInteraction::absoluteNormalForce_ absoluteNormalForce_\endlink.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::computeForce()
{
    NormalForceInteraction::computeNormalForce();
    FrictionForceInteraction::computeFrictionForce();
    AdhesiveForceInteraction::computeAdhesionForce();
}

/*!
 * \details Returns the elastic energy stored in the Interaction, adding up 
 * contributions from the normal, frictional and adhesive interaction
 * The elastic (=potential) energy is defined as the energy gained by separating P_ and I_.
 * \return the computed elastic energy.
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Mdouble
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::getElasticEnergy() const
{
    return NormalForceInteraction::getElasticEnergy() + FrictionForceInteraction::getElasticEnergy() +
           AdhesiveForceInteraction::getElasticEnergy();
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::actionsAfterTimeStep()
{
    AdhesiveForceInteraction::actionsAfterTimeStep();
}

/*!
 * \brief returns the overlap at which the repulsive elastic force equals a given adhesive force; to be implemented in the normal force
 */
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
Mdouble
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::getElasticEnergyAtEquilibrium(
        Mdouble adhesiveForce) const
{
    return NormalForceInteraction::getElasticEnergyAtEquilibrium(adhesiveForce);
};


template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::getMPIInteraction(
        void* interactionDataArray, unsigned int index) const
{
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction> interactionData;
    
    interactionData.copyFromInteraction(this);
    
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* array = static_cast<MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>*>(interactionDataArray);
    array[index] = interactionData;
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::setMPIInteraction(
        void* interactionDataArray, unsigned int index, const bool resetPointers)
{
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* array = static_cast<MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>*>(interactionDataArray);
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction> interactionData = array[index];
    
    //Copy the interaction data into the interaction
    interactionData.copyToInteraction(this, resetPointers);
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::getInteractionDetails(
        void* interactionDataArray, unsigned int index, unsigned int& identificationP, unsigned int& identificationI,
        bool& isWallInteraction, unsigned& timeStamp)
{
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* array = static_cast<MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>*>(interactionDataArray);
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction> interactionData = array[index];
    
    identificationP = interactionData.P;
    identificationI = interactionData.I;
    isWallInteraction = interactionData.isWallInteraction;
    timeStamp = interactionData.timeStamp;
}


template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void*
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::createMPIInteractionDataArray(
        unsigned int numberOfInteractions) const
{
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* array = new MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>[numberOfInteractions];
    void* interactionArray = static_cast<void*>(array);
    return interactionArray;
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void
Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::deleteMPIInteractionDataArray(
        void* dataArray)
{
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* array = static_cast<MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>*>(dataArray);
    delete[] array;
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::createMPIType()
{
#ifdef MERCURY_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
    MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction> dummyHistoryInteraction;

    communicator.createMercuryMPIType(dummyHistoryInteraction, MercuryMPIType::INTERACTION);
#endif
}

#endif

