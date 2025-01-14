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

#ifndef BASESPECIES_H
#define BASESPECIES_H

#include "BaseObject.h"
#include "Math/ExtendedMath.h"
#include "Math/Helpers.h"
#include "Species/NormalForceSpecies/BaseNormalForce.h"
#include "Species/FrictionForceSpecies/BaseFrictionForce.h"
#include "Species/AdhesiveForceSpecies/BaseAdhesiveForce.h"

class SpeciesHandler;// derived from BaseHandler<ParticleSpecies>
//class BaseParticle; //
class BaseInteractable;

class BaseInteraction;

/*!
 * \brief BaseSpecies is the class from which all other species are derived.
 * \details A base species is almost empty, with the exception of a pointer to 
 * the particleHandler. All other properties get added in the inherited classes.
 * See Species for more details.
 */
//Note the getVelocity can for some Species be dependent on which point on the Species is meant.
class BaseSpecies : public BaseObject
{
public:
    
    ///\brief The default constructor.
    BaseSpecies();

    BaseSpecies(BaseNormalForce* normalForce, BaseFrictionForce* frictionForce_, BaseAdhesiveForce* adhesiveForce);

    ///\brief The copy constructor.
    BaseSpecies(const BaseSpecies& p);
    
    ///\brief The default destructor.
    ~BaseSpecies();
    
    /*!
     * \brief Creates a deep copy of the object from which it is called.
     * \details Creates a deep copy of the Species, or MixedSpecies from which 
     * it is called. As this depends on the template parameters of Species, the 
     * definition of this function is in the Species class.  
     * It is defined as a virtual function here to allow the function 
     * to be called from a BaseSpecies pointer (which is the kind of pointer 
     * used for MixedSpecies).
     * 
     * To create a copy of a MixedSpecies, use
     * > BaseSpecies* mixedSpecies->copy();
     * To create a copy of a Species, use
     * > ParticleSpecies* species->copy();
     */
    virtual BaseSpecies* copy() const = 0;
    
    /*!
     * Copies the content of this into the species bs, if they are of the same type.
     */
    virtual void copyInto(BaseSpecies* s) const =0;
    
    ///\brief Sets the pointer to the handler to which this species belongs.
    void setHandler(SpeciesHandler* handler);
    
    ///\brief Returns the pointer to the handler to which this species belongs.
    SpeciesHandler* getHandler() const;
    
    ///\brief Returns the harmonic mean of two variables.
    static Mdouble average(Mdouble a, Mdouble b);

    ///\brief Returns the harmonic mean of two variables, returning inf if either is inf.
    static Mdouble averageInf(Mdouble a, Mdouble b);
    
    ///\brief creates default values for mixed species
    /*!
     * \details returns the largest separation distance (negative overlap) at which 
     * (adhesive) short-range forces can occur (needed for contact detection).
     * Defined in each of the AdhesiveForceSpecies
     * It is defined as a virtual function here to allow the function 
     * to be called from a BaseSpecies pointer (which is the kind of pointer 
     * used for MixedSpecies).
     */
    virtual void mixAll(BaseSpecies* S, BaseSpecies* T) = 0;
    
//setters and getters
    
    ///\brief Returns true if torques (i.e. angular degrees of freedom) have to be calculated.
    /*!
     * \details returns true if torques have to be calculated. This is currently 
     * only true if a tangential force is applied, but can be true for 
     * non-spherical objects even if no tangential force is applied .
     * Defined in FrictionForceSpecies.
     * It is defined as a virtual function here to allow the function 
     * to be called from a BaseSpecies pointer (which is the kind of pointer 
     * used for MixedSpecies).
     */
    virtual bool getUseAngularDOFs() const = 0;
    
    ///\brief returns new Interaction object.
    /*!
     * \details As each Species has its own Interaction type, getNewInteraction 
     * can be used to access the right type of Interaction
     * It is defined as a virtual function here to allow the function 
     * to be called from a BaseSpecies pointer (which is the kind of pointer 
     * used for MixedSpecies).
     */
    virtual BaseInteraction* getNewInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp) const = 0;
    
    //TODO add documentation 
    virtual BaseInteraction* getEmptyInteraction() const = 0;
    
    virtual void deleteEmptyInteraction(BaseInteraction* interaction) const = 0;

    ///\brief returns the largest separation distance at which adhesive short-range forces can occur.
    /*!
     * \details returns the largest separation distance (negative overlap) at which
     * (adhesive) short-range forces can occur (needed for contact detection).
     * Defined in each of the AdhesiveForceSpecies
     * It is defined as a virtual function here to allow the function
     * to be called from a BaseSpecies pointer (which is the kind of pointer
     * used for MixedSpecies).
     */
    Mdouble getInteractionDistance() const {return interactionDistance_;}

    BaseNormalForce* getNormalForce() const {return normalForce_;}

    BaseFrictionForce* getFrictionForce() const {return frictionForce_;}

    BaseAdhesiveForce* getAdhesiveForce() const {return adhesiveForce_;}

    void setInteractionDistance(Mdouble interactionDistance);

protected:

    /*!
     * \brief Sets the boolean constantRestitution_.
     */
    void write(std::ostream& os) const override;

    /*!
     * \param[in] is input stream (typically the restart file)
     */
    void read(std::istream& is) override;

    /*!
     * \brief A pointer to the normal force parameters
     * \details This pointer is used by the Interaction's to get a pointer to the species
     * The pointer is set in the constructors of SPecies and MixedSpecies
     */
    BaseNormalForce* normalForce_;

    /*!
     * \brief A pointer to the friction force parameters
     * \details This pointer is used by the Interaction's to get a pointer to the species
     * The pointer is set in the constructors of SPecies and MixedSpecies
     */
    BaseFrictionForce* frictionForce_;

    /*!
     * \brief A pointer to the adhesive force parameters
     * \details This pointer is used by the Interaction's to get a pointer to the species
     * The pointer is set in the constructors of SPecies and MixedSpecies
     */
    BaseAdhesiveForce* adhesiveForce_;

private:
    /*!
     * \brief A pointer to the handler to which this species belongs. It is 
     * initialized to nullptr and gets set when SpeciesHandler::addObject() 
     * is called. 
     */
    SpeciesHandler* handler_;

    /**
     * Returns the distance between particles of this species below which adhesive forces can occur (needed for contact detection)
     * set by the adhesive species
     */
    Mdouble interactionDistance_;
};

#endif
