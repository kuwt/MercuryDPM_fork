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

#ifndef SPECIES_H
#define SPECIES_H

#include "MixedSpecies.h"
#include "ParticleSpecies.h"
#include "Species/FrictionForceSpecies/EmptyFrictionSpecies.h"
#include "Species/AdhesiveForceSpecies/EmptyAdhesiveSpecies.h"
#include "Interactions/Interaction.h"

class BaseInteraction;

/*!
 * \class Species
 * \brief Contains material and contact force properties.
 * \details The species class contains both material and contact force properties. 
 * As there are many types of contact forces, the class is templated to allow for different force models.
 * 
 * This is done in a diamond inheritance structure: First, four kinds of Species are created:
 * - ParticleSpecies: Contains material properties like density.
 * - NormalForceSpecies: Describes the normal contact force, therefore contains properties such as stiffness and dissipation.
 * - FrictionForceSpecies: Describes the tangential contact force, therefore contains properties such as SlidingFrictionCoefficient and SlidingDissipation.
 * - AdhesiveForceSpecies: Describes the short-range normal contact force, which are typically adhesive. Therefore, it contains properties such as interactionDistance and adhesiveStiffness.
 * A full Species object is then derived by inheriting from all of the above.
 *  \dot
 *  digraph example {
 *      node [shape=record, fontname=Helvetica, fontsize=10];
 *      a [ label="class BaseSpecies" URL="\ref BaseSpecies"];
 *      f [ label="class ParticleSpecies" URL="\ref ParticleSpecies"];
 *      b [ label="class AdhesiveForceSpecies" URL="\ref AdhesiveForceSpecies"];
 *      c [ label="class FrictionForceSpecies" URL="\ref FrictionForceSpecies"];
 *      d [ label="class NormalForceSpecies" URL="\ref NormalForceSpecies"];
 *      e [ label="class Species" URL="\ref Species"];
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
 * There is only one ParticleSpecies; the three force species can be added using templates. 
 * E.g. to create a species with linear elastic-dissipative normal forces, sliding 
 * friction forces, and linear reversible short-range adhesion forces, use
 * >  Species<LinearViscoelasticNormalSpecies, SlidingFrictionSpecies, LinearReversibleAdhesiveSpecies> species;
 * EmptyFrictionSpecies and EmptyAdhesiveSpecies are the default template parameters, 
 * therefore they don't have to be explicitly named. A species with only linear 
 * elastic-dissipative normal forces and sliding friction tangential forces is created using
 * >  Species<LinearViscoelasticNormalSpecies, SlidingFrictionSpecies> species;
 * Also, to avoid that the user has to use templates, many combination of forces 
 * are predefined as a typedef, with the names of the force laws concatenated, e.g.
 * >  typedef Species<LinearViscoelasticNormalSpecies, SlidingFrictionSpecies> LinearViscoelasticSlidingFrictionSpecies;
 * Therefore to create such a Species, one has to type
 * >  LinearViscoelasticSlidingFrictionSpecies species;
 * 
 * If your code contains multiple Species, then MixedSpecies have to be defined, 
 * which contain the description of the forces acting between two species (it 
 * does not need a ParticleSpecies, since it is never used to assign mass to an object).  
 * Note, we assume that all Species have to be of the same type when multiple 
 * Species are used.
 *  \dot
 *  digraph example {
 *      node [shape=record, fontname=Helvetica, fontsize=10];
 *      a [ label="class BaseSpecies" URL="\ref BaseSpecies"];
 *      b [ label="class AdhesiveForceSpecies" URL="\ref AdhesiveForceSpecies"];
 *      c [ label="class FrictionForceSpecies" URL="\ref FrictionForceSpecies"];
 *      d [ label="class NormalForceSpecies" URL="\ref NormalForceSpecies"];
 *      e [ label="class MixedSpecies" URL="\ref MixedSpecies"];
 *      a -> b [ arrowhead="open" ];
 *      a -> c [ arrowhead="open" ];
 *      a -> d [ arrowhead="open" ];
 *      b -> e [ arrowhead="open" ];
 *      c -> e [ arrowhead="open" ];
 *      d -> e [ arrowhead="open" ];
 *  }
 *  \enddot
 * 
 * Since these MixedSpecies don't require a set of particle properties (such as 
 * density), they are created without a ParticleSpecies. The MixedSpecies is
 * predefined as a typedef inside the Species definition, and can be accessed using
 * >  LinearViscoelasticSlidingFrictionSpecies::MixedSpeciesType mixedSpecies;
 * 
 * When a contact between two particles is determined, an Interaction object is 
 * created by Species::getNewInteraction, as the type of Interaction depends on 
 * the Species type. See Interaction for details. 
 * 
 * \class NormalForceSpecies 
 * \brief Defines a contact force parallel to the contact normal.
 * \details Has to contain the member functions read, write, getBaseName, mix.
 * 
 * \class FrictionForceSpecies 
 * \brief Defines a contact force orthogonal to the contact normal.
 * \details Has to contain the member function getUseAngularDOFs.
 * Also requires read, write, getBaseName, mix.
 * 
 * \class AdhesiveForceSpecies 
 * \brief Defines a short-range (non-contact) force parallel to the contact normal, usually adhesive.
 * \details Has to contain the member function getInteractionDistance.
 * Also requires read, write, getBaseName, mix.
 */
template<class NormalForceSpecies, class FrictionForceSpecies = EmptyFrictionSpecies, class AdhesiveForceSpecies = EmptyAdhesiveSpecies>
class Species final
        : public ParticleSpecies, public NormalForceSpecies, public FrictionForceSpecies, public AdhesiveForceSpecies
{
public:
    ///\brief The correct MixedSpecies type for this Species
    /*!
     * The mixed species between two objects of type 
     * > Species<NormalForceSpecies,FrictionForceSpecies,AdhesiveForceSpecies>
     * is 
     * > MixedSpecies<NormalForceSpecies,FrictionForceSpecies,AdhesiveForceSpecies>
     * To access this MixedSpecies type, you can use this typedef.
     */
    typedef MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies> MixedSpeciesType;
    typedef Interaction<typename NormalForceSpecies::InteractionType, typename FrictionForceSpecies::InteractionType, typename AdhesiveForceSpecies::InteractionType> InteractionType;
    
    
    ///\brief The default constructor.
    Species();
    
    ///\brief The default copy constructor.
    Species(const Species& s);
    
    ///\brief The default destructor.
    virtual ~Species();
    
    /*!
     * \brief Creates a deep copy of the Species from which it is called.
     */
    Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>* copy() const final;
    
    /*!
     * Copies the content of this into the species bs, if they are of the same type.
     */
    void copyInto(BaseSpecies* bs) const final;
    
    /*!
     * \brief Creates a new MixedSpecies with the same force properties as the 
     * Species from which it is called.
     */
    BaseSpecies* copyMixed() const final;
    
    /// Called by SpeciesHandler::readAndAddObject
    void read(std::istream& is) final;
    
    /// \brief Writes the Species properties to an output stream.
    void write(std::ostream& os) const final;
    
    ///\brief Returns the name of the Species as it is used in the restart file. 
    std::string getName() const final;
    
    /*!
     * \brief When a contact between two particles is determined, an Interaction 
     * object is created, as the type of Interaction depends on the Species type.
     */
    BaseInteraction*
    getNewInteraction(BaseInteractable* const P, BaseInteractable* const I, unsigned timeStamp) const final;
    
    //used to create a dummy for MPI purposes (I need a prototype of the interaction)
    BaseInteraction* getEmptyInteraction() const final;
    
    void deleteEmptyInteraction(BaseInteraction* interaction) const final;
    
    /*!
     * \brief Returns true if torques have to be calculated.
     */
    bool getUseAngularDOFs() const final;
    
    /*!
     * \brief This function should not be called.
     */
    void mixAll(BaseSpecies* const S, BaseSpecies* const T) final;

    void actionsAfterTimeStep(BaseParticle* particle) const final {
        NormalForceSpecies::actionsAfterTimeStep(particle);
    }
};

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::Species()
        : ParticleSpecies(this,this,this), NormalForceSpecies(), FrictionForceSpecies(), AdhesiveForceSpecies()
{
    normalForce_ = this;
    frictionForce_ = this;
    adhesiveForce_ = this;
    normalForce_->setBaseSpecies(this);
    frictionForce_->setBaseSpecies(this);
    adhesiveForce_->setBaseSpecies(this);
    logger(DEBUG, "Species::Species() finished");
}

/*!
 * \param[in] s the species that is copied
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::Species(const Species& s)
        : ParticleSpecies(s), NormalForceSpecies(s), FrictionForceSpecies(s), AdhesiveForceSpecies(s)
{
    normalForce_ = this;
    frictionForce_ = this;
    adhesiveForce_ = this;
    normalForce_->setBaseSpecies(this);
    frictionForce_->setBaseSpecies(this);
    adhesiveForce_->setBaseSpecies(this);
    logger(DEBUG, "Species::Species(const Species &p) finished");
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::~Species()
{
    logger(DEBUG, "Species::~Species() finished");
}

/*!
 * \details Creates a deep copy of the Species. As this depends on the template 
 * parameters of Species, the definition of this function is in the Species class.  
 * 
 * To create a copy of a Species, use
 * > species->copy();
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>*
Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::copy() const
{
    return new Species(*this);
}

///Species copy method. It calls to copy constructor of this Species, useful for polymorphism
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
BaseSpecies* Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::copyMixed() const
{
    return new MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>(*this);
}

/*!
 * Useful for polymorphism:
 *    speciesHandler.getObject(i)->copyInto(bs);
 * creates a deep copy (i.e. also copies properties of the derived species), whereas
 *    bs = speciesHandler.getObject(i);
 * would only create a shallow copy.
 *
 * Note: this also copies the index and id. Is that a good idea or not? TW
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::copyInto(BaseSpecies* bs) const
{
    if (bs == nullptr)
    {
        logger(WARN, "Error in %::copyInto: cannot copy into a nullptr");
        return;
    }
    auto* s = dynamic_cast<Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>*>(bs);
    if (s == nullptr)
    {
        logger(WARN, "Error in %::copyInto: copying of species failed", getName());
        return;
    }
    *s = *this;
}

/*!
 * \details It prints human readable Species information to the output stream, typically to Files::restartFile::fstream_
 * The basic species information is written in ParticleSpecies::write; 
 * then the three force types write additional information to the stream.
 * \param[out] os output stream (typically the restart file)
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::write(std::ostream& os) const
{
    ParticleSpecies::write(os);
    NormalForceSpecies::write(os);
    FrictionForceSpecies::write(os);
    AdhesiveForceSpecies::write(os);
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::read(std::istream& is)
{
    ParticleSpecies::read(is);
    NormalForceSpecies::read(is);
    FrictionForceSpecies::read(is);
    AdhesiveForceSpecies::read(is);
    // ensure that interaction distance is recomputed after restarting
    AdhesiveForceSpecies::setInteractionDistance();
}

/*!
 * \details Returns the name of the Species as it is used in the restart file. 
 * The name of the species is a concatenation of the names of the three force  
 * components, e.g.
 * > Species<LinearViscoelasticNormalSpecies,SlidingFrictionSpecies,ReversibleAdhesiveSpecies> species;
 * > std::cout << species.getName();
 * will output "LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies".
 * The EmptyFrictionSpecies and the EmptyAdhesiveSpecies return empty strings, such that 
 * > Species<LinearViscoelasticNormalSpecies> species;
 * > std::cout << species.getName();
 * will output "LinearViscoelasticSpecies".
 * 
 * \return The name of the Species.
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
std::string Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getName() const
{
    return NormalForceSpecies::getBaseName()
           + FrictionForceSpecies::getBaseName()
           + AdhesiveForceSpecies::getBaseName() + "Species";
}

/*!
 * \details The input parameters of this function are directly passed into the constructor for the new interaction.
 * See Interaction for details.
 * \param[in] P first of the two objects that interact
 * \param[in] I second of the two objects that interact
 * \param[in] timeStamp current value of DPMBase::time_
 * \return pointer to the newly created Interaction.
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
BaseInteraction*
Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getNewInteraction(BaseInteractable* const P,
                                                                                           BaseInteractable* const I,
                                                                                           unsigned timeStamp) const
{
    return new InteractionType(P, I, timeStamp);
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
BaseInteraction* Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getEmptyInteraction() const
{
    return new Interaction<typename NormalForceSpecies::InteractionType, typename FrictionForceSpecies::InteractionType, typename AdhesiveForceSpecies::InteractionType>();
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::deleteEmptyInteraction(
        BaseInteraction* interaction) const
{
    Interaction<typename NormalForceSpecies::InteractionType, typename FrictionForceSpecies::InteractionType, typename AdhesiveForceSpecies::InteractionType>* interactionDestroyer;
    interactionDestroyer = dynamic_cast<Interaction<typename NormalForceSpecies::InteractionType, typename FrictionForceSpecies::InteractionType, typename AdhesiveForceSpecies::InteractionType>*>(interaction);
    delete interactionDestroyer;
}

/*!
 * \details Returns true for any FrictionForceSpecies except EmptyFrictionSpecies, 
 * because for spherical particles, torques are only caused by tangential forces. 
 * See SpeciesHandler::useAngularDOFs for more details
 * \return true iff torques have to be calculated
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
bool Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getUseAngularDOFs() const
{
    return FrictionForceSpecies::getUseAngularDOFs();
}

/*!
 * \details This function should not be called as the mix of two Species results 
 * in a MixedSpecies, not a Species. However, I don't know how to avoid creating 
 * this function, as BaseSpecies pointers are used for the MixedSpecies, so 
 * BaseSpecies must have this function defined.
 * \param[in] S the first  of the two species whose properties are mixed to create the new species
 * \param[in] T the second of the two species whose properties are mixed to create the new species
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::mixAll(BaseSpecies* const S UNUSED,
                                                                                     BaseSpecies* const T UNUSED)
{
    logger(ERROR, "%::mix() This function should not be called", getName());
}

#endif
