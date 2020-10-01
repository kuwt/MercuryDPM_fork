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

#ifndef MIXEDSPECIES_H
#define MIXEDSPECIES_H

#include "Species/FrictionForceSpecies/EmptyFrictionSpecies.h"
#include "Species/AdhesiveForceSpecies/EmptyAdhesiveSpecies.h"
#include "Interactions/Interaction.h"

class BaseInteraction;

/*!
 * \class MixedSpecies
 * \brief Contains contact force properties for contacts between 
 * particles with two different species.
 * \details See Species for details. 
*/
template<class NormalForceSpecies, class FrictionForceSpecies = EmptyFrictionSpecies, class AdhesiveForceSpecies = EmptyAdhesiveSpecies>
class MixedSpecies final : public BaseSpecies, public NormalForceSpecies, public FrictionForceSpecies, public AdhesiveForceSpecies
{
public:
    
    ///\brief The default constructor.
    MixedSpecies();
    
    ///\brief The default copy constructor.
    MixedSpecies(const MixedSpecies& s);
    
    ///\brief Creates a mixed species with the same force properties as a Species.
    MixedSpecies(const Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>& s);
    
    ///\brief The default destructor.
    virtual ~MixedSpecies();
    
    /*!
     * \brief Creates a deep copy of the MixedSpecies from which it is called.
     */
    MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>* copy() const final;
    
    /*!
     * \brief Copies the content of this into the species bs, if they are of the same type.
     */
    void copyInto(BaseSpecies* bs) const final;
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is) final;
    
    /// \brief Writes the MixedSpecies properties to an output stream.
    void write(std::ostream& os) const final;
    
    ///\brief Returns the name of the MixedSpecies as it is used in the restart file. 
    std::string getName() const final;
    
    /*!
     * \brief When a contact between two particles is determined, an Interaction 
     * object is created, as the type of Interaction depends on the MixedSpecies type.
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
     * \brief sets the MixedSpecies properties by mixing the properties of  
     * two particle species
     */
    void mixAll(BaseSpecies* const S, BaseSpecies* const T) final;
};

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::MixedSpecies()
        : BaseSpecies(this,this,this), NormalForceSpecies(), FrictionForceSpecies(), AdhesiveForceSpecies()
{
    logger(DEBUG, "MixedSpecies::MixedSpecies() finished");
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::MixedSpecies(const MixedSpecies& s)
        : BaseSpecies(s), NormalForceSpecies(s), FrictionForceSpecies(s), AdhesiveForceSpecies(s)
{
    normalForce_ = this;
    frictionForce_ = this;
    adhesiveForce_ = this;
    normalForce_->setBaseSpecies(this);
    frictionForce_->setBaseSpecies(this);
    adhesiveForce_->setBaseSpecies(this);
    logger(DEBUG, "MixedSpecies::MixedSpecies(const MixedSpecies &p) finished");
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::MixedSpecies(
        const Species<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>& s)
        : BaseSpecies(s), NormalForceSpecies(s), FrictionForceSpecies(s), AdhesiveForceSpecies(s)
{
    normalForce_ = this;
    frictionForce_ = this;
    adhesiveForce_ = this;
    normalForce_->setBaseSpecies(this);
    frictionForce_->setBaseSpecies(this);
    adhesiveForce_->setBaseSpecies(this);
    logger(DEBUG, "MixedSpecies::MixedSpecies(const MixedSpecies &p) finished");
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::~MixedSpecies()
{
    logger(DEBUG, "MixedSpecies::~MixedSpecies() finished");
}

///MixedSpecies copy method. It calls to copy constructor of this MixedSpecies, useful for polymorphism
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>*
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::copy() const
{
    return new MixedSpecies(*this);
}

/*!
 * Useful for polymorphism:
 *    speciesHandler.getObject(i)->copyInto(bs);
 * creates a deep copy (i.e. also copies properties of the derived species), whereas
 *    bs = speciesHandler.getObject(i);
 * would only create a shallow copy.
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::copyInto(BaseSpecies* bs) const
{
    if (bs == nullptr)
    {
        logger(WARN, "Error in %::copyInto: cannot copy into a nullptr");
        return;
    }
    auto* s = dynamic_cast<MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>*>(bs);
    if (s == nullptr)
    {
        logger(WARN, "Error in %::copyInto: copying of % failed", getName(), s->getName());
        return;
    }
    *s = *this;
}

/*!
 * \details It prints human readable MixedSpecies information to the output stream, 
 * typically to Files::restartFile::fstream_.
 * The basic species information is written in ParticleSpecies::write; 
 * then the three force types write additional information to the stream.
 * \param[out] os output stream (typically the restart file)
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::write(std::ostream& os) const
{
    os << getName();
    os << " idA " << BaseObject::getId();
    os << " idB " << BaseObject::getIndex();
    BaseSpecies::write(os);
    NormalForceSpecies::write(os);
    FrictionForceSpecies::write(os);
    AdhesiveForceSpecies::write(os);
}

/*!
 * Called by SpeciesHandler::readAndAddObject
 * \param[in] is input stream (typically the restart file)
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>
::read(std::istream& is)
{
    //note: name is already read by SpeciesHandler::readAndAddObject
    std::string dummy;
    unsigned int id, index;
    is >> dummy >> id;
    is >> dummy >> index;
    BaseObject::setId(id);
    BaseObject::setIndex(index);
    BaseSpecies::read(is);
    NormalForceSpecies::read(is);
    FrictionForceSpecies::read(is);
    AdhesiveForceSpecies::read(is);
}

/*!
 * \details Returns the name of the MixedSpecies as it is used in the restart file. 
 * The name of the species is a concatenation of the names of the three force  
 * components, e.g.
 * > MixedSpecies<LinearViscoelasticNormalSpecies,SlidingFrictionSpecies,ReversibleAdhesiveSpecies> species;
 * > std::cout << species.getName();
 * will output "LinearViscoelasticSlidingFrictionReversibleAdhesiveMixedSpecies".
 * The EmptyFrictionSpecies and the EmptyAdhesiveSpecies return empty strings, such that 
 * > MixedSpecies<LinearViscoelasticNormalSpecies> species;
 * > std::cout << species.getName();
 * will output "LinearViscoelasticMixedSpecies".
 * 
 * \return The name of the MixedSpecies.
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
std::string MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getName() const
{
    return NormalForceSpecies::getBaseName()
           + FrictionForceSpecies::getBaseName()
           + AdhesiveForceSpecies::getBaseName() + "MixedSpecies";
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
BaseInteraction* MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getNewInteraction(
        BaseInteractable* const P, BaseInteractable* const I, unsigned timeStamp) const
{
    // JMFT: memory is allocated here, so it will need to be freed later
    return new Interaction<typename NormalForceSpecies::InteractionType, typename FrictionForceSpecies::InteractionType, typename AdhesiveForceSpecies::InteractionType>(
            P, I, timeStamp);
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
BaseInteraction*
MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getEmptyInteraction() const
{
    return new Interaction<typename NormalForceSpecies::InteractionType, typename FrictionForceSpecies::InteractionType, typename AdhesiveForceSpecies::InteractionType>();
}

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::deleteEmptyInteraction(
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
bool MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::getUseAngularDOFs() const
{
    return FrictionForceSpecies::getUseAngularDOFs();
}

/*!
 * \details Uses the harmonic mean for most properties. Calls the mix function 
 * for each of the force species from which MixedSpecies is derived.
 * \param[in] S the first  of two species whose properties are mixed to create the new species
 * \param[in] T the second of two species whose properties are mixed to create the new species
 */
template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
void MixedSpecies<NormalForceSpecies, FrictionForceSpecies, AdhesiveForceSpecies>::mixAll(BaseSpecies* const S,
                                                                                          BaseSpecies* const T)
{
    logger.assert_always(T!= nullptr && S!= nullptr,"Arguments of mixAll cannot be null pointers");

    logger.assert_always(S->getNormalForce()->getConstantRestitution() == T->getNormalForce()->getConstantRestitution(), "mixing two LinearPlasticViscoelasticNormalSpecies, but only one has constantRestitution");
    NormalForceSpecies::setConstantRestitution(S->getNormalForce()->getConstantRestitution());

    const auto TN = dynamic_cast<NormalForceSpecies*> (T);
    const auto TF = dynamic_cast<FrictionForceSpecies*> (T);
    const auto TA = dynamic_cast<AdhesiveForceSpecies*> (T);
    logger.assert_always(TN!= nullptr && TF!= nullptr && TA!= nullptr,
            "Cannot mix two species of different type (% and %)",S->getName(),T->getName());

    const auto SN = dynamic_cast<NormalForceSpecies*> (S);
    const auto SF = dynamic_cast<FrictionForceSpecies*> (S);
    const auto SA = dynamic_cast<AdhesiveForceSpecies*> (S);
    logger.assert_always(SN!= nullptr && SF!= nullptr && SA!= nullptr,
            "Cannot mix two species of different type (% and %)",S->getName(),T->getName());

    NormalForceSpecies::mix(SN,TN);
    FrictionForceSpecies::mix(SF,TF);
    AdhesiveForceSpecies::mix(SA,TA);
}


#endif
