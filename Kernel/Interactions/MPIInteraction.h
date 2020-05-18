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

#ifndef MPIINTERACTION_H
#define MPIINTERACTION_H

#include "Interactions/NormalForceInteractions/LinearPlasticViscoelasticInteraction.h"
#include "Interactions/FrictionForceInteractions/SlidingFrictionInteraction.h"
#include "Interactions/FrictionForceInteractions/FrictionInteraction.h"
#include "Interactions/AdhesiveForceInteractions/LiquidBridgeWilletInteraction.h"
#include "Interactions/AdhesiveForceInteractions/LiquidMigrationWilletInteraction.h"
#include "Interactions/AdhesiveForceInteractions/IrreversibleAdhesiveInteraction.h"
#include "Interactions/AdhesiveForceInteractions/BondedInteraction.h"
#include "Logger.h"
#include "MpiDataClass.h"

template<class NormalForceSpecies, class FrictionForceSpecies, class AdhesiveForceSpecies>
class Interaction;

//A monster of a template class, but it is rather universal for all interactions, saves a lot of typing
//In terms of memory it is not as efficient as declaring a class for every interaction type, this due to extra padding concerning
//the empty class. the extra memory used in the transmission is less than 8 bytes per interaction.
//Below all is used to create the MPIHistoryInteraction class which is a lean version to the actual interaction, containing only the history parameters
template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
class MPIInteraction
{
public:
    MPIInteraction();
    
    unsigned int P;
    unsigned int I;
    unsigned int speciesId;
    bool isWallInteraction;
    unsigned timeStamp;
    Vec3D force;
    Vec3D torque;
    
    
    typename std::conditional<
            //If true, enable sliding
            std::is_base_of<SlidingFrictionInteraction, FrictionForceInteraction>::value, Vec3D, Empty>::type
            slidingSpring;
    
    //Friction
    typename std::conditional<
            //if true, enable rolling
            std::is_base_of<FrictionInteraction, FrictionForceInteraction>::value, Vec3D, Empty>::type rollingSpring;
    
    typename std::conditional<
            //if true, enable torsion
            std::is_base_of<FrictionInteraction, FrictionForceInteraction>::value, Vec3D, Empty>::type torsionSpring;
    
    //Contact
    typename std::conditional<
            //if true, enablel wasInContact
            (std::is_base_of<LiquidMigrationWilletInteraction, AdhesiveForceInteraction>::value
             || std::is_base_of<LiquidBridgeWilletInteraction, AdhesiveForceInteraction>::value
             ||
             std::is_base_of<IrreversibleAdhesiveInteraction, AdhesiveForceInteraction>::value), bool, Empty>::type wasInContact;
    
    //Bonded
    typename std::conditional<
            //if true, enable bonded
            std::is_base_of<BondedInteraction, AdhesiveForceInteraction>::value, bool, Empty>::type bonded;
    
    //Liquidbridge
    typename std::conditional<
            //if true, enable liquidbridgeVolume
            std::is_base_of<LiquidMigrationWilletInteraction, AdhesiveForceInteraction>::value, Mdouble, Empty>::type liquidbridgeVolume;

    //MaxOverlap
    typename std::conditional<
        //if true, enable max overlap
        std::is_base_of<LinearPlasticViscoelasticInteraction, NormalForceInteraction>::value, Mdouble, Empty>::type maxOverlap;    
    
    void copyFromInteraction(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction);
    
    void copyToInteraction(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction,
            const bool resetPointers);
    
    //Sliding
    template<class DUMMY= FrictionForceInteraction>
    typename std::enable_if<std::is_base_of<SlidingFrictionInteraction, DUMMY>::value, void>::type
    getSlidingSpring(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        slidingSpring = interaction->FrictionForceInteraction::getSlidingSpring();
    }
    
    //Sliding
    template<class DUMMY= FrictionForceInteraction>
    typename std::enable_if<std::is_base_of<SlidingFrictionInteraction, DUMMY>::value, void>::type
    setSlidingSpring(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        interaction->FrictionForceInteraction::setSlidingSpring(slidingSpring);
    }
    
    //No sliding
    template<class DUMMY = FrictionForceInteraction>
    typename std::enable_if<!(std::is_base_of<SlidingFrictionInteraction, DUMMY>::value), void>::type
    getSlidingSpring(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //No sliding
    template<class DUMMY = FrictionForceInteraction>
    typename std::enable_if<!(std::is_base_of<SlidingFrictionInteraction, DUMMY>::value), void>::type
    setSlidingSpring(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    
    //Friction
    template<class DUMMY= FrictionForceInteraction>
    typename std::enable_if<std::is_base_of<FrictionInteraction, DUMMY>::value, void>::type
    getFrictionSprings(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        rollingSpring = interaction->FrictionForceInteraction::getRollingSpring();
        torsionSpring = interaction->FrictionForceInteraction::getTorsionSpring();
    }
    
    //Friction
    template<class DUMMY= FrictionForceInteraction>
    typename std::enable_if<std::is_base_of<FrictionInteraction, DUMMY>::value, void>::type
    setFrictionSprings(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        interaction->FrictionForceInteraction::setRollingSpring(rollingSpring);
        interaction->FrictionForceInteraction::setTorsionSpring(torsionSpring);
    }
    
    //No friction
    template<class DUMMY = FrictionForceInteraction>
    typename std::enable_if<!(std::is_base_of<FrictionInteraction, DUMMY>::value), void>::type
    getFrictionSprings(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //No friction
    template<class DUMMY = FrictionForceInteraction>
    typename std::enable_if<!(std::is_base_of<FrictionInteraction, DUMMY>::value), void>::type
    setFrictionSprings(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //Contact
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<(
            std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value
            || std::is_base_of<LiquidBridgeWilletInteraction, DUMMY>::value
            || std::is_base_of<IrreversibleAdhesiveInteraction, DUMMY>::value), void>::type
    getWasInContact(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        wasInContact = interaction->AdhesiveForceInteraction::getWasInContact();
    }
    
    //Contact
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<(
            std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value
            || std::is_base_of<LiquidBridgeWilletInteraction, DUMMY>::value
            || std::is_base_of<IrreversibleAdhesiveInteraction, DUMMY>::value), void>::type
    setWasInContact(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        interaction->AdhesiveForceInteraction::setWasInContact(wasInContact);
    }
    
    //No contact
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<!(
            std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value
            || std::is_base_of<LiquidBridgeWilletInteraction, DUMMY>::value
            || std::is_base_of<IrreversibleAdhesiveInteraction, DUMMY>::value), void>::type
    getWasInContact(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //No contact
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<!(
            std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value
            || std::is_base_of<LiquidBridgeWilletInteraction, DUMMY>::value
            || std::is_base_of<IrreversibleAdhesiveInteraction, DUMMY>::value), void>::type
    setWasInContact(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //Bonded
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<(std::is_base_of<BondedInteraction, DUMMY>::value), void>::type
    getBonded(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        bonded = interaction->AdhesiveForceInteraction::getBonded();
    }
    
    //Bonded
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<(std::is_base_of<BondedInteraction, DUMMY>::value), void>::type
    setBonded(Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        interaction->AdhesiveForceInteraction::setBonded(bonded);
    }
    
    //No bonded
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<!(std::is_base_of<BondedInteraction, DUMMY>::value), void>::type
    getBonded(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //No bonded
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<!(std::is_base_of<BondedInteraction, DUMMY>::value), void>::type
    setBonded(Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //Liquidbridge
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<(std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value), void>::type
    getLiquidBridge(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        liquidbridgeVolume = interaction->AdhesiveForceInteraction::getLiquidBridgeVolume();
    }
    
    //Liquidbridge
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<(std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value), void>::type
    setLiquidBridge(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        interaction->AdhesiveForceInteraction::setLiquidBridgeVolume(liquidbridgeVolume);
    }
    
    //No Liquidbridge
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<!(std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value), void>::type
    getLiquidBridge(
            const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
    //No Liquidbridge
    template<class DUMMY = AdhesiveForceInteraction>
    typename std::enable_if<!(std::is_base_of<LiquidMigrationWilletInteraction, DUMMY>::value), void>::type
    setLiquidBridge(
            Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }

    //Max overlap
    template<class DUMMY = NormalForceInteraction>
    typename std::enable_if<(std::is_base_of<LinearPlasticViscoelasticInteraction, DUMMY>::value), void>::type
    getMaximumOverlap(
                const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        maxOverlap = interaction->NormalForceInteraction::getMaxOverlap();
    }

    //Max overlap
    template<class DUMMY = NormalForceInteraction>
    typename std::enable_if<(std::is_base_of<LinearPlasticViscoelasticInteraction, DUMMY>::value), void>::type
    setMaximumOverlap(
                Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
        interaction->NormalForceInteraction::setMaxOverlap(maxOverlap);
    }

    //Max overlap
    template<class DUMMY = NormalForceInteraction>
    typename std::enable_if<!(std::is_base_of<LinearPlasticViscoelasticInteraction, DUMMY>::value), void>::type
    getMaximumOverlap(
                const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }

    //Max overlap
    template<class DUMMY = NormalForceInteraction>
    typename std::enable_if<!(std::is_base_of<LinearPlasticViscoelasticInteraction, DUMMY>::value), void>::type
    setMaximumOverlap(
                Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
    {
    }
    
};

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::MPIInteraction()
{
}

template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::copyFromInteraction(
        const Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction)
{
    //logger.assert(interaction->getP()->getId(),"Trying to copy an unreal interaction: P is not defined");
    //logger.assert(interaction->getP()->getId(),"Trying to copy an unreal interaction: I is not defined");
    
    if (interaction->getP() == nullptr) logger(WARN, "P is not defined!!");
    if (interaction->getI() == nullptr) logger(WARN, "I is not defined!!");
    P = interaction->getP()->getId();
    I = interaction->getI()->getId();
    
    if (dynamic_cast<const BaseParticle*>(interaction->getI()) == nullptr)
    {
        isWallInteraction = true;
    }
    else
    {
        isWallInteraction = false;
    }
    
    timeStamp = interaction->getTimeStamp();
    force = interaction->getForce();
    torque = interaction->getTorque();
    getSlidingSpring(interaction);
    getFrictionSprings(interaction);
    getWasInContact(interaction);
    getBonded(interaction);
    getLiquidBridge(interaction);
    getMaximumOverlap(interaction);
}


template<class NormalForceInteraction, class FrictionForceInteraction, class AdhesiveForceInteraction>
void MPIInteraction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>::copyToInteraction(
        Interaction<NormalForceInteraction, FrictionForceInteraction, AdhesiveForceInteraction>* interaction,
        const bool resetPointers)
{
    //Basic interaction values
    BaseInteraction* basicInteraction = static_cast<BaseInteraction*>(interaction);
    basicInteraction->setBasicMPIInteractionValues(P, I, timeStamp, force, torque, isWallInteraction, resetPointers);
    
    //Specific history interaction values, interaction type denpendent
    setSlidingSpring(interaction);
    setFrictionSprings(interaction);
    setWasInContact(interaction);
    setBonded(interaction);
    setLiquidBridge(interaction);
    setMaximumOverlap(interaction);
}


#endif
