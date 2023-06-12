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


#include "LiquidMigrationLSInteraction.h"
#include "Species/AdhesiveForceSpecies/LiquidMigrationLSSpecies.h"
#include "Particles/LiquidFilmParticle.h"
#include "InteractionHandler.h"
#include "DPMBase.h"
#include "cmath"
#include <iomanip>
#include <fstream>
#include <Interactions/FrictionForceInteractions/SlidingFrictionInteraction.h>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
LiquidMigrationLSInteraction::LiquidMigrationLSInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                   unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    liquidBridgeVolume_ = 0.0;
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    logger(INFO, "LiquidMigrationLSInteraction::LiquidMigrationLSInteraction() finished");
#endif
}

//used for mpi
LiquidMigrationLSInteraction::LiquidMigrationLSInteraction()
        : BaseInteraction()
{
    liquidBridgeVolume_ = 0.0;
    wasInContact_ = false;
#ifdef DEBUG_CONSTRUCTOR
    logger(INFO, "LiquidMigrationLSInteraction::LiquidMigrationLSInteraction() finished");
#endif
}

/*!
 * \param[in] p
 */
LiquidMigrationLSInteraction::LiquidMigrationLSInteraction(const LiquidMigrationLSInteraction& p)
        : BaseInteraction(p)
{
    liquidBridgeVolume_ = p.liquidBridgeVolume_;
    wasInContact_ = p.wasInContact_;
#ifdef DEBUG_CONSTRUCTOR
    logger(INFO, "LiquidMigrationLSInteraction::LiquidMigrationLSInteraction(const LiquidMigrationLSInteraction &p finished");
#endif
}

/*!
 * 
 */
LiquidMigrationLSInteraction::~LiquidMigrationLSInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    logger(INFO, "LiquidMigrationLSInteraction::~LiquidMigrationLSInteraction() finished");
#endif
}

void LiquidMigrationLSInteraction::actionsOnErase()
{
    rupture();
};

/*!
 * \param[in,out] os
 */
void LiquidMigrationLSInteraction::write(std::ostream& os) const
{
    os
            << " wasInContact " << wasInContact_
            << " liquidBridgeVolume " << liquidBridgeVolume_;
}

/*!
 * \param[in,out] is
 */
void LiquidMigrationLSInteraction::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> wasInContact_;
    is >> dummy >> liquidBridgeVolume_;
}

/*!
 * The contact model is based on the description given in (Zhang and Wu, 2020):
 * https://www.sciencedirect.com/science/article/pii/S0032591019311258
 * This model includes a normal capillary force, a normal lubrication force and a tangential lubrication force.
 * The capillary force is valid when particles are within rupture distance.
 * The formulation is in (Lian and Seville, 2015):
 * https://doi.org/10.1016/j.cis.2015.11.003.
 * The lubrication force is valid when particles are within limiting distance.
 * Note: the tangential lubrication force is added to sliding friction force through SlidingFrictionInteraction::addTangentialForce().
 */
void LiquidMigrationLSInteraction::computeAdhesionForce()
{
    // Adding no capillary force for liquid bridge volume = 0
    if (getLiquidBridgeVolume() == 0) return;
    
    if (getOverlap() >= 0)
    {
        // if particles are in contact
        const LiquidMigrationLSSpecies* species = getSpecies();
        const Mdouble effectiveRadius = getEffectiveRadius();
        const Mdouble fdotn = -2.0 * species->getSurfaceTension() * effectiveRadius *
                              (2.0 * constants::pi - 2.082 * pow(getHalfFillingAngle(), 0.9197)) *
                              (std::cos(species->getContactAngle()) - 0.16 * getHalfFillingAngle() *
                              (1.0 + getHalfFillingAngle()) * std::tan(species->getContactAngle())); //max value of capillary;
        addForce(getNormal() * fdotn);
    }
    else if (wasInContact_)
    {
        // if particles are not in contact, but within their interaction distance
        const LiquidMigrationLSSpecies* species = getSpecies();
        const Mdouble effectiveRadius = getEffectiveRadius();
        const Mdouble s_nondim = -getOverlap()/getRuptureDistance() * (1.0 + 0.5 * species->getContactAngle() *
                            std::cos(species->getContactAngle()) * std::cos(getHalfFillingAngle() + species->getContactAngle()) ); //dimensionless S
        const Mdouble fdotnc = -2.0 * species->getSurfaceTension() * effectiveRadius *
                                (2.0 * constants::pi - 2.082 * pow(getHalfFillingAngle(), 0.9197)) *
                                (std::cos(species->getContactAngle()) - 0.16 * getHalfFillingAngle() *
                                (1.0 + getHalfFillingAngle()) * std::tan(species->getContactAngle()))
                                /(1.0 + 0.8438 * pow(getHalfFillingAngle(), -0.724) * s_nondim +
                                1.3393 * pow(getHalfFillingAngle(), -1.459) * pow(s_nondim, 2.0));//capillary force

        Mdouble fdotnl = 0.0;//normal lubrication force in value
        Vec3D fdottl = Vec3D(0.0, 0.0, 0.0);//tangential lubrication force in vector

        if(-getOverlap() > getLimitingDistance())
        {
            fdotnl = -6.0 * constants::pi * species->getViscosity() * pow(effectiveRadius, 2.0)
                                    * getNormalRelativeVelocity()/(-getOverlap());

            Vec3D tangentialRelativeVelocity = getRelativeVelocity() - getNormal() * getNormalRelativeVelocity();

            fdottl = -6.0 * constants::pi * species->getViscosity() * effectiveRadius * tangentialRelativeVelocity
                     * (8.0/15.0 * std::log(effectiveRadius/(-getOverlap())) + 0.9588);

        }
        addForce(getNormal() * (fdotnc + fdotnl) + fdottl);

        //add tangential lubrication force to sliding friction force
        auto slidingFrictionInteraction = dynamic_cast<SlidingFrictionInteraction*>(this);
        slidingFrictionInteraction->addTangentialForce(fdottl);
    }
}

/// test if particle needs to be ruptured
void LiquidMigrationLSInteraction::actionsAfterTimeStep()
{
    if (wasInContact_)
    {
        if (-getOverlap() >= getRuptureDistance())
        {
            rupture();
        }
    }
    else
    {
        if (getOverlap() >= 0)
        {
            form();
        }
    }
}

void LiquidMigrationLSInteraction::form()
{
    //form a bridge
    //todo: extend to neighbours

    wasInContact_ = true;
    const LiquidMigrationLSSpecies* species = getSpecies();
    LiquidFilmParticle* IParticle = dynamic_cast<LiquidFilmParticle*>(getI());
    LiquidFilmParticle* PParticle = dynamic_cast<LiquidFilmParticle*>(getP());
    if (IParticle == nullptr) //if I is a wall
    {
        //do not form bridge if the volume is below minimum
        if (PParticle->getLiquidVolume() < species->getLiquidBridgeVolumeMin())
        {
            return;
        }
            //if below max bridge volume, move all liquid from film to volume
        else if (PParticle->getLiquidVolume() <= species->getLiquidBridgeVolumeMax())
        {
            liquidBridgeVolume_ = PParticle->getLiquidVolume();
            PParticle->setLiquidVolume(0.0);
        }
            //if above max bridge volume, fill the liquid bridge and keep the rest of the liquid in the particle
        else
        {
            liquidBridgeVolume_ = species->getLiquidBridgeVolumeMax();
            PParticle->setLiquidVolume(PParticle->getLiquidVolume() - species->getLiquidBridgeVolumeMax());
        }
//        if (liquidBridgeVolume_) logger(INFO,"Forming liquid bridge of volume % between particles % and wall %",liquidBridgeVolume_,getP()->getId(),getI()->getId());
    }
    else if (PParticle == nullptr) //if P is a wall
    {
        logger(ERROR,"Should not happen");
    }
    else //if P and I are particles
    {
        //if I is a ghost particle, apply volume change only to real particles (this removes the possibility that contacts are established only on one side of the periodic boundary; the same problem could exist for ruptures!)
        LiquidFilmParticle* IParticleReal;
        if (IParticle->getPeriodicFromParticle())
        {
            IParticleReal = dynamic_cast<LiquidFilmParticle*>(IParticle->getPeriodicFromParticle());
        }
        else
        {
            IParticleReal = IParticle;
        }

        Mdouble distributableLiquidVolume =
                PParticle->getLiquidVolume() + IParticleReal->getLiquidVolume();
        //assign all liquid of the contacting particles to the bridge,
        //if the total volume does not exceed LiquidBridgeVolumeMax
        ///\todo: maybe we need to check ghost particles?
        if (distributableLiquidVolume <= species->getLiquidBridgeVolumeMin())
        {
            return;
        }
        else if (distributableLiquidVolume <= species->getLiquidBridgeVolumeMax())
        {
            liquidBridgeVolume_ = distributableLiquidVolume;
            if (!IParticle->getPeriodicFromParticle() ||
                PParticle->getIndex() < IParticle->getPeriodicFromParticle()->getIndex())
            {
                PParticle->setLiquidVolume(0.0);
                IParticleReal->setLiquidVolume(0.0);
            }
        }
        else //if the total volume exceeds LiquidBridgeVolumeMax, only distribute the max value
        {
            liquidBridgeVolume_ = species->getLiquidBridgeVolumeMax();
            Mdouble pFraction =
                    PParticle->getLiquidVolume() / distributableLiquidVolume;
            if (!IParticle->getPeriodicFromParticle() ||
                PParticle->getIndex() < IParticle->getPeriodicFromParticle()->getIndex())
            {
                PParticle->addLiquidVolume(-pFraction * species->getLiquidBridgeVolumeMax());
                IParticleReal->addLiquidVolume(-(1.0 - pFraction) * species->getLiquidBridgeVolumeMax());
            }
        }
//        if (liquidBridgeVolume_) logger(INFO,"Forming liquid bridge of volume % between particles % and % (MPI %%, V % %, overlap %)",liquidBridgeVolume_,getP()->getId(),getI()->getId(),PParticle->isMPIParticle(),IParticle->isMPIParticle(),PParticle->getLiquidVolume(),IParticle->getLiquidVolume(),getOverlap());
    }
}

int LiquidMigrationLSInteraction::getNumberOfContacts(BaseInteractable* interactable)
{
    int numContacts = 0;
    for (auto i : interactable->getInteractions())
    {
        LiquidMigrationLSInteraction* j = dynamic_cast<LiquidMigrationLSInteraction*>(i);
        LiquidFilmParticle* jIParticle = dynamic_cast<LiquidFilmParticle*>(i->getI());
        if (j != this && jIParticle != nullptr && j->getLiquidBridgeVolume() != 0.0)
            numContacts++;
    }
    return numContacts;
}

void LiquidMigrationLSInteraction::rupture()
{
    // remove the contact history
    wasInContact_ = false;

    //if the bridge is already empty, do nothing
    if (getLiquidBridgeVolume() == 0.0)
        return;

    //else rupture a bridge
    const LiquidMigrationLSSpecies* species = getSpecies();
    LiquidFilmParticle* IParticle = dynamic_cast<LiquidFilmParticle*>(getI());
    LiquidFilmParticle* PParticle = dynamic_cast<LiquidFilmParticle*>(getP());
    if (IParticle == nullptr) //if I is a wall
    {
        int numContactsP = 0;
        for (auto i : getP()->getInteractions())
        {
            LiquidMigrationLSInteraction* j = dynamic_cast<LiquidMigrationLSInteraction*>(i);
            if (j != this && j != nullptr && j->getLiquidBridgeVolume() != 0.0)
            {
                numContactsP++;
            }

        }
//        logger(INFO,"Rupturing liquid bridge of volume % between particles % and wall % (numContacts %)",liquidBridgeVolume_,getP()->getId(),getI()->getId(),numContactsP);
        if (numContactsP > 0)
        {
            // Updating new volume with distribution less than critical volume
            Mdouble ExcessBridgeVolume = 0.0;
            Mdouble newVolume = liquidBridgeVolume_ * species->getDistributionCoefficient() / (numContactsP);

            for (auto i : getP()->getInteractions())
            {
                LiquidMigrationLSInteraction* j =
                        dynamic_cast<LiquidMigrationLSInteraction*>(i);
                if (j != this && j != nullptr && j->getLiquidBridgeVolume() != 0.0)
                {
                    j->setLiquidBridgeVolume(
                            j->getLiquidBridgeVolume() + newVolume);
                    if (j->getLiquidBridgeVolume() >=
                        species->getLiquidBridgeVolumeMax())
                    {
                        ExcessBridgeVolume +=
                                j->getLiquidBridgeVolume()
                                - species->getLiquidBridgeVolumeMax();
                        j->setLiquidBridgeVolume(
                                species->getLiquidBridgeVolumeMax());
                    }
                }
            }
            PParticle->setLiquidVolume(
                    PParticle->getLiquidVolume() + ExcessBridgeVolume +
                    liquidBridgeVolume_ * (1 - species->getDistributionCoefficient()));

            setLiquidBridgeVolume(0.0);
        }
        else
        {
            PParticle->setLiquidVolume(PParticle->getLiquidVolume() + liquidBridgeVolume_);
        }
        liquidBridgeVolume_ = 0.0;
        //for (auto i : *getHandler())
        //{
        //    LiquidMigrationLSInteraction* j =
        //           dynamic_cast<LiquidMigrationLSInteraction*>(i);
        //}
        //for (auto i : getHandler()->getDPMBase()->particleHandler)
        //{
        //    LiquidFilmParticle* j = dynamic_cast<LiquidFilmParticle*>(i);
        //}
    }
    else if (PParticle == nullptr) //if P is a wall
    {
        logger(ERROR,"this should not happen");
    }
    else //if P and I are particles
    {
        //count interaction partners of p (this contact, ghosts and interaction without liquid bridge dont count)
        int numContactsP = getNumberOfContacts(getP());
//        logger(INFO,"Rupturing liquid bridge of volume % between particles % and % (numContacts %, MPI %%, overlap %)",liquidBridgeVolume_,getP()->getId(),getI()->getId(),numContactsP,PParticle->isMPIParticle(),IParticle->isMPIParticle(),getOverlap());
        if (numContactsP < 1) //if P has only one contact (the one that gets ruptured), pass the fluid into it)
        {
            PParticle->addLiquidVolume(0.5 * liquidBridgeVolume_);
        }
        else //if P has only multiple contacts pass the fluid into it)
        {
            // move fluid to neighbouring contacts
            Mdouble perContactVolume = 0.5 * liquidBridgeVolume_ * species->getDistributionCoefficient() / numContactsP;
            for (auto i : getP()->getInteractions())
            {
                LiquidMigrationLSInteraction* j =
                        dynamic_cast<LiquidMigrationLSInteraction*>(i);
                LiquidFilmParticle* jIParticle =
                        dynamic_cast<LiquidFilmParticle*>(i->getI());
                if (j != this && jIParticle != nullptr && j->getLiquidBridgeVolume() != 0.0)
                {
                    Mdouble excessVolume =
                            perContactVolume - (species->getLiquidBridgeVolumeMax() - j->getLiquidBridgeVolume());
                    if (excessVolume < 0.0)
                    {
                        j->addLiquidBridgeVolume(perContactVolume);
                        PParticle->addLiquidVolume(
                                0.5 * liquidBridgeVolume_ * (1 - species->getDistributionCoefficient()) / numContactsP);
                    }
                    else
                    {
                        j->setLiquidBridgeVolume(species->getLiquidBridgeVolumeMax());
                        PParticle->addLiquidVolume(excessVolume);
                        PParticle->addLiquidVolume(
                                0.5 * liquidBridgeVolume_ * (1 - species->getDistributionCoefficient()) / numContactsP);
                    }
                }
            }
            //Mdouble PParticle_fin = PParticle->getLiquidVolume();

        }

        //count interaction partners of i (ghosts and interaction without liquid bridge dont count)
        int numContactsI = getNumberOfContacts(getI());
        if (numContactsI < 1)
        {
            IParticle->addLiquidVolume(0.5 * liquidBridgeVolume_);
            //std::cout << "added " << 0.5 * liquidBridgeVolume_
            //   << " to particle #" << IParticle->getIndex()
            //<< ", V=" << IParticle->getLiquidVolume() << std::endl;
        }
        else
        {
            // move fluid to neighbouring contacts
            Mdouble perContactVolume =
                    0.5 * liquidBridgeVolume_ * species->getDistributionCoefficient() / (numContactsI);
            for (BaseInteraction* i : getI()->getInteractions())
            {
                LiquidMigrationLSInteraction* j = dynamic_cast<LiquidMigrationLSInteraction*>(i);
                LiquidFilmParticle* jIParticle = dynamic_cast<LiquidFilmParticle*>(i->getI());
                if (j != this && jIParticle != nullptr && j->getLiquidBridgeVolume() != 0.0)
                {
                    Mdouble excessVolume =
                            perContactVolume
                            - species->getLiquidBridgeVolumeMax()
                            + j->getLiquidBridgeVolume();
                    if (excessVolume < 0.0)
                    {
                        j->addLiquidBridgeVolume(perContactVolume);
                        IParticle->addLiquidVolume(
                                0.5 * liquidBridgeVolume_ * (1 - species->getDistributionCoefficient()) /
                                (numContactsI));
                    }
                    else
                    {
                        //std::cout << "excess i-Volume " << excessVolume << std::endl;
                        j->setLiquidBridgeVolume(species->getLiquidBridgeVolumeMax());
                        IParticle->addLiquidVolume(excessVolume);
                        IParticle->addLiquidVolume(
                                0.5 * liquidBridgeVolume_ * (1 - species->getDistributionCoefficient()) /
                                (numContactsI));
                    }
                    //~ std::cout << "added " << perContactVolume 
                    //~ << " to i-contact #" << j->getIndex() 
                    //~ << ", V=" << j->getLiquidBridgeVolume() << std::endl;
                }
            }
        }

        //~ std::cout << "ruptured liquid bridge #" << getId()
        //~ << " between particles #"
        //~ << PParticle->getIndex()
        //~ << " and #"
        //~ << IParticle->getIndex()
        //~ << ": V" << liquidBridgeVolume_
        //~ << " -> Vp" << PParticle->getLiquidVolume()
        //~ << " Vi" << IParticle->getLiquidVolume()
        //~ << " Np" << numContactsP
        //~ << " Ni" << numContactsI
        //~ << std::endl;

        //to balance the added volume, remove the liquid from the bridge
        liquidBridgeVolume_ = 0.0;

    }
}

/*!
 * \return Mdouble
 */
Mdouble LiquidMigrationLSInteraction::getElasticEnergy() const
{
    ///\todo TW
    return 0.0;
}

/*!
 * \return const LiquidMigrationLSSpecies*
 */
const LiquidMigrationLSSpecies* LiquidMigrationLSInteraction::getSpecies() const
{
    return static_cast<const LiquidMigrationLSSpecies*>(getBaseSpecies()->getAdhesiveForce());

}

/*!
 * \return std::string
 */
std::string LiquidMigrationLSInteraction::getBaseName() const
{
    return "LiquidMigrationLS";
}

Mdouble LiquidMigrationLSInteraction::getLiquidBridgeVolume() const
{
    return liquidBridgeVolume_;
}

void LiquidMigrationLSInteraction::setLiquidBridgeVolume(Mdouble liquidBridgeVolume)
{
    liquidBridgeVolume_ = liquidBridgeVolume;
}

void LiquidMigrationLSInteraction::addLiquidBridgeVolume(Mdouble liquidBridgeVolume)
{
    liquidBridgeVolume_ += liquidBridgeVolume;
}

bool LiquidMigrationLSInteraction::getWasInContact() const
{
    return wasInContact_;
}

void LiquidMigrationLSInteraction::setWasInContact(bool wasInContact)
{
    wasInContact_ = wasInContact;
}

Mdouble LiquidMigrationLSInteraction::getRuptureDistance()
{
    const LiquidMigrationLSSpecies* species = getSpecies();
    return (1.0 + 0.5 * species->getContactAngle()) * (cbrt(liquidBridgeVolume_) + 0.1 * pow(liquidBridgeVolume_, 2.0/3.0));
}

Mdouble LiquidMigrationLSInteraction::getHalfFillingAngle()
{
    const LiquidMigrationLSSpecies* species = getSpecies();
    return pow(liquidBridgeVolume_, 0.26) + (0.285 - 0.4594 * species->getContactAngle()) * pow(liquidBridgeVolume_, (0.6327 + 0.897 * species->getContactAngle()));
}

Mdouble LiquidMigrationLSInteraction::getLimitingDistance()
{
    return getRuptureDistance()/4.0;
}

unsigned LiquidMigrationLSInteraction::getNumberOfFieldsVTK() const
{
    return 1;
}

std::string LiquidMigrationLSInteraction::getTypeVTK(unsigned i) const
{
    return "Float32";
}

std::string LiquidMigrationLSInteraction::getNameVTK(unsigned i) const
{
    return "liquidBridgeRadius";
    
}

std::vector<Mdouble> LiquidMigrationLSInteraction::getFieldVTK(unsigned i) const
{
    return std::vector<Mdouble>(1, cbrt(liquidBridgeVolume_));
}

Mdouble LiquidMigrationLSInteraction::getTotalLiquidFilmVolume(ParticleHandler& h)
{
    Mdouble volume = 0;
    for (auto p : h)
    {
        auto l = dynamic_cast<LiquidFilmParticle*>(p);
        if (l)
        {
            volume += l->getLiquidVolume();
        }
    }
    return volume;
}

Mdouble LiquidMigrationLSInteraction::getTotalLiquidBridgeVolume(InteractionHandler& h)
{
    Mdouble volume = 0;
    for (auto i : h)
    {
        auto l = dynamic_cast<LiquidMigrationLSInteraction*>(i);
        if (l)
        {
            volume += l->getLiquidBridgeVolume();
        }
    }
    return volume;
}
