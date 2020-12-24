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


#include "ChargedBondedInteraction.h"
#include "Species/AdhesiveForceSpecies/ChargedBondedSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
//#include <cassert>

///\todo Clean up this file by using the logger instead of cout, //cout, cerr and assert, and by motivating why the commented out code needs to be here.
///\todo Complete the documentation of these methods

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
ChargedBondedInteraction::ChargedBondedInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    //ensuring that, by default, particles are not 'bonded'
    //i.e. they will not unintentionally 'stick' to any overlapping particles!
    bonded_ = false;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ChargedBondedInteraction::ChargedBondedInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
ChargedBondedInteraction::ChargedBondedInteraction(const ChargedBondedInteraction& p)
        : BaseInteraction(p)
{
    //carrying the history parameter over for copied particles to ensure that any bonded particles
    //remain bonded!
    bonded_ = p.bonded_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ChargedBondedInteraction::ChargedBondedInteraction(const ChargedBondedInteraction &p finished"<<std::endl;
#endif
}

ChargedBondedInteraction::ChargedBondedInteraction()
{
#ifdef MERCURY_USE_MPI
    logger(FATAL,"ChargedBondedInteractions are currently not implemented in parallel MercuryDPM");
#endif
}

/*!
 *
 */
ChargedBondedInteraction::~ChargedBondedInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"ChargedBondedInteraction::ChargedBondedInteractionaction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] os
 */
void ChargedBondedInteraction::write(std::ostream& os  UNUSED) const
{
    os << " bonded " << bonded_;
}

/*!
 * \param[in] is
 */
void ChargedBondedInteraction::read(std::istream& is  UNUSED)
{
    std::string dummy;
    //logger(INFO,"ChargedBondedSpecies %",dummy);
    is >> dummy >> bonded_;
}

/*!
 *
 */
void ChargedBondedInteraction::computeAdhesionForce()
{
    
    const ChargedBondedSpecies* species = getSpecies();
    //std::cout << getSpecies()->getCharge() << std::endl;
    
    //creating local parameters to store the charges of both particles
    //involved in the interaction to allow for quick calculation
    const auto pSpecies = dynamic_cast<const ChargedBondedSpecies*>(getP()->getSpecies());
    const auto iSpecies = dynamic_cast<const ChargedBondedSpecies*>(getI()->getSpecies());
    logger.assert(pSpecies,"No ChargedBondedSpecies");
    logger.assert(iSpecies,"No ChargedBondedSpecies");
    const int pCharge = pSpecies->getCharge();
    const int iCharge = iSpecies->getCharge();
    
    //similarly, creating local parameters to store the relevant stiffness
    //and max force values
    const Mdouble k = species->getAdhesionStiffness();
    const Mdouble fMax = species->getAdhesionForceMax();
    
    const Mdouble kWaals = species->getVanDerWaalsStiffness();
    const Mdouble fMaxWaals = species->getVanDerWaalsForceMax();
    const Mdouble rWaals = fMaxWaals / kWaals;
    
    
    //First, adding bonded force if applicable
    if (bonded_ && getOverlap() >= 0)
    {
        addForce(getNormal() * (-species->getBondForceMax()
                                - species->getBondDissipation() * getNormalRelativeVelocity()));
        return;
    }
    
    
    //determining which of the three possible cases for force based on charge -
    //repulsive (like charges), attractive (unlike charges) or none (1 or more uncharged) -
    //is relevant for the current combination of particle charges...
    //(Note that the charge set function contains a safety check that means charge can only be
    // +/- 1, i.e. the expressions used below should not produce errors!
    
    //case 1 - 1 or more particles has no charge, i.e. no EM force between them
    if ((pCharge == 0) or (iCharge == 0))
    {
        //No need to write anything here as nothing needs to be returned!
        //std::cout << "no charge" << std::endl;
    }
        //case 2: unlike charges --> attractive force
    else if (pCharge == -iCharge)
    {
        //std::cout << "dissimilar charge" << std::endl;
        //in the case of particles with opposing charges, simply implements the
        //standard attractive force used for adhesive particle interactions
        if (getOverlap() >= 0)
        {
            addForce(getNormal() * -fMax);
            addForce(getNormal() * -fMaxWaals);
        }
        else if (getOverlap() >= -rWaals)
        {
            addForce(getNormal() * (-kWaals * getOverlap() - fMaxWaals));
            addForce(getNormal() * (-k * getOverlap() - fMax));
        }
        else
        {
            addForce(getNormal() * (-k * getOverlap() - fMax));
        }
    }
        //case 3: like charges --> repulsive force
    else if (pCharge == iCharge)
    {
        //std::cout << "similar charge" << std::endl;
        //in the case of particles with like charges, reverse the direction of the force applied
        //such that particles repel one another
        if (getOverlap() >= 0)
        {
            addForce(getNormal() * +fMax);
            addForce(getNormal() * -fMaxWaals);
            ///TW: I added the vdW force; KW, why was that force not active? Note, this change also shows up in energy
        }
        else if (getOverlap() >= -rWaals)
        {
            addForce(getNormal() * (-kWaals * getOverlap() - fMaxWaals));
            addForce(getNormal() * (+k * getOverlap() + fMax));
            //std::cout << "Waals = " << getNormal() * (-kWaals * getOverlap() - fMaxWaals) << std::endl;
        }
        else
        {
            addForce(getNormal() * (+k * getOverlap() + fMax));
        }
    }
        //if none of the above are satisfied, something must have gone very wrong!
    else
    {
        std::cerr << "Error: Particle charge has erroneous value" << std::endl;
        exit(-1);
    }
}

/*!
 * \details Elastic (=potential) energy is defined as the energy gained by separating two interactables.
 * As it costs energy to separate adhesive interactables, the elastic energy is negative.
 * \return the elastic energy stored in the adhesive spring.
 */
Mdouble ChargedBondedInteraction::getElasticEnergy() const
{
    const ChargedBondedSpecies* species = getSpecies();
    const auto pSpecies = static_cast<const ChargedBondedSpecies*>(getP()->getSpecies()->getAdhesiveForce());
    const auto iSpecies = static_cast<const ChargedBondedSpecies*>(getI()->getSpecies()->getAdhesiveForce());
    logger.assert(pSpecies,"No ChargedBondedSpecies");
    logger.assert(iSpecies,"No ChargedBondedSpecies");
    const int pCharge = pSpecies->getCharge();
    const int iCharge = iSpecies->getCharge();
    
    const Mdouble k = species->getAdhesionStiffness();
    const Mdouble fMax = species->getAdhesionForceMax();
    const Mdouble r = getBaseSpecies()->getInteractionDistance();
    
    const Mdouble kWaals = species->getVanDerWaalsStiffness();
    const Mdouble fMaxWaals = species->getVanDerWaalsForceMax();
    const Mdouble rWaals = (fMaxWaals == 0) ? 0 : (fMaxWaals / kWaals);
    
    
    //First, adding bonded force if applicable
    if (bonded_ && getOverlap() >= 0)
    {
        //comment to ignore BondForce
        Mdouble elasticEnergyAtEquilibrium = getElasticEnergyAtEquilibrium(species->getBondForceMax());
        return -species->getBondForceMax() * getOverlap() + elasticEnergyAtEquilibrium;
    }
    
    Mdouble elasticEnergy = 0.0;
    if ((pCharge != 0) && (iCharge != 0))
    {
        if (pCharge == -iCharge)
        {
            if (getOverlap() >= 0)
            {
                elasticEnergy -= (0.5 * rWaals + getOverlap()) * fMaxWaals
                                 + (0.5 * r + getOverlap()) * fMax;
            }
            else if (getOverlap() >= -rWaals)
            {
                elasticEnergy -= (0.5 * kWaals * mathsFunc::square(getOverlap() + rWaals)) +
                                 (0.5 * k * mathsFunc::square(getOverlap() + r));
            }
            else
            {
                elasticEnergy -= (0.5 * k * mathsFunc::square(getOverlap() + r));
            }
        }
            //case 3: like charges --> repulsive force
        else if (pCharge == iCharge)
        {
            if (getOverlap() >= 0)
            {
                elasticEnergy -= (0.5 * rWaals + getOverlap()) * fMaxWaals
                                 - (0.5 * r + getOverlap()) * fMax;
            }
            else if (getOverlap() >= -rWaals)
            {
                elasticEnergy -= (0.5 * kWaals * mathsFunc::square(getOverlap() + rWaals)) -
                                 (0.5 * k * mathsFunc::square(getOverlap() + r));
            }
            else
            {
                elasticEnergy += (0.5 * k * mathsFunc::square(getOverlap() + r));
            }
        }
        else
        {
            logger(ERROR, "Particle charge has erroneous value");
        }
    }
    return elasticEnergy;
}

/*!
 * \return a constant pointer to an instance of this class.
 */
const ChargedBondedSpecies* ChargedBondedInteraction::getSpecies() const
{
    return static_cast<const ChargedBondedSpecies*> (getBaseSpecies()->getAdhesiveForce()); //downcast
}

/*!
 * \return std::string
 */
std::string ChargedBondedInteraction::getBaseName() const
{
    return "ChargedBonded";
}

/*!
 * \details Used to set the private variable 'bonded' to true, thus allowing the user to choose to fix a given pair
 * of interacting, overlapping particles together
 */
void ChargedBondedInteraction::bond()
{
    bonded_ = true;
}

/*!
 * \details Used to set the private variable 'bonded' to false, thus allowing the user to choose to separate (or 'unbond')
 * a given pair of interacting, overlapping particles together which were previously fixed (bonded) together
 * Useful, for example, in implementing breakage mechanisms.
 */
void ChargedBondedInteraction::unbond()
{
    bonded_ = false;
}
