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

#include "ChargedBondedSpecies.h"
#include "Logger.h"
#include "Species/BaseSpecies.h"

/*!
 * Default constructor for charged species. Sets default values for all relevant parameters.
 * Note: if the stiffness of particles is left as zero, no force will be felt during interaction with
 * other particles
 * \param[in] s the species that is copied
 */
ChargedBondedSpecies::ChargedBondedSpecies()
{
    //setting adhesion values initially to zero such that, by default, particles do not experience
    //long range forces
    adhesionForceMax_ = 0;
    adhesionStiffness_ = 1;
    //similarly, setting bond properties to zero such that, by default, particles cannot be 'bondd'...
    bondForceMax_ = 0;
    //...and do not impart any excess dissipation!
    bondDissipation_ = 0;
    //Setting also parameters corresponding to the van der Waals force to zero by default
    vanDerWaalsForceMax_ = 0;
    vanDerWaalsStiffness_ = 1;
    charge_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ChargedBondedSpecies::ChargedBondedSpecies() finished"<<std::endl;
#endif
}

/*!
 * Copy constructor for charged species
 * \param[in] s the species that is copied
 */
ChargedBondedSpecies::ChargedBondedSpecies(const ChargedBondedSpecies& s)
{
    adhesionForceMax_ = s.adhesionForceMax_;
    adhesionStiffness_ = s.adhesionStiffness_;
    bondForceMax_ = s.bondForceMax_;
    bondDissipation_ = s.bondDissipation_;
    vanDerWaalsForceMax_ = s.vanDerWaalsForceMax_;
    vanDerWaalsStiffness_ = s.vanDerWaalsStiffness_;
    charge_ = s.charge_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ChargedBondedSpecies::ChargedBondedSpecies(const ChargedBondedSpecies &p) finished"<<std::endl;
#endif
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void ChargedBondedSpecies::write(std::ostream& os) const
{
    os << " adhesionForceMax " << adhesionForceMax_;
    os << " adhesionStiffness " << adhesionStiffness_;
    os << " charge " << charge_;
    os << " bondForceMax " << bondForceMax_;
    os << " bondDissipation " << bondDissipation_;
    os << " vanDerWaalsForceMax " << vanDerWaalsForceMax_;
    os << " vanDerWaalsStiffness " << vanDerWaalsStiffness_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void ChargedBondedSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> adhesionForceMax_;
    is >> dummy >> adhesionStiffness_;
    is >> dummy >> charge_;
    is >> dummy >> bondForceMax_;
    is >> dummy >> bondDissipation_;
    is >> dummy >> vanDerWaalsForceMax_;
    is >> dummy >> vanDerWaalsStiffness_;
}

/*!
 * Returns the name of thebcurrent species
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string ChargedBondedSpecies::getBaseName() const
{
    return "ChargedBonded";
}

/*!
 * \details Calling this function acts to assign values of maximum adhesion force and stiffness
 * for mixed species interactions.
 * For all parameters we assume that the harmonic mean of the parameters of the
 * original two species is a sensible default.
 *
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void ChargedBondedSpecies::mix(ChargedBondedSpecies* const S, ChargedBondedSpecies* const T)
{
    //'mixing' properties relating to the charged force interactions between particles
    setAdhesionForceMax(BaseSpecies::average(S->getAdhesionForceMax(), T->getAdhesionForceMax()));
    setAdhesionStiffness(BaseSpecies::average(S->getAdhesionStiffness(), T->getAdhesionStiffness()));
    
    //ensuring that, in addition, bond properties are also 'mixed'
    bondForceMax_ = BaseSpecies::average(S->getBondForceMax(), T->getBondForceMax());
    bondDissipation_ = BaseSpecies::average(S->getBondDissipation(), T->getBondDissipation());
    
    //and that the van der Waals force also is computed correctly
    vanDerWaalsForceMax_ = BaseSpecies::average(S->getVanDerWaalsForceMax(), T->getVanDerWaalsForceMax());
    vanDerWaalsStiffness_ = BaseSpecies::average(S->getVanDerWaalsStiffness(), T->getVanDerWaalsStiffness());
    
    charge_ = 0; //note mixedSpecies dont need charge, it's a particle property.
}

///\return the maximum separation distance below which adhesive forces can occur (needed for contact detection)
void ChargedBondedSpecies::setInteractionDistance()
{
    logger.assert(adhesionStiffness_ != 0.0,"ChargedBondedSpecies::getInteractionDistance(): adhesionStiffness cannot be zero");
    getBaseSpecies()->setInteractionDistance(adhesionForceMax_ / adhesionStiffness_);
}


///Allows the spring constant to be changed
void ChargedBondedSpecies::setAdhesionStiffness(Mdouble new_k0)
{
    if (new_k0 >= 0) {
        adhesionStiffness_ = new_k0;
        setInteractionDistance();
    } else {
        logger(ERROR, "Error in setAdhesionStiffness");
    }
}

///Allows the spring constant to be accessed
Mdouble ChargedBondedSpecies::getAdhesionStiffness() const
{
    return adhesionStiffness_;
}

///Allows the spring constant to be changed
void ChargedBondedSpecies::setAdhesionForceMax(Mdouble new_f0)
{
    if (new_f0 >= 0) {
        adhesionForceMax_ = new_f0;
        setInteractionDistance();
    } else {
        logger(ERROR, "Error in setAdhesionForceMax");
    }
}

///Allows the spring constant to be accessed
Mdouble ChargedBondedSpecies::getAdhesionForceMax() const
{
    return adhesionForceMax_;
}

//overwrites the baseSpecies version of this function, allowing particles' charges to be accessed from the
//program running
int ChargedBondedSpecies::getCharge() const
{
    return charge_;
}

//allows the user to manually set the charge of a particle
void ChargedBondedSpecies::setCharge(int newCharge)
{
    //making sure that the user can only enter 1 (positive charge) -1 (negative charge)
    //or zero (no charge)
    if (newCharge == 0 || newCharge == 1 || newCharge == -1)
        charge_ = newCharge;
    else
    {
        std::cerr << "Error in setCharge - charge must be +1, -1 or zero" << std::endl;
        exit(-1);
    }
}

//*********************************************************************************************************************
//******************************ADDING ADDITIONAL VARIABLES FOR BOND INTERACTIONS**************************************
//*********************************************************************************************************************

///\details Allows the spring constant to be changed
///This means that we can alter the strength with which particles are "bondd" together and, thus, the mean
///distance by which bondd particles overlap by. As such, we can alter the geometry of multi-particle structures formed.
void ChargedBondedSpecies::setBondForceMax(Mdouble new_f0)
{
    if (new_f0 >= 0)
    {
        bondForceMax_ = new_f0;
    }
    else
    {
        std::cerr << "Error in setBondForceMax" << std::endl;
        exit(-1);
    }
}

///Allows the the maximal force of the 'bond' used to bond particles together to be accessed
Mdouble ChargedBondedSpecies::getBondForceMax() const
{
    return bondForceMax_;
}

///Allows the value of the dissipation between bondd particles to be changed. This can be used to eliminate
///oscillations which may arise due to the forces holding the particles together competing with the usual
///forces which act to push particles in contact apart!
void ChargedBondedSpecies::setBondDissipation(Mdouble disp)
{
    if (disp >= 0)
        bondDissipation_ = disp;
    else
    {
        std::cerr << "Error in setAdhesionDissipation" << std::endl;
        exit(-1);
    }
}

///Allows the value of the dissipation between bondd particles to be accessed. This can be used to eliminate
///oscillations which may arise due to the forces holding the particles together competing with the usual
///forces which act to push particles in contact apart!
Mdouble ChargedBondedSpecies::getBondDissipation() const
{
    return bondDissipation_;
}


//*********************************************************************************************************************
//****************************ADDING ADDITIONAL VARIABLES VAN DER WAALS INTERACTIONS***********************************
//*********************************************************************************************************************
//Adding parameters to recreate a highly simplified (but relatively computationally efficient!)
// van der Waals-like force at short distances

//A function to set the maximum strength of the van der Waals force.
//Note that this should, by definition, be great enough to overcome the maximal repulsive force experienced by
//a particle, as the net force must be attractive in order to correctly represent vand der Waals
void ChargedBondedSpecies::setVanDerWaalsForceMax(Mdouble fMax)
{
    if (fMax < 0)
    {
        std::cerr << "Error in setVanDerWaalsForceMax." << std::endl;
        exit(-1);
    }
    else
    {
        vanDerWaalsForceMax_ = fMax;
    }
}

//a function to return the value of the van der Waals maximal force
Mdouble ChargedBondedSpecies::getVanDerWaalsForceMax() const
{
    return vanDerWaalsForceMax_;
}

//A function to set the stiffness and hence - when combined with the maximal van der Waals force - *range*
//of the van der Waals force applied to closely interacting particles.
void ChargedBondedSpecies::setVanDerWaalsStiffness(Mdouble stiffness)
{
    vanDerWaalsStiffness_ = stiffness;
}

//a function to return the value of the van der Waals maximal force
Mdouble ChargedBondedSpecies::getVanDerWaalsStiffness() const
{
    return vanDerWaalsStiffness_;
}
