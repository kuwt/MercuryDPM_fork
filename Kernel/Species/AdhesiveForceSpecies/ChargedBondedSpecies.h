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

#ifndef ChargedBondedSpecies_H
#define ChargedBondedSpecies_H

#include "Species/AdhesiveForceSpecies/BaseAdhesiveForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/AdhesiveForceInteractions/ChargedBondedInteraction.h"

/*!
 * \brief ChargedBondedSpecies contains the parameters used to describe a linear reversible short-range force.
 * \details See ChargedBondedInteraction::computeForce for a description of the force law.
 */
class ChargedBondedSpecies : public BaseAdhesiveForce
{
public:
    ///\brief The correct Interaction type for this AdhesiveForceSpecies
    typedef ChargedBondedInteraction InteractionType;
    
    ///\brief The default constructor.
    ChargedBondedSpecies();
    
    ///\brief The default constructor.
    ChargedBondedSpecies(const ChargedBondedSpecies& s);
    
    ///\brief The default constructor.
    ~ChargedBondedSpecies() = default;
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;
    
    ///\brief creates default values for mixed species
    void mix(ChargedBondedSpecies* S, ChargedBondedSpecies* T);

//adhesion-specific functions
    
    ///\brief returns the largest separation distance at which adhesive short-range forces can occur.
    void setInteractionDistance();

//setters and getters
    ///\brief Allows the spring constant to be changed
    void setAdhesionStiffness(Mdouble new_k0);
    
    ///\brief Allows the spring constant to be accessed
    Mdouble getAdhesionStiffness() const;
    
    ///\brief Allows the spring constant to be changed
    void setAdhesionForceMax(Mdouble new_f0);
    
    ///\brief Allows the spring constant to be accessed
    Mdouble getAdhesionForceMax() const;
    
    //function to allows the charge of a particle to be easily accessed
    int getCharge() const;
    
    //allows the user to set the charge possessed by a particle species
    void setCharge(int newCharge);

//*********************************************************************************************************************
//******************************ADDING ADDITIONAL FUNCTIONS FOR BOND INTERACTIONS**************************************
//*********************************************************************************************************************
    
    ///\brief Allows the spring constant for the BOND to be changed
    ///(Do not confuse with the charged interaction strength!)
    void setBondForceMax(Mdouble new_f0);
    
    ///\brief Allows the maximal force for 'bonding' particles together to be accessed
    Mdouble getBondForceMax() const;
    
    ///\brief Allows the additional dissipation used to damp oscillations between bondd particles to be changed
    void setBondDissipation(Mdouble disp);
    
    ///\brief Allows the additional dissipation used to damp oscillations between bondd particles to be accessed
    Mdouble getBondDissipation() const;

private:
    ///\brief stiffness of linear adhesion force
    Mdouble adhesionStiffness_;
    
    ///\brief adhesion force at zero overlap
    Mdouble adhesionForceMax_;
    
    //creating a simple integer flag to denote whether a particle is charged or not
    //the boolean will be 0 if the particle has no charge, 1 if positively charged
    //and -1 if negatively charged
    int charge_;

//*********************************************************************************************************************
//******************************ADDING ADDITIONAL VARIABLES FOR BOND INTERACTIONS**************************************
//*********************************************************************************************************************
    
    ///\brief The maximal force which acts to bind together particles which are "bondd" into a single body
    Mdouble bondForceMax_;
    
    ///\brief dissipation in bond
    ///\details the additional dissipation used to 'damp' oscillations between bondd particles
    Mdouble bondDissipation_;


//*********************************************************************************************************************
//****************************ADDING ADDITIONAL VARIABLES VAN DER WAALS INTERACTIONS***********************************
//*********************************************************************************************************************
    //Adding parameters to recreate a highly simplified (but relatively computationally efficient!)
    // van der Waals-like force at short distances
    
    //The maximal strength of the van der Waals force.
    //Note that this should, by definition, be great enough to overcome the maximal repulsive force experienced by
    //a particle, as the net force must be attractive in order to correctly represent vand der Waals
    Mdouble vanDerWaalsForceMax_;
    
    //The stiffness used to determine the 'reach' of the van der Waals force applied to particles.
    Mdouble vanDerWaalsStiffness_;

public:
    //Declaring the relevant get and set functions for the van der Waals stiffness and maximal force
    Mdouble getVanDerWaalsStiffness() const;
    
    Mdouble getVanDerWaalsForceMax() const;
    
    void setVanDerWaalsStiffness(Mdouble);
    
    void setVanDerWaalsForceMax(Mdouble);
    
    
};

#endif
