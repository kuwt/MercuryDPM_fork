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

#ifndef SolidifyingLiquidMigrationWilletSpecies_h
#define SolidifyingLiquidMigrationWilletSpecies_h

#include "Species/AdhesiveForceSpecies/LiquidMigrationWilletSpecies.h"
#include "Math/ExtendedMath.h"
#include "Interactions/AdhesiveForceInteractions/SolidifyingLiquidMigrationWilletInteraction.h"
#include "../BaseSpecies.h"

/*!
 * \brief SolidifyingLiquidMigrationWilletSpecies contains the parameters used to describe a linear irreversible short-range force.
 * \details See SolidifyingLiquidMigrationWilletInteraction::computeForce for a description of the force law.
 */
class SolidifyingLiquidMigrationWilletSpecies : public LiquidMigrationWilletSpecies
{
public:
    ///\brief The correct Interaction type for this AdhesiveForceSpecies
    typedef SolidifyingLiquidMigrationWilletInteraction InteractionType;

    ///\brief The default constructor.
    SolidifyingLiquidMigrationWilletSpecies();

    ///\brief The default copy constructor.
    SolidifyingLiquidMigrationWilletSpecies(const SolidifyingLiquidMigrationWilletSpecies& s);

    ///\brief The default destructor.
    ~SolidifyingLiquidMigrationWilletSpecies();

    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);

    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;

    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;

    ///\brief creates default values for mixed species
    void mix(SolidifyingLiquidMigrationWilletSpecies* S, SolidifyingLiquidMigrationWilletSpecies* T);

//setters and getters
    // Allows the spring constant to be changed
    void setBondDissipation(Mdouble bondDissipation) {
        logger.assert_debug(bondDissipation >= 0, "Error in setBondDissipation");
        bondDissipation_ = bondDissipation;
    }

    // Allows the spring constant to be accessed
    Mdouble getBondDissipation() const {
        return bondDissipation_;
    }

    // Allows bondForceMax_ to be changed
    void setBondForceMax(Mdouble bondForceMax) {
        logger.assert_debug(bondForceMax >= 0, "Error in setBondForceMax");
        bondForceMax_ = bondForceMax;
    }

    // Allows bondForceMax_ to be accessed
    Mdouble getBondForceMax() const {
        return bondForceMax_;
    }

    // Allows bondForceMax_ to be changed
    void setSolidFraction(Mdouble solidFraction) {
        logger.assert_debug(solidFraction >= 0, "Error in setBondForceMax");
        solidFraction_ = solidFraction;
    }

    // Allows bondForceMax_ to be accessed
    Mdouble getSolidFraction() const {
        return solidFraction_;
    }

private:
    ///\brief dissipation in bond
    Mdouble bondDissipation_;

    ///\brief adhesion force at zero overlap
    Mdouble bondForceMax_;

    ///Different from bondInteraction:
    Mdouble solidFraction_;
};

#endif
