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

#include "SolidifyingLiquidMigrationWilletSpecies.h"
#include "Logger.h"

/*!
 * Default constructor for gluable species. Sets default values for all relevant parameters.
 * Note: if the stiffness of particles is left as zero, no force will be felt during interaction with
 * other particles
 * \param[in] s the species that is copied
 */
SolidifyingLiquidMigrationWilletSpecies::SolidifyingLiquidMigrationWilletSpecies()
{
    bondForceMax_ = 0;
    bondDissipation_ = 0;
    solidFraction_ = 0;
}

/*!
 * \param[in] the species that is copied
 */
SolidifyingLiquidMigrationWilletSpecies::SolidifyingLiquidMigrationWilletSpecies(const SolidifyingLiquidMigrationWilletSpecies& s) : LiquidMigrationWilletSpecies(s)
{
    bondForceMax_ = s.bondForceMax_;
    bondDissipation_ = s.bondDissipation_;
    solidFraction_ = s.solidFraction_;
}

SolidifyingLiquidMigrationWilletSpecies::~SolidifyingLiquidMigrationWilletSpecies()
{
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string SolidifyingLiquidMigrationWilletSpecies::getBaseName() const
{
    return "Solidifying" + LiquidMigrationWilletSpecies::getBaseName();
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void SolidifyingLiquidMigrationWilletSpecies::write(std::ostream& os) const
{
    LiquidMigrationWilletSpecies::write(os);
    os << " bondForceMax " << bondForceMax_;
    os << " bondDissipation " << bondDissipation_;
    os << " solidFraction " << solidFraction_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void SolidifyingLiquidMigrationWilletSpecies::read(std::istream& is)
{
    LiquidMigrationWilletSpecies::read(is);
    std::string dummy;
    is >> dummy >> bondForceMax_;
    is >> dummy >> bondDissipation_;
    is >> dummy >> solidFraction_;
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void SolidifyingLiquidMigrationWilletSpecies::mix(SolidifyingLiquidMigrationWilletSpecies* const S, SolidifyingLiquidMigrationWilletSpecies* const T)
{
    LiquidMigrationWilletSpecies::mix(S,T);
    bondForceMax_ = BaseSpecies::average(S->getBondForceMax(), T->getBondForceMax());
    bondDissipation_ = BaseSpecies::average(S->getBondDissipation(), T->getBondDissipation());
}

