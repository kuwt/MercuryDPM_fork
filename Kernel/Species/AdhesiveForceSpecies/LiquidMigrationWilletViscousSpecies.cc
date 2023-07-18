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

#include "LiquidMigrationWilletViscousSpecies.h"
#include "Logger.h"

LiquidMigrationWilletViscousSpecies::LiquidMigrationWilletViscousSpecies()
        : LiquidMigrationWilletSpecies()
{
    viscosity_ = 0.0;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidMigrationWilletViscousSpecies::LiquidMigrationWilletViscousSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
LiquidMigrationWilletViscousSpecies::LiquidMigrationWilletViscousSpecies(const LiquidMigrationWilletViscousSpecies& s)
        : LiquidMigrationWilletSpecies(s)
{
    viscosity_ = s.viscosity_;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidMigrationWilletViscousSpecies::LiquidMigrationWilletViscousSpecies(const LiquidMigrationWilletViscousSpecies &p) finished"<<std::endl;
#endif
}

LiquidMigrationWilletViscousSpecies::~LiquidMigrationWilletViscousSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LiquidMigrationWilletViscousSpecies::~LiquidMigrationWilletViscousSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void LiquidMigrationWilletViscousSpecies::write(std::ostream& os) const
{
    LiquidMigrationWilletSpecies::write(os);
    os << " viscosity_ " << viscosity_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void LiquidMigrationWilletViscousSpecies::read(std::istream& is)
{
    LiquidMigrationWilletSpecies::read(is);
    std::string dummy;
    is >> dummy >> viscosity_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string LiquidMigrationWilletViscousSpecies::getBaseName() const
{
    return "LiquidMigrationWilletViscous";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void LiquidMigrationWilletViscousSpecies::mix(LiquidMigrationWilletViscousSpecies* const S, LiquidMigrationWilletViscousSpecies* const T)
{
    LiquidMigrationWilletSpecies::mix(S, T);
    viscosity_ = BaseSpecies::average(S->getViscosity(), T->getViscosity());
}


/*!
 * \param[in] viscosity the viscosity of the liquid.
 */
void LiquidMigrationWilletViscousSpecies::setViscosity(Mdouble viscosity)
{
    if (viscosity >= 0)
    {
        viscosity_ = viscosity;
    }
    else
    {
        std::cerr << "Error in setViscosity" << std::endl;
        exit(-1);
    }
}

/*!
 * \return the the viscosity of the liquid.
 */
Mdouble LiquidMigrationWilletViscousSpecies::getViscosity() const
{
    return viscosity_;
}


