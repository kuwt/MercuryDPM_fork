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

#include "RegimeSinterSpecies.h"
#include "Logger.h"

RegimeSinterSpecies::RegimeSinterSpecies()
{
//    adhesionForceMax_ = 0;
//    adhesionStiffness_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"RegimeSinterSpecies::RegimeSinterSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] s the species that is copied
 */
RegimeSinterSpecies::RegimeSinterSpecies(const RegimeSinterSpecies& s)
{

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"RegimeSinterSpecies::RegimeSinterSpecies(const RegimeSinterSpecies &p) finished"<<std::endl;
#endif
}

RegimeSinterSpecies::~RegimeSinterSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"RegimeSinterSpecies::~RegimeSinterSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void RegimeSinterSpecies::write(std::ostream& os) const
{
//    os << " adhesionForceMax " << adhesionForceMax_;
//    os << " adhesionStiffness " << adhesionStiffness_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void RegimeSinterSpecies::read(std::istream& is)
{
//    std::string dummy;
//    is >> dummy >> adhesionForceMax_;
//    is >> dummy >> adhesionStiffness_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string RegimeSinterSpecies::getBaseName() const
{
    return "RegimeSinter";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void RegimeSinterSpecies::mix(RegimeSinterSpecies* const S, RegimeSinterSpecies* const T)
{
//    adhesionForceMax_ = BaseSpecies::average(S->getAdhesionForceMax(), T->getAdhesionForceMax());
//    adhesionStiffness_ = BaseSpecies::average(S->getAdhesionStiffness(), T->getAdhesionStiffness());
}