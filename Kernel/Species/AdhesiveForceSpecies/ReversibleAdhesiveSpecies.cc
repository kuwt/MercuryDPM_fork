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

#include "ReversibleAdhesiveSpecies.h"
#include "Logger.h"
#include "Species/BaseSpecies.h"

ReversibleAdhesiveSpecies::ReversibleAdhesiveSpecies()
{
    adhesionForceMax_ = 0;
    adhesionStiffness_ = 0;
    //interaction distance not yet determined
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ReversibleAdhesiveSpecies::ReversibleAdhesiveSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] s the species that is copied
 */
ReversibleAdhesiveSpecies::ReversibleAdhesiveSpecies(const ReversibleAdhesiveSpecies& s)
{
    adhesionForceMax_ = s.adhesionForceMax_;
    adhesionStiffness_ = s.adhesionStiffness_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"ReversibleAdhesiveSpecies::ReversibleAdhesiveSpecies(const ReversibleAdhesiveSpecies &p) finished"<<std::endl;
#endif
}

ReversibleAdhesiveSpecies::~ReversibleAdhesiveSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"ReversibleAdhesiveSpecies::~ReversibleAdhesiveSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void ReversibleAdhesiveSpecies::write(std::ostream& os) const
{
    os << " adhesionForceMax " << adhesionForceMax_;
    os << " adhesionStiffness " << adhesionStiffness_;
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void ReversibleAdhesiveSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> adhesionForceMax_;
    is >> dummy >> adhesionStiffness_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string ReversibleAdhesiveSpecies::getBaseName() const
{
    return "ReversibleAdhesive";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void ReversibleAdhesiveSpecies::mix(ReversibleAdhesiveSpecies* const S, ReversibleAdhesiveSpecies* const T)
{
    adhesionForceMax_ = BaseSpecies::average(S->getAdhesionForceMax(), T->getAdhesionForceMax());
    adhesionStiffness_ = BaseSpecies::average(S->getAdhesionStiffness(), T->getAdhesionStiffness());
}

///\return the maximum separation distance below which adhesive forces can occur (needed for contact detection)
void ReversibleAdhesiveSpecies::setInteractionDistance()
{
    logger.assert(adhesionStiffness_ != 0.0,"ReversibleAdhesiveSpecies::getInteractionDistance(): adhesionStiffness cannot be zero");
    getBaseSpecies()->setInteractionDistance(adhesionForceMax_ / adhesionStiffness_);
}

///Allows the spring constant to be changed
void ReversibleAdhesiveSpecies::setAdhesionStiffness(Mdouble adhesionStiffness)
{
    logger.assert(adhesionStiffness >= 0,"Error in setAdhesionStiffness");
    adhesionStiffness_ = adhesionStiffness;
    setInteractionDistance();
}

///Allows the spring constant to be accessed
Mdouble ReversibleAdhesiveSpecies::getAdhesionStiffness() const
{
    return adhesionStiffness_;
}

///Allows the spring constant to be changed
void ReversibleAdhesiveSpecies::setAdhesionForceMax(Mdouble adhesionForceMax)
{
    logger.assert(adhesionForceMax >= 0,"Error in setAdhesionForceMax");
    adhesionForceMax_ = adhesionForceMax;
    if (adhesionStiffness_!=0) setInteractionDistance();
}

///Allows the spring constant to be accessed
Mdouble ReversibleAdhesiveSpecies::getAdhesionForceMax() const
{
    return adhesionForceMax_;
}
