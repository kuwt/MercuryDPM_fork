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

#include "LinearSinterSpecies.h"

LinearSinterSpecies::LinearSinterSpecies()
    : LinearPlasticViscoelasticNormalSpecies()
{
    neckGrowthRate_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterSpecies::LinearSinterSpecies() finished"<<std::endl;
#endif
}

LinearSinterSpecies::LinearSinterSpecies(const LinearSinterSpecies &s)
    : LinearPlasticViscoelasticNormalSpecies(s)
{
    neckGrowthRate_=s.neckGrowthRate_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterSpecies::LinearSinterSpecies(const LinearSinterSpecies &p) finished"<<std::endl;
#endif
}

LinearSinterSpecies::~LinearSinterSpecies()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearSinterSpecies::~LinearSinterSpecies() finished"<<std::endl;
#endif   
}

//the name is set such that the full name does not extend
std::string LinearSinterSpecies::getBaseName() const
{
    return "Sinter";
}

Mdouble LinearSinterSpecies::getNeckGrowthRate() const
{
    return neckGrowthRate_;
}

void LinearSinterSpecies::setNeckGrowthRate(Mdouble neckGrowthRate)
{
    neckGrowthRate_=neckGrowthRate;
}

