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
#include "Dipole.h"
#include "Multipole.h"
#include "Math/ExtendedMath.h"
#include "Math/NumericalVector.h"
#include <complex>
#include <vector>

Dipole::Dipole(int p, NumericalVector<>* squaredFactorials, Vec3D location, Vec3D velocity, Mdouble strength) :
        Multipole(p, squaredFactorials, location),
        velocity_(velocity),
        strength_(strength)
{
}

void Dipole::computeMultipoleExpansion()
{
    size_t nTerms = (p_ + 1) * (p_ + 1);
    NumericalVector<std::complex<Mdouble>> multipoleExpansionCoefficients(nTerms);
    
    //Calculate dipole coefficients for spherical harmonics
    Mdouble s1 = strength_ * velocity_.getComponent(1);
    Mdouble s2 = strength_ * velocity_.getComponent(2);
    Mdouble s3 = strength_ * velocity_.getComponent(3);
    
    multipoleExpansionCoefficients[1] = 1.0 / std::sqrt(2.0) * (-s1 + s2 / constants::i);
    multipoleExpansionCoefficients[2] = s3;
    multipoleExpansionCoefficients[3] = -1.0 / (sqrt(2.0)) * (s1 + s2 / constants::i);
    
    multipoleExpansionCoefficients_ = multipoleExpansionCoefficients;
}

