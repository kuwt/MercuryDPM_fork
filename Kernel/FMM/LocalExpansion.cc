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
#include "LocalExpansion.h"
#include "Math/ExtendedMath.h"
#include "Math/NumericalVector.h"
#include "Math/Vector.h"
#include "Multipole.h"
#include <complex>
#include <vector>

LocalExpansion::LocalExpansion(int p, NumericalVector<>* squaredFactorials, Vec3D location) :
        p_(p),
        squaredFactorials_(squaredFactorials),
        location_(location)
{
}

LocalExpansion::~LocalExpansion()
= default;

void LocalExpansion::initialiseLocalExpansion()
{
    size_t nTerms = 0.5 * (p_ + 1) * (2 * p_ + 2);
    localExpansionCoefficients_(nTerms);
}

NumericalVector<std::complex<Mdouble>> LocalExpansion::translateLocalExpansion(Vec3D location)
{
    int nTerms = 0.5 * (p_ + 1) * (2 * p_ + 2);
    NumericalVector<std::complex<Mdouble>> translatedLocalExpansionCoefficients(nTerms);
    //std::cout << "size: " << nTerms << std::endl;
    
    //compute angles and distance in new framework
    //todo: fix this rubble with quaternions
    Mdouble rho = 1.0;
    Mdouble alpha = 1.0;
    Mdouble beta = 1.0;
    
    
    //Compute translated local expansion coefficients
    NumericalVector<std::complex<Mdouble>> sphericalHarmonics = sphericalHarmonics::sphericalHarmonics(p_, alpha, beta);
    
    for (int j = 0; j <= p_; j++)
    {
        for (int k = -j; k <= j; k++)
        {
            std::complex<Mdouble> result = {0.0, 0.0};
            int location = j * j + (k + j);
            for (int n = j; n <= p_; n++)
            {
                for (int m = (k - j + n); m <= (k - n + j); m++)
                {
                    int location_O = n * n + (m + n);
                    int location_A1 = (n - j) * (n - j) + ((m - k) + (n - j));
                    int location_A2 = location;
                    int location_Y = location_A1;
                    int location_A3 = location_O;
                    std::complex<Mdouble> J = std::pow(constants::i, std::abs(m) - std::abs(m - k) - std::abs(k));
                    /*				std::cout << "location_A1: " << location_O << std::endl;
                                    std::cout << "locationExpansion: " << localExpansionCoefficients_[location_O] << std::endl;
                                    std::cout << "A1: " << (*squaredFactorials_)(location_A1) << std::endl;
                                    std::cout << "A2: " << (*squaredFactorials_)(location_A2) << std::endl;
                                    std::cout << "A3: " << (*squaredFactorials_)(location_A3) << std::endl;
                                    std::cout << "Spherical: " << sphericalHarmonics[location_Y] << std::endl;*/
                    //\todo TW note: a warning says += cannot be done here
                    result += localExpansionCoefficients_[location_O] * J * (*squaredFactorials_)(location_A1) *
                              (*squaredFactorials_)(location_A2) * sphericalHarmonics[location_Y] *
                              std::pow(rho, n - j) / (std::pow(-1, n + j) * (*squaredFactorials_)(location_A3));
                }
            }
            //std::cout << "location: " << location << std::endl;
            translatedLocalExpansionCoefficients[location] = result;
        }
    }
    return translatedLocalExpansionCoefficients;
}

void LocalExpansion::addLocalExpansionCoefficients(NumericalVector<std::complex<Mdouble>> localExpansionCoefficients)
{
    if (localExpansionCoefficients.size() > localExpansionCoefficients_.size())
    {
        std::cout << "Multipole expansion coefficient sizes are not correct." << std::endl;
        std::exit(-1);
    }
    
    localExpansionCoefficients_ += localExpansionCoefficients;
    
}
