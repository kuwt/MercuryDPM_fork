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
#include "Math/ExtendedMath.h"
#include "Math/NumericalVector.h"
#include "Math/Vector.h"
#include "Multipole.h"
#include <cmath>
#include <iostream>
#include <vector>

Multipole::Multipole(int p, NumericalVector<>* squaredFactorials, Vec3D location) :
        p_(p),
        squaredFactorials_(squaredFactorials),
        location_(location)
{
}

Multipole::~Multipole()
= default;

void Multipole::computeMultipoleExpansion()
{
    int nTerms = 0.5 * (p_ + 1) * (2 * p_ + 2);
    multipoleExpansionCoefficients_(nTerms);
}

NumericalVector<std::complex<Mdouble>> Multipole::TranslateMultipoleExpansionTo(Vec3D location)
{
    //todo: Find a better name for A/squaredFactorials
    int nTerms = 0.5 * (p_ + 1) * (2 * p_ + 2);
    NumericalVector<std::complex<Mdouble>> translatedMultipoleCoefficients(nTerms);
    
    //Check if a multipole expansion has taken place
    if (multipoleExpansionCoefficients_.size() == 0)
    {
        std::cout << "Multipole is not yet expanded." << std::endl;
        std::exit(-1);
    }
    
    //Determine rho, alpha and beta
    Vec3D distance = location_ - location;
    
    //Todo: Add a quarternion step in here with distance as input, to avoid NaN values in the angles.
    Mdouble rho = 1.0;
    Mdouble alpha = 1.0;
    Mdouble beta = 1.0;

/*	std::cout << "rho=" << rho << std::endl;
	std::cout << "alpha=" << alpha << std::endl;
	std::cout << "beta=" << beta <<std::endl;*/
    
    //Calculate spherical harmonics for alpha and beta
    NumericalVector<std::complex<Mdouble>> sphericalHarmonics = sphericalHarmonics::sphericalHarmonics(p_, alpha, beta);
    
    //Compute the transfered multipole coefficients
    for (int j = 0; j <= p_; j++)
    {
        for (int k = -j; k <= j; k++)
        {
            std::complex<Mdouble> result = 0.0;
            int location = j * j + (k + j);
            for (int n = 0; n <= j; n++)
            {
                int a = std::max(k + n - j, -n);
                int b = std::min(k + j - n, n);
                for (int m = a; m <= b; m++)
                {
                    int location_O = (j - n) * (j - n) + ((k - m) + (j - n));
                    int location_A1 = n * n + (m + n);
                    int location_Y = n * n + (-m + n);
                    int location_A2 = location_O;
                    int location_A3 = location;
                    result += multipoleExpansionCoefficients_[location_O] *
                              std::pow(constants::i, (std::abs(k) - std::abs(m) - std::abs(k - m)))
                              * (*squaredFactorials_)(location_A1) * (*squaredFactorials_)(location_A2) *
                              std::pow(rho, n) * sphericalHarmonics[location_Y] / (*squaredFactorials_)(location_A3);
                }
            }
            translatedMultipoleCoefficients[location] = result;
        }
    }
    
    return translatedMultipoleCoefficients;
}

NumericalVector<std::complex<Mdouble>> Multipole::convertMultipoleToLocal(Vec3D location)
{
    int nTerms = 0.5 * (p_ + 1) * (2 * p_ + 2);
    NumericalVector<std::complex<Mdouble>> localExpansionCoefficients(nTerms);
    
    //Compute alpha and beta;
    //Todo: use quarternions to do this shite
    Mdouble rho = 1.0;
    Mdouble alpha = 1.0;
    Mdouble beta = 1.0;
    
    NumericalVector<std::complex<Mdouble>> sphericalHarmonics = sphericalHarmonics::sphericalHarmonics(2 * p_, alpha,
                                                                                                       beta);
    
    for (int j = 0; j <= p_; j++)
    {
        for (int k = -j; k <= j; k++)
        {
            std::complex<Mdouble> result = 0.0;
            int location = j * j + (k + j);
            for (int n = 0; n <= p_; n++)
            {
                for (int m = -n; m <= n; m++)
                {
                    int location_A1 = n * n + (m + n);
                    int location_A2 = j * j + (j + k);
                    int location_A3 = (j + n) * (j + n) + ((m - k) + (j + n));
                    int location_Y = location_A3;
                    int location_O = location_A1;
                    std::complex<Mdouble> J = std::pow(constants::i, std::abs(k - m) - std::abs(k) - std::abs(m));
                    //\todo TW note: a warning says += cannot be done here
                    result += multipoleExpansionCoefficients_[location_O] * J * (*squaredFactorials_)(location_A1) *
                              (*squaredFactorials_)(location_A2) * sphericalHarmonics[location_Y] /
                              ((*squaredFactorials_)(location_A3) * std::pow(rho, (j + n + 1)));
                    
                }
            }
            localExpansionCoefficients[location] = result;
        }
    }
    return localExpansionCoefficients;
}

/// Adds multipole coefficients to an existing multipole
/// todo: remove this function; it should not be required anymore
void Multipole::addMultipoleCoefficients(NumericalVector<std::complex<Mdouble>> multipoleExpansionCoefficients)
{
    if (multipoleExpansionCoefficients.size() > multipoleExpansionCoefficients_.size())
    {
        std::cout << "Multipole expansion coefficient sizes are not correct." << std::endl;
        std::exit(-1);
    }
    
    multipoleExpansionCoefficients_ += multipoleExpansionCoefficients;
}
