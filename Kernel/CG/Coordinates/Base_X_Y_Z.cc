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

#include "Base_X_Y_Z.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

Mdouble Base_X_Y_Z::getGaussPrefactor(Mdouble width, Mdouble cutoff)
{
    Mdouble prefactor = 1.0 / (constants::sqrt_2 * constants::sqrt_pi * width);
    return prefactor / erf(cutoff / (constants::sqrt_2 * width));
}

Mdouble Base_X_Y_Z::getGaussIntegralPrefactor(Mdouble distance, Mdouble width, Mdouble cutoff)
{
    Mdouble widthSqrt2 = width * constants::sqrt_2;
    Mdouble a = -cutoff;
    Mdouble b = cutoff + distance;
    return 0.5 / (
            +erf(b / widthSqrt2) * b
            + widthSqrt2 / constants::sqrt_pi * exp(-mathsFunc::square(b / widthSqrt2))
            - erf(a / widthSqrt2) * a
            - widthSqrt2 / constants::sqrt_pi * exp(-mathsFunc::square(a / widthSqrt2))
    );
}

/*!
 * \details The volume is computed as
 * \f[volume=\int_0^1\sum_{i=1}^n c_i r^i 2 dr = 2 \sum_{i=1}^n c_i/(i+1) \f]
 */
void Base_X_Y_Z::normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff)
{
    Mdouble volume = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); i++)
        volume += coefficients[i] / static_cast<Mdouble>(i + 1);
    volume *= 2.0 * cutoff;
    for (double& coefficient : coefficients)
        coefficient /= volume;
    //logger(INFO,"Volume %",volume);
}

const unsigned Base_X_Y_Z::countVariables()
{
    return 1;
}

//template<StatType T>
//double NORMALIZED_POLYNOMIAL<T>::get_volume()
//{
//    double volume = 0;
//    unsigned int N = coefficients.size();
//    if (dim == 3)
//    {
//        for (unsigned int i = 0; i < N; i++)
//            volume += coefficients[i] / (2. + N - i);
//        volume *= 4. * constants::pi;
//    }
//    else if (dim == 2)
//    {
//        std::cerr << "dim=2 is not working yet" << std::endl;
//        exit(-1);
//        for (unsigned int i = 0; i < coefficients.size(); i++)
//            volume += coefficients[i] / (1. + N - i);
//        volume *= 2. * constants::pi;
//    }
//    else if (dim == 1)
//    {
//        std::cerr << "dim=1 is not working yet" << std::endl;
//        exit(-1);
//        for (unsigned int i = 0; i < coefficients.size(); i++)
//            volume += coefficients[i] / (0. + N - i);
//        volume *= 2.;
//    }
//    else
//    {
//        std::cerr << "Error in get_volume: dim=" << dim << std::endl;
//        exit(-1);
//    }
//    return volume;
//}
