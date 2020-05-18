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

#include "Base_XY_XZ_YZ.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"

using namespace CGCoordinates;

Mdouble Base_XY_XZ_YZ::getGaussPrefactor(Mdouble width, Mdouble cutoff)
{
    //Wolfram alpha: integrate(x*exp(-x^2/(2w^2)),{x,0,c})/integrate(x*exp(-x^2/(2w^2)),{x,0,inf})=1-e^(-c^2/(2 w^2))
    Mdouble prefactor = 1.0 / (constants::sqrt_2 * constants::sqrt_pi * width);
    return mathsFunc::square(prefactor) / (1.0 - exp(-0.5 * mathsFunc::square(cutoff / width)));
}

Mdouble Base_XY_XZ_YZ::getGaussIntegralPrefactor(Mdouble distance, Mdouble width, Mdouble cutoff)
{
    Mdouble widthSqrt2 = width * constants::sqrt_2;
    Mdouble a = -cutoff;
    Mdouble b = cutoff + distance;
    //1D prefactor
    Mdouble prefactor_ = 1.0 / (widthSqrt2 * constants::sqrt_pi);
    prefactor_ /= erf(cutoff / (widthSqrt2));
    return prefactor_ * 0.5 / (
            +erf(b / widthSqrt2) * b
            + widthSqrt2 / constants::sqrt_pi * exp(-mathsFunc::square(b / widthSqrt2))
            - erf(a / widthSqrt2) * a
            - widthSqrt2 / constants::sqrt_pi * exp(-mathsFunc::square(a / widthSqrt2))
    );
}

/*!
 * \details The volume is computed as
 * \f[volume=\int_0^1\sum_{i=1}^n c_i r^i 2 pi r dr = 2 pi \sum_{i=1}^n c_i/(i+2) \f]
 * with 2 pi r the circumference of a circle.
 */
void Base_XY_XZ_YZ::normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff)
{
    Mdouble volume = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); i++)
        volume += coefficients[i] / static_cast<Mdouble>(i + 2);
    volume *= 2.0 * constants::pi * mathsFunc::square(cutoff);
    for (double& coefficient : coefficients)
        coefficient /= volume;
}

const unsigned Base_XY_XZ_YZ::countVariables()
{
    return 2;
}
