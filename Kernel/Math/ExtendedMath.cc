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

#include "NumericalVector.h"
#include "ExtendedMath.h"
#include "Quaternion.h"
#include <cmath>
#include "Logger.h"
////This has the defintion of quiet nan
//#include <limits>
//#include <iostream>
//#include <sys/stat.h>
#include <iomanip>
//#include <cmath>
//#include <sstream>
//#include <cstdlib>
//#include <limits>
//#include <string>


// sin(x) = x - x^3/3! + x^5/5! - x^7/7! + ...
Mdouble mathsFunc::sin(Mdouble x)
{
    Mdouble N = floor(x / (2.0 * constants::pi) + 0.5);
    x -= N * 2.0 * constants::pi;
    Mdouble Sum = 0;
    Mdouble Power = x;
    Mdouble Sign = 1;
    const Mdouble x2 = x * x;
    Mdouble Fact = 1.0;
    for (unsigned int i = 1; i < 25; i += 2)
    {
        Sum += Sign * Power / Fact;
        Power *= x2;
        Fact *= (i + 1) * (i + 2);
        Sign *= -1.0;
    }
    return Sum;
}

// cos(x) = 1 - x^2/2! + x^4/4! - x^6/6! + ...
Mdouble mathsFunc::cos(Mdouble x)
{
    Mdouble N = floor(x / (2.0 * constants::pi) + 0.5);
    x -= N * 2.0 * constants::pi;
    Mdouble Sum = 1.0;
    Mdouble Power = 1;
    Mdouble Sign = 1;
    const Mdouble x2 = x * x;
    Mdouble Fact = 1.0;
    for (unsigned int i = 2; i < 25; i += 2)
    {
        Power *= x2;
        Fact *= i * (i - 1);
        Sign *= -1.0;
        Sum += Sign * Power / Fact;
    }
    return Sum;
}

//from: http://www.codeproject.com/Tips/311714/Natural-Logarithms-and-Exponent
Mdouble mathsFunc::exp(Mdouble Exponent)
{
    Mdouble X, P, Frac, I, L;
    X = Exponent;
    Frac = X;
    P = (1.0 + X);
    I = 1.0;
    
    do
    {
        I++;
        Frac *= (X / I);
        L = P;
        P += Frac;
    } while (L != P);
    
    return P;
}

///\todo check if this function works
Mdouble mathsFunc::log(Mdouble Power)
{
    Mdouble N, P, L, R, A, E;
    E = 2.71828182845905;
    P = Power;
    N = 0.0;
    
    // This speeds up the convergence by calculating the integral
    while (P >= E)
    {
        P /= E;
        N++;
    }
    N += (P / E);
    P = Power;
    do
    {
        A = N;
        L = (P / (exp(N - 1.0)));
        R = ((N - 1.0) * E);
        N = ((L + R) / E);
    } while (N < A);
    //} while (N != A);
    
    return N;
}

/**
 * This is the gamma function, gives 'exact' answers for the half integer values
 * This is done using the recussion relation and the known values for 1 and 0.5
 * Note, return NaN for non-half integer values
 */

Mdouble mathsFunc::gamma(Mdouble gamma_in)
{
    const Mdouble ep = 1e-5;
    
    if (gamma_in > 1.0 + ep)
    {
        return ((gamma_in - 1) * gamma(gamma_in - 1));
    }
    else
    {
        
        if ((gamma_in - ep < 1.0) && (gamma_in + ep > 1.0))
            return 1.0;
        else if ((gamma_in - ep < 0.5) && (gamma_in + ep > 0.5))
            return constants::sqrt_pi;
        else
            return std::numeric_limits<Mdouble>::quiet_NaN();
    }
} //end func gamma

/*!
 * Computes the beta function for real numbers, based on cmath's approximation of the ln of the gamma function.
 * See https://en.wikipedia.org/wiki/Beta_function for more details
 * \param z first Mdouble argument to compute the beta function of
 * \param w second Mdouble argument to compute the beta function of
 * \return  the value of beta(z,w) as an Mdouble
 */
Mdouble mathsFunc::beta(Mdouble z, Mdouble w)
{
    return std::exp(std::lgamma(z) + std::lgamma(w) - std::lgamma(z + w));
}

/**
 * This is a chi_squared function return the value x and degrees of freedom k
 */
Mdouble mathsFunc::chi_squared(const Mdouble x, const unsigned int k)
{
    
    Mdouble prefactor = pow(2, k / 2.0) * gamma(k / 2.0);
    Mdouble mainfactor = pow(x, k / 2.0 - 1) * exp(x / -2.0);
    
    return mainfactor / prefactor;
    
}

/**This calculates the probability based on a chi squared test
 * First we calculated the  cumulative chi_squared function.
 * This is the function which actually gives the probability back
 * It is calculated by calling the normal chi_squared function and using the trapezoidal rule.
 * The final results is 1-the cumulative chi_squared function
 */
Mdouble mathsFunc::chi_squared_prob(const Mdouble x_max, const unsigned int k)
{

//The current value was picked by tried were it stopped effect the 4 d.p.
    const int num_steps_per_unit = 100;
    Mdouble sum = 0;
    Mdouble x = 0;
    long int num_steps = static_cast<int>(num_steps_per_unit * x_max);
//Use trapezional rule, but ignoring the ends
    for (int i = 0; i < num_steps; i++)
    {
        x = x_max / num_steps * (i + 0.5);
        sum = sum + chi_squared(x, k);
    }
    return 1.0 - sum * x_max / num_steps;
    
}

Mdouble mathsFunc::goldenSectionSearch(Mdouble(* function)(const Mdouble), Mdouble min, Mdouble cur, Mdouble max,
                                       Mdouble endCondition, Mdouble curVal)
{
    if (std::abs(max - min) < endCondition)
    {
        return 0.5 * (min + max);
    }
    std::cout << "Min=" << min << " Max=" << max << " diff=" << max - min << std::endl;
    Mdouble resphi = 2 - 0.5 * (1 + std::sqrt(5));
    Mdouble x;
    if (max - cur > cur - min)
    {
        x = cur + resphi * (max - cur);
    }
    else
    {
        x = cur - resphi * (cur - min);
    }
    if (std::isnan(curVal))
        curVal = function(cur);
    Mdouble xVal = function(x);
    if (xVal < curVal)
    {
        if (max - cur > cur - min)
        {
            return goldenSectionSearch(function, cur, x, max, endCondition, xVal);
        }
        else
        {
            return goldenSectionSearch(function, min, x, cur, endCondition, xVal);
        }
    }
    else
    {
        if (max - cur > cur - min)
        {
            return goldenSectionSearch(function, min, cur, x, endCondition, curVal);
        }
        else
        {
            return goldenSectionSearch(function, x, cur, max, endCondition, curVal);
        }
    }
}

bool mathsFunc::isEqual(Mdouble v1, Mdouble v2, Mdouble absError)
{
    return std::abs(v1 - v2) <= absError;
}


bool mathsFunc::isEqual(Vec3D v1, Vec3D v2, Mdouble absError)
{
    return isEqual(v1.X, v2.X, absError) && isEqual(v1.Y, v2.Y, absError) && isEqual(v1.Z, v2.Z, absError);
}

bool mathsFunc::isEqual(Matrix3D m1, Matrix3D m2, Mdouble absError)
{
    return (isEqual(m1.XX, m2.XX, absError)
            && isEqual(m1.XY, m2.XY, absError)
            && isEqual(m1.XZ, m2.XZ, absError)
            && isEqual(m1.YX, m2.YX, absError)
            && isEqual(m1.YY, m2.YY, absError)
            && isEqual(m1.YZ, m2.YZ, absError)
            && isEqual(m1.ZX, m2.ZX, absError)
            && isEqual(m1.ZY, m2.ZY, absError)
            && isEqual(m1.ZZ, m2.ZZ, absError));
}

bool mathsFunc::isEqual(MatrixSymmetric3D m1, MatrixSymmetric3D m2, Mdouble absError)
{
    return (isEqual(m1.XX, m2.XX, absError)
            && isEqual(m1.XY, m2.XY, absError)
            && isEqual(m1.XZ, m2.XZ, absError)
            && isEqual(m1.YY, m2.YY, absError)
            && isEqual(m1.YZ, m2.YZ, absError)
            && isEqual(m1.ZZ, m2.ZZ, absError));
}

bool mathsFunc::isEqual(Quaternion v1, Quaternion v2, double absError)
{
    return isEqual(v1.getComponent(0), v2.getComponent(0), absError) &&
           isEqual(v1.getComponent(1), v2.getComponent(1), absError) &&
           isEqual(v1.getComponent(2), v2.getComponent(2), absError) &&
           isEqual(v1.getComponent(3), v2.getComponent(3), absError);
}

Mdouble mathsFunc::chebyshev(Mdouble x, const Mdouble coef[], int N)
{
    const Mdouble* p = coef;
    Mdouble b0 = *p++;
    Mdouble b1 = 0, b2;
    int i = N - 1;
    
    logger.assert(i > 0, "i is greater than 0");
    do
    {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + *p++;
    } while (--i);
    
    return (0.5 * (b0 - b2));
}

Mdouble mathsFunc::I0_exp(Mdouble x)
{
    // Coefficients for [0..8]
    const Mdouble A[] =
            {
                    -4.415341646479339379501E-18,
                    3.330794518822238097831E-17,
                    -2.431279846547954693591E-16,
                    1.715391285555133030611E-15,
                    -1.168533287799345168081E-14,
                    7.676185498604935616881E-14,
                    -4.856446783111929460901E-13,
                    2.955052663129639834611E-12,
                    -1.726826291441555707231E-11,
                    9.675809035373236912241E-11,
                    -5.189795601635262906661E-10,
                    2.659823724682386650351E-9,
                    -1.300025009986248042121E-8,
                    6.046995022541918949321E-8,
                    -2.670793853940611733911E-7,
                    1.117387539120103718151E-6,
                    -4.416738358458750563591E-6,
                    1.644844807072889708931E-5,
                    -5.754195010082103703981E-5,
                    1.885028850958416557291E-4,
                    -5.763755745385823658851E-4,
                    1.639475616941335798421E-3,
                    -4.324309995050575944301E-3,
                    1.054646039459499831831E-2,
                    -2.373741480589946881561E-2,
                    4.930528423967070848781E-2,
                    -9.490109704804764442101E-2,
                    1.716209015222087753491E-1,
                    -3.046826723431983986831E-1,
                    6.767952744094760849951E-1
            };
    
    // Coefficients for [8..infinity]
    const Mdouble B[] =
            {
                    -7.233180487874753954561E-18,
                    -4.830504485944182071261E-18,
                    4.465621420296759999011E-17,
                    3.461222867697461093101E-17,
                    -2.827623980516583484941E-16,
                    -3.425485619677219134621E-16,
                    1.772560133056526383601E-15,
                    3.811680669352622420751E-15,
                    -9.554846698828307648701E-15,
                    -4.150569347287222086631E-14,
                    1.540086217521409826911E-14,
                    3.852778382742142701141E-13,
                    7.180124451383666233671E-13,
                    -1.794178531506806117781E-12,
                    -1.321581184044771311881E-11,
                    -3.149916527963241364541E-11,
                    1.188914710784643834241E-11,
                    4.940602388224969589101E-10,
                    3.396232025708386345151E-9,
                    2.266668990498178064591E-8,
                    2.048918589469063741831E-7,
                    2.891370520834756482971E-6,
                    6.889758346916823984261E-5,
                    3.369116478255694089901E-3,
                    8.044904110141088316081E-1
            };
    
    if (x < 0)
        x = -x;
    
    if (x <= 8.0)
    {
        Mdouble y = (x / 2.0) - 2.0;
        return (chebyshev(y, A, 30));
    }
    
    return (chebyshev(32.0 / x - 2.0, B, 25) / sqrt(x));
    
}

Mdouble mathsFunc::I0(Mdouble x)
{
    if (x < 0)
        x = -x;
    return exp(x) * I0_exp(x);
}

// This function computes all n order and positive order m associated Legendre polynomials using recursive formulations
// The polynomials are evaluated at x.
// They are organised as follows: P_0^0, P_1^0, P_1^1, P_2^0, P_2^1, ...
NumericalVector<> sphericalHarmonics::associatedLegendrePolynomials(int n, Mdouble x)
{
    //Given n and m, we only have to compute P_n^(|m|)
    //The function will return all these P values for theta
    std::size_t nTerms = 0.5 * (n + 1) * (n + 2);
    NumericalVector<> polynomials(nTerms);
    
    size_t location_current;
    size_t location_previous;
    Mdouble temp;
    
    polynomials(0) = 1; //P_0^0 = 1;
    for (int l = 1; l <= n; l++)
    {
        //first compute P_l^l
        location_current = 0.5 * (l + 1) * (l + 2) - 1;
        location_previous = location_current - (l + 1);
        polynomials(location_current) = -(2.0 * (l - 1.0) + 1.0) * std::sqrt(1.0 - x * x) *
                                        polynomials(location_previous); // Recursive formula from wiki
        
        //second, compute P_l^(l-1) based on P_(l-1)^(l-1)
        polynomials(location_current - 1) =
                x * (2.0 * (l - 1) + 1) * polynomials(location_previous); // Recursive formula from wiki
    }
    
    //thirdly, compute other values
    for (int l = 2; l <= n; l++)
    {
        for (int m = (l - 2); m >= 0; m--)
        {
            location_current = (0.5 * (l + 1) * (l + 2) - 1) - l + m;
            temp = polynomials(location_current + 2) +
                   2.0 * (m + 1) * x / sqrt(1 - x * x) * polynomials(location_current + 1);
            polynomials(location_current) = temp / ((m - l) * (l + m + 1)); // variation on greengard eqn 3.34 from wiki
        }
    }
    
    
    return polynomials;
}


// This function computes up to p order spherical harmonics as function of theta and phi
// They are organised as follows: Y_0^0, Y_1^-1, Y_1^0, Y_1^1, Y_2^-2, ...
NumericalVector<std::complex<Mdouble>> sphericalHarmonics::sphericalHarmonics(int p, Mdouble theta, Mdouble phi)
{
    std::size_t nTerms = 0.5 * (p + 1) * (2 * p + 2);
    NumericalVector<std::complex<Mdouble>> Y(nTerms);
    NumericalVector<> polynomials = associatedLegendrePolynomials(p, std::cos(theta));
    
    //Compute the spherical harmonics
    for (int n = 0; n <= p; n++)
    {
        for (int mt = -n; mt <= n; mt++)
        {
            Mdouble m = mt;
            Mdouble m_abs = std::abs(mt);
            std::size_t location_current = n * n + (m + n); //n^2 is begin of Y_n^-n
            std::size_t location_polynomial = 0.5 * n * (n + 1) + m_abs;
            int fact1 = mathsFunc::factorial(n - m_abs);
            int fact2 = mathsFunc::factorial(n + m_abs);
            Mdouble fact = 1.0 * fact1 / fact2;
            std::complex<Mdouble> value =
                    std::sqrt(fact) * polynomials(location_polynomial) * std::exp(constants::i * m * phi);
            Y(location_current) = value;
        }
    }
    return Y;
}

//Compute all squaredFactorials (see eqn 5.23 in a short course on fast multipole methods) up to order p
// They are organised as follows: A_0^0, A_1^-1, A_1^0, A_1^1, A_2^-2, ...
NumericalVector<> sphericalHarmonics::computeSquaredFactorialValues(int p)
{
    std::size_t nTerms = 0.5 * (p + 1) * (2 * p + 2);
    NumericalVector<> squaredFactorials(nTerms);
    
    for (int n = 0; n <= p; n++)
    {
        for (int m = -n; m <= n; m++)
        {
            std::size_t location = n * n + (m + n); //(n^2 is begin of Y_n)
            int fact1 = mathsFunc::factorial(n - m);
            int fact2 = mathsFunc::factorial(n + m);
            Mdouble fact = fact1 * fact2;
            squaredFactorials(location) = std::pow(-1.0, n) / std::sqrt(fact);
        }
    }
    
    return squaredFactorials;
}



