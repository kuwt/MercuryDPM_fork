//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include <algorithm>
#include "NurbsUtils.h"
#include "Logger.h"

namespace NurbsUtils
{

bool isKnotVectorMonotonic(const std::vector<double>& knots)
{
    return std::is_sorted(knots.begin(), knots.end());
}

bool close(double a, double b, double eps)
{
    return (std::abs(a - b) < eps) ? true : false;
}

int findSpan(int degree, const std::vector<double>& knots, double u)
{
    // index of last control point
    int n = static_cast<int>(knots.size()) - degree - 2;
    
    // For u that is equal to last knot value
    if (close(u, knots[n + 1]))
    {
        return n;
    }
    
    // For values of u that lies outside the domain
    if (u > knots[n + 1])
    {
        return n;
    }
    if (u < knots[degree])
    {
        return degree;
    }
    
    // Binary search
    // TODO: Replace this with std::lower_bound
    int low = degree;
    int high = n + 1;
    int mid = (int) std::floor((low + high) / 2.0);
    while (u < knots[mid] || u >= knots[mid + 1])
    {
        if (u < knots[mid])
        {
            high = mid;
        }
        else
        {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return mid;
}

double bsplineOneBasis(int i, int deg, const std::vector<double>& U, double u)
{
    int m = static_cast<int>(U.size()) - 1;
    // Special case
    if ((i == 0 && close(u, U[0])) || (i == m - deg - 1 && close(u, U[m])))
    {
        return 1.0;
    }
    // Local Property
    if (u < U[i] || u >= U[i + deg + 1])
    {
        return 0.0;
    }
    // Initialize zeroth-degree functions
    std::vector<double> N;
    N.resize(deg + 1);
    for (int j = 0; j <= deg; j++)
    {
        N[j] = (u >= U[i + j] && u < U[i + j + 1]) ? 1.0 : 0.0;
    }
    // Compute triangular table
    for (int k = 1; k <= deg; k++)
    {
        double saved = (close(N[0], 0.0)) ? 0.0
                                          : ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
        for (int j = 0; j < deg - k + 1; j++)
        {
            double Uleft = U[i + j + 1];
            double Uright = U[i + j + k + 1];
            if (close(N[j + 1], 0.0))
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                double temp = N[j + 1] / (Uright - Uleft);
                N[j] = saved + (Uright - u) * temp;
                saved = (u - Uleft) * temp;
            }
        }
    }
    return N[0];
}

void bsplineBasis(int deg, int span, const std::vector<double>& knots, double u,
                  std::vector<double>& N)
{
    N.clear();
    N.resize(deg + 1, 0.0);
    std::vector<double> left, right;
    left.resize(deg + 1, 0.0);
    right.resize(deg + 1, 0.0);
    
    N[0] = 1.0;
    
    for (int j = 1; j <= deg; j++)
    {
        left[j] = (u - knots[span + 1 - j]);
        right[j] = knots[span + j] - u;
        Mdouble saved = 0.0;
        for (int r = 0; r < j; r++)
        {
            const Mdouble temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}
    
}
