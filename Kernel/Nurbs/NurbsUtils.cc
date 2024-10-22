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

void bsplineDerBasis(int deg, int span, const std::vector<double>& knots, double u,
                     int nDers, std::vector<std::vector<double>> &ders) {

    std::vector<double> left, right;
    left.resize(deg + 1, 0.0);
    right.resize(deg + 1, 0.0);
    double saved = 0.0, temp = 0.0;

    array2<double> ndu(deg + 1, deg + 1);
    ndu(0, 0) = 1.0;

    for (int j = 1; j <= deg; j++) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        saved = 0.0;

        for (int r = 0; r < j; r++) {
            // Lower triangle
            ndu(j, r) = right[r + 1] + left[j - r];
            temp = ndu(r, j - 1) / ndu(j, r);
            // Upper triangle
            ndu(r, j) = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        ndu(j, j) = saved;
    }

    ders.clear();
    ders.resize(nDers + 1);
    for (int i = 0; i < ders.size(); i++) {
        ders[i].resize(deg + 1, 0.0);
    }

    for (int j = 0; j <= deg; j++) {
        ders[0][j] = ndu(j, deg);
    }

    array2<double> a(2, deg + 1);

    for (int r = 0; r <= deg; r++) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;

        for (int k = 1; k <= nDers; k++) {
            double d = 0.0;
            int rk = r - k;
            int pk = deg - k;
            int j1 = 0;
            int j2 = 0;

            if (r >= k) {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }

            if (rk >= -1) {
                j1 = 1;
            }
            else {
                j1 = -rk;
            }

            if (r - 1 <= pk) {
                j2 = k - 1;
            }
            else {
                j2 = deg - r;
            }

            for (int j = j1; j <= j2; j++) {
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }

            if (r <= pk) {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }


            ders[k][r] = d;

            int temp = s1;
            s1 = s2;
            s2 = temp;
        }
    }

    double fac = static_cast<double>(deg);
    for (int k = 1; k <= nDers; k++) {
        for (int j = 0; j <= deg; j++) {
            ders[k][j] *= fac;
        }
        fac *= static_cast<double>(deg - k);
    }
}

std::vector<Mdouble> createUniformKnotVector(unsigned int numberOfControlPoints, unsigned int degree, bool clampedAtStart, bool clampedAtEnd)
{
    // Total number of knots
    unsigned int numberOfKnots = numberOfControlPoints + degree + 1;
    std::vector<Mdouble> knots;
    knots.reserve(numberOfKnots);
    
    // Create uniform from 0 to 1
    for (int i = 0; i < numberOfKnots; i++)
    {
        knots.push_back(i / static_cast<Mdouble>(numberOfKnots - 1));
    }
    
    // When clamped at start, first degree+1 values are set to 0.
    if (clampedAtStart)
    {
        for (int i = 0; i <= degree; i++)
        {
            knots[i] = 0.0;
        }
    }
    
    // When clamped at end, last degree+1 values are set to 1
    if (clampedAtEnd)
    {
        for (int i = 0; i <= degree; i++)
        {
            knots[knots.size() - 1 - i] = 1.0;
        }
    }
    
    return knots;
}

void normalizeKnotVector(std::vector<Mdouble>& knots)
{
    // Reset the interval to [0,1]
    const Mdouble minK = knots.front();
    const Mdouble maxK = knots.back();
    
    // Already in normalized form
    if (close(minK, 0.0) && close(maxK, 1.0))
        return;
    
    for (Mdouble& k : knots)
    {
        k = (k - minK) / (maxK - minK);
    }
}

void extendKnotVector(std::vector<Mdouble>& knots, unsigned int degree, unsigned int numStart, unsigned int numEnd, bool forceBothEndsUniform)
{
    // Usually it should already be in normalized form, but just to be sure.
    // Needed because otherwise comparing to uniform step size doesn't work.
    normalizeKnotVector(knots);
    
    Mdouble uniformStepSize = 1.0 / (knots.size() - 1);
    
    // This only works properly when the original start and end knots are uniformly spaced.
    // Therefore force the degree+1 number of knots at the start and end to have uniform step size
    // and remember if they weren't already uniform, so when needed a message can be logged.
    
    // Since the original knot vector might change a bit, only do stuff when there are actually knots to be added.
    if (numStart > 0 || forceBothEndsUniform)
    {
        // Check if first degree+1 knots are uniform.
        // Update: also an additional number of knots added to the end should be uniform (the ones "wrapping around")
        bool startsOfUniform = true;
        for (int i = 0; i <= degree + numEnd; i++)
        {
            if (!close(knots[i], i * uniformStepSize))
            {
                startsOfUniform = false;
                knots[i] = i * uniformStepSize;
            }
        }
    
        // When needed, let user know that original knot vector has been edited.
        if (!startsOfUniform)
        {
            logger(INFO, "Start of knot vector (first degree+1 values) has been changed to be uniform. The shape might be slightly affected. ");
        }
    
        // Insert knots at start
        for (int i = 1; i <= numStart; i++)
        {
            knots.insert(knots.begin(), -i * uniformStepSize);
        }
    }
    
    // Since the original knot vector might change a bit, only do stuff when there are actually knots to be added.
    if (numEnd > 0 || forceBothEndsUniform)
    {
        // Check if last degree+1 knots are uniform.
        // Update: also an additional number of knots added to the start should be uniform (the ones "wrapping around")
        bool endsOfUniform = true;
        for (int i = 0; i <= degree + numStart; i++)
        {
            if (!close(knots[knots.size() - 1 - i], (1.0 - i * uniformStepSize)))
            {
                endsOfUniform = false;
                knots[knots.size() - 1 - i] = 1.0 - i * uniformStepSize;
            }
        }
    
        // When needed, let user know that original knot vector has been edited.
        if (!endsOfUniform)
        {
            logger(INFO, "End of knot vector (last degree+1 values) has been changed to be uniform. The shape might be slightly affected. ");
        }
    
        // Add knots to end
        for (int i = 1; i <= numEnd; i++)
        {
            knots.push_back(1.0 + i * uniformStepSize);
        }
    }
    
    // Reset to interval [0, 1]
    normalizeKnotVector(knots);
}

Vec3D evaluate(Mdouble u, Mdouble v, std::vector<Mdouble> knotsU, std::vector<Mdouble> knotsV,
               std::vector<std::vector<Vec3D>> controlPoints, std::vector<std::vector<Mdouble>> weights)
{
    unsigned int degreeU = knotsU.size() - controlPoints.size() - 1;
    unsigned int degreeV = knotsV.size() - controlPoints[0].size() - 1;
    
    Vec3D point = {0,0,0};
    double temp = 0;
    
    // Find span and non-zero basis functions
    int spanU = findSpan(degreeU, knotsU, u);
    int spanV = findSpan(degreeV, knotsV, v);
    std::vector<double> Nu, Nv;
    bsplineBasis(degreeU, spanU, knotsU, u, Nu);
    bsplineBasis(degreeV, spanV, knotsV, v, Nv);
    // make linear combination
    for (int l = 0; l <= degreeV; ++l)
    {
        for (int k = 0; k <= degreeU; ++k)
        {
            double weight = Nv[l] * Nu[k] * weights[spanU - degreeU + k][spanV - degreeV + l];
            point += weight * controlPoints[spanU - degreeU + k][spanV - degreeV + l];
            temp += weight;
        }
    }
    point /= temp;
    return point;
}

}
