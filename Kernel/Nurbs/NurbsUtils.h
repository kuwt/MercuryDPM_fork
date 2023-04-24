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

#ifndef MERCURYDPM_NURBSUTILS_H
#define MERCURYDPM_NURBSUTILS_H

#include <vector>
#include <algorithm>
#include "Math/Vector.h"

namespace NurbsUtils
{


/**
// A simple class for representing 2D runtime arrays.
*/
template <typename T>
class array2 {
public:
    array2(size_t nRows, size_t nCols, T fillValue = 0.0)
    : rows(nRows), cols(nCols) {
        data.resize(rows * cols, fillValue);
    }

    T operator()(size_t row, size_t col) const {
        return data[row*cols + col];
    }

    T& operator()(size_t row, size_t col) {
        return data[row*cols + col];
    }

private:
    size_t rows, cols;
    std::vector<T> data;
};

bool isKnotVectorMonotonic(const std::vector<double>& knots);


bool close(double a, double b, double eps = std::numeric_limits<double>::epsilon());

/**
Find the span of the given parameter in the knot vector.
@param[in] degree Degree of the curve.
@param[in] knots Knot vector of the curve.
@param[in] u Parameter value.
@return Span index into the knot vector such that (span - 1) < u <= span
*/
int findSpan(int degree, const std::vector<double>& knots, double u);


/**
Compute a single B-spline basis function
@param[in] i The ith basis function to compute.
@param[in] deg Degree of the basis function.
@param[in] U Knot vector corresponding to the basis functions.
@param[in] u Parameter to evaluate the basis functions at.
@return The value of the ith basis function at u.
*/
double bsplineOneBasis(int i, int deg, const std::vector<double>& U, double u);

/**
// Compute all non-zero B-spline basis functions
@param[in] deg Degree of the basis function.
@param[in] span Index obtained from findSpan() corresponding the u and knots.
@param[in] knots Knot vector corresponding to the basis functions.
@param[in] u Parameter to evaluate the basis functions at.
@param[in, out] N Values of (deg+1) non-zero basis functions.
*/
void bsplineBasis(int deg, int span, const std::vector<double>& knots, double u, std::vector<double>& N);

/**
// Compute all non-zero derivatives of B-spline basis functions
@param[in] deg Degree of the basis function.
@param[in] span Index obtained from findSpan() corresponding the u and knots.
@param[in] knots Knot vector corresponding to the basis functions.
@param[in] u Parameter to evaluate the basis functions at.
@param[in] nDers Number of derivatives to compute (nDers <= deg)
@param[in, out] ders Values of non-zero derivatives of basis functions.
*/
void bsplineDerBasis(int deg, int span, const std::vector<double>& knots, double u,
                     int nDers, std::vector<std::vector<double>>& ders);

/*!
 * \brief Creates a uniform (clamped) knot vector
 * @param numberOfControlPoints The number of control points the knot vector should correspond to
 * @param degree The degree the knot vector should correspond to
 * @param clampedAtStart Whether or not the knot vector should be clamped at the start
 * @param clampedAtEnd Whether or not the knot vector should be clamped at the end
 * @return A uniform (clamped) knot vector
 */
std::vector<Mdouble> createUniformKnotVector(unsigned int numberOfControlPoints, unsigned int degree, bool clampedAtStart, bool clampedAtEnd);

/*!
 * \brief Resets the knot vector to the interval [0, 1]
 * @param knots The knot vector to be normalized
 */
void normalizeKnotVector(std::vector<Mdouble>& knots);

/*!
 * \brief Extends the knot vector for when control points have been added at the start or end
 * @param knots The knot vector to be extended
 * @param degree The nurbs degree
 * @param numStart Number to be added at the start
 * @param numEnd Number to be added at the end
 * @param forceBothEndsUniform When only adding to one end, the other end can remain untouched or it can be forced to be made uniform
 */
void extendKnotVector(std::vector<Mdouble>& knots, unsigned int degree, unsigned int numStart, unsigned int numEnd, bool forceBothEndsUniform = false);

/*!
 * \brief Evaluate point on a NURBS surface
 * @param u Parameter to evaluate the surface at
 * @param v Parameter to evaluate the surface at
 * @param knotsU Knot vector in u-direction
 * @param knotsV Knot vector in v-direction
 * @param controlPoints 2D vector of control points
 * @param weights 2D vector of weights
 * @return Resulting point on the surface at (u, v)
 */
Vec3D evaluate(Mdouble u, Mdouble v, std::vector<Mdouble> knotsU, std::vector<Mdouble> knotsV,
               std::vector<std::vector<Vec3D>> controlPoints, std::vector<std::vector<Mdouble>> weights);
}

#endif //MERCURYDPM_NURBSUTILS_H
