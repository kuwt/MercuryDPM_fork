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
#ifndef Polynomial_H
#define Polynomial_H

#include <GeneralDefine.h>
#include "IntegralType.h"
#include <iostream>
#include <vector>

class BaseParticle;

class BaseInteraction;

class DPMBase;

class Coordinates;

namespace CGFunctions
{

/*!
 * \brief PolynomialType is used to define how files are opened random fixed-particle bottom
 * \todo add user-defined PolynomialType USER
 */
enum class PolynomialType : unsigned char
{
    HEAVISIDE = 0,
    LINEAR = 1,
    LUCY = 2
};

std::ostream& operator<<(std::ostream& os, PolynomialType type);

std::istream& operator>>(std::istream& is, PolynomialType& type);

/*!
 * \brief Defines the position of the CGPoint (e.g. x, y, z) and the parameters 
 * of a polynomial coarse-graining function (width and cutoff).
 * \details The class is derived from a CGCoordinates class, thus it contains 
 * the position of the CGPoint in the non-averaged directions. 
 * It further contains the parameters of a polynomial coarse-graining function 
 * (coefficients and cutoff), as well as internal variables (not user-defined, 
 * but set by the code), like the integralPrefactor. 
 * It also contains all functions that only depend on the coordinate and 
 * function type; e.g. evaluateCGFunction and evaluateCGIntegral. 
 * 
 * The polynomialType_ is used to define what polynomial is actually used; 
 * currently only Heaviside, linear and Lucy polynomials are possible. 
 * 
 * To simplify notation, the Heaviside, Linear and Lucy classed are derived from 
 * Polynomial, which set the polynomial type explicitly. Thus, the Polynomial 
 * class should not be used directly. 
 * 
 * See CGFunctions and Gauss for more details.
 */
template<class Coordinates>
class Polynomial
{
public:
    
    typedef Coordinates CoordinatesType;
    
    /*!
     * \brief Default constructor, sets all parameters to zero.
     */
    Polynomial();
    
    /*!
     * \brief Copy constructor. It copies all objects the class contains.
     * \param[in] p the class that has to be copied
     */
    Polynomial(const Polynomial& p) = default;
    
    /*!
     * \brief Destructor, it simply destructs the PolynomialCoordinates and all the objects it contains.
     */
    ~Polynomial() = default;
    
    /*!
     * \brief Writes class content into an output stream, usually a stat file
     * \param[out] os output stream
     */
    void write(std::ostream& os) const;
    
    /*!
     * \brief
     */
    void setPolynomialType(PolynomialType polynomialType);
    
    /*!
     * \brief Set the cutoff radius
     */
    void setWidth(Mdouble width);
    
    Mdouble getWidth() const;
    
    /*!
     * \brief Set the standard deviation
     */
    void setStandardDeviation(Mdouble std);
    
    /*!
     * \brief
     */
    void setCutoff(Mdouble cutoff);
    
    /*!
     * \brief
     */
    Mdouble getCutoff() const;
    
    /*!
     * \brief
     */
    void computeCoefficients();
    
    /*!
     * \brief Evaluates the coarse-graining function
     */
    Mdouble evaluateCGFunction(const Vec3D& position, const Coordinates& r);
    
    /*!
     * \brief Evaluates the line integral needed for the calculation of stresses.
     */
    Mdouble
    evaluateCGIntegral(const BaseInteraction& i, const Coordinates& r, IntegralType type = IntegralType::I_TO_P);
    
    /*!
     * \brief Evaluates the line integral needed for the calculation of stresses
     * for 1D CGCoordinates.
     */
    Mdouble
    evaluateCGIntegral1D(const BaseInteraction& i, const Coordinates& r, IntegralType type = IntegralType::I_TO_P);
    
    std::vector<Mdouble> getCoefficients();
    
    Vec3D evaluateCGFunctionDerivatives(const Vec3D& position, const Coordinates& r);
    
    Mdouble evaluateCGFunctionDerivativeWithFD(const Vec3D& position, const Coordinates& r, const int i);
    
    /*!
     * \brief Returns the finite difference step size used to evaluate derivatives
     * of the CG function
     */
    Mdouble getEps() const;
    
    /*!
     * \brief Sets the finite difference step size used to evaluate derivatives
     * of the CG function.
     */
    void setEps(Mdouble eps);

protected:
    
    /*!
     * Stores the polynomial coefficients,starting with the highest order
     * \f$c_i\f$ of \f$\sum_i c_i (|pos|/c)^(n-i)\f$.
     * The coefficients are normalised such that
     * \f$\int \phi(\vec r,\vec r_i) d\vec r = 1\f$;
     * see Polynomial<Coordinates>::evaluateCGFunction for details.
     *
     * Internal variables that should not be set by the user.
     * \todo Make variables internal.
     */
    std::vector<Mdouble> coefficients_;
    
    /*!
     * The polynomial type. The user should choose the Polynomial type by
     * choosing one of the classes derived from Polynomial:
     * Heaviside, Linear, or Lucy
     */
    PolynomialType polynomialType_;
    
    Mdouble normalLength_;
    
    Vec3D normal_;
    
    unsigned currentInteraction_;
    
    /*!
     * cutoff radius of the polynomial cg function, beyond which the function is
     * zero.
     */
    Mdouble cutoff_;
    
    /*!
     * Finite Difference step size, used to compute the derivative of the CG function
     */
    Mdouble eps_;
    
};
    
} //namespace CGFunctions
#include "Polynomial.hcc"

#endif
