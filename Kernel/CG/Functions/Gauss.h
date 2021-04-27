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
#ifndef Gauss_H
#define Gauss_H

#include <GeneralDefine.h>
#include <iostream>
#include "IntegralType.h"
#include "CG/Coordinates/R.h"
#include "CG/Coordinates/X.h"
#include "CG/Coordinates/Y.h"
#include "CG/Coordinates/Z.h"
#include "CG/Coordinates/RZ.h"
#include "CG/Coordinates/XY.h"
#include "CG/Coordinates/XZ.h"
#include "CG/Coordinates/YZ.h"
#include "CG/Coordinates/XYZ.h"
#include "CG/Coordinates/O.h"

class BaseParticle;

class BaseInteraction;

class DPMBase;

class Coordinates;

/*!
 * \namespace CGFunctions
 * \brief Contains base classes of CGPoint; CGPoint is always templated with one 
 * of these classes; these classes contain the position of the CGPoint and the 
 * parameters of the coarse-graining function (width, cutoff, ...).
 * \details See member class Gauss for more details.
 */
namespace CGFunctions
{

/*!
 * \class Gauss 
 * \brief Defines the position of the CGPoint (e.g. x, y, z) and the parameters 
 * of the Gauss coarse-graining function (width and cutoff).
 * \details The class is derived from a CGCoordinates class, thus it contains 
 * the position of the CGPoint in the non-averaged directions. It also contains 
 * the parameters of the Gauss coarse-graining function (width and cutoff), as 
 * well as internal variables (not user-defined, but set by the code), like the 
 * prefactor of the Gauss function (which depends on width, cutoff, and the 
 * coordinate type). 
 * It also contains all functions that only depend on the coordinate and 
 * function type; e.g. evaluateCGFunction and evaluateCGIntegral. 
 * 
 * See CGFunctions for more details.
 * 
 * \class Coordinates
 * \brief Template argument; use a member class of CGCoordinates to instantiate.
 */
template<class Coordinates>
class Gauss
{
public:
    
    typedef Coordinates CoordinatesType;
    
    /*!
     * \brief Default constructor, it simply creates an empty GaussCoordinates.
     */
    Gauss();
    
    /*!
     * \brief Copy constructor. It copies the GaussCoordinates and all objects it contains.
     */
    Gauss(const Gauss& c) = default;
    
    /*!
     * \brief Destructor, does nothing, as no new'ed objects are used.
     */
    ~Gauss() = default;
    
    /*!
     * \brief Writes class content into an output stream, usually a stat file
     */
    void write(std::ostream& os) const;
    
    /*!
     * \brief
     */
    //void setStandardDeviation(Mdouble standardDeviation);
    
    /*!
     * \brief Sets the width of the coarse-graining function.
     */
    void setWidth(Mdouble width);
    
    /*!
     * \brief Sets the standard deviation of the coarse-graining function.
     */
    void setStandardDeviation(Mdouble std);
    
    /*!
     * \brief Sets the width and cutoff of the coarse-graining function.
     */
    void setWidthAndCutoff(Mdouble width, Mdouble cutoff);
    
    /*!
     * \brief Returns the width of the coarse-graining function.
     */
    Mdouble getWidth() const;
    
    /*!
     * \brief Returns the cutoff of the coarse-graining function.
     */
    Mdouble getCutoff() const;
    
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
    
    /*!
     * \brief Evaluates the coarse-graining function.
     */
    Mdouble evaluateCGFunction(const Vec3D& position, const Coordinates r);
    
    Vec3D evaluateCGFunctionDerivatives(const Vec3D& position, const Coordinates& r);
    
    Mdouble evaluateCGFunctionDerivativeWithFD(const Vec3D& position, const Coordinates& r, const int i);
    
    Mdouble evaluateCylindricalCGFunction(const Vec3D& position, const CGCoordinates::R r)
    { return constants::NaN; }
    
    /*!
     * \brief Evaluates the line integral needed for the calculation of stresses.
     */
    Mdouble
    evaluateCGIntegral(const BaseInteraction& i, const Coordinates r, IntegralType type = IntegralType::I_TO_P);
    
    Mdouble evaluateCylindricalCGIntegral(const BaseInteraction& c, const CGCoordinates::R r,
                                          IntegralType type = IntegralType::I_TO_P)
    { return constants::NaN; }
    
    /*!
     * \brief Evaluates the line integral needed for the calculation of stresses
     * for 1D CGCoordinates.
     */
    Mdouble
    evaluateCGIntegral1D(const BaseInteraction& i, const Coordinates r, IntegralType type = IntegralType::I_TO_P);

protected:
    
    /*!
     * The coarse-graining width.
     * It is 1.0 by default, unless modified by the user.
     * \todo TW I thought of implementing width_, cutoff_ and prefactor_ as
     * const Mdouble&, as theses parameters are not set per CGPoint,
     * but per CG object (i.e. all points share the same values).
     * However, this site seems to discourage this:
     * http://stackoverflow.com/questions/12387239/reference-member-variables-as-class-members
     * Anyone has an opinion on this? \@dducks \@thorntonar
     */
    Mdouble width_;
    
    /*!
     * The cutoff of the coarse-graining function.
     * It is 3.0 by default, unless modified by the user.
     */
    Mdouble cutoff_;
    
    /*!
     * The prefactor of the coarse-graining function.
     * This parameter is internal, thus cannot be set directly by the user.
     */
    Mdouble prefactor_;
    
    /*!
     * The prefactor of the coarse-graining integral.
     * Depends on the value of currentDistance_.
     * This parameter is internal, thus cannot be set directly by the user.
     */
    Mdouble integralPrefactor_;
    
    Mdouble normalLength_;
    
    Vec3D normal_;
    
    /*!
     * The length of the branch vector of the current interaction.
     * If this parameter is changed, the integralPrefactor_ has to be recomputed.
     * This parameter is internal, thus cannot be set directly by the user.
     */
    unsigned currentInteraction_;
    
    /*!
     * Finite Difference step size, used to compute the derivative of the CG function
     */
    Mdouble eps_;
    
};

///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::O> GaussO;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::X> GaussX;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::Y> GaussY;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::Z> GaussZ;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::YZ> GaussYZ;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::XZ> GaussXZ;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::XY> GaussXY;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::XYZ> GaussXYZ;
///Defines a short notation for the Gaussian CGFunction templated with a certain CGCoordinate.
typedef CGFunctions::Gauss<CGCoordinates::R> GaussR;
    
} //namespace CGFunctions
#include "Gauss.hcc"

#endif
