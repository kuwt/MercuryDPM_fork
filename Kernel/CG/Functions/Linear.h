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
#ifndef Linear_H
#define Linear_H

#include "Polynomial.h"

namespace CGFunctions
{

/*!
 * \brief A specialisation of Polynomials for PolynomialType::Linear.
 * See Polynomial for details.
 */
template<class Coordinates>
class Linear : public Polynomial<Coordinates>
{
public:
    
    /*!
     * \brief Default constructor, simply sets the PolynomialType.
     * \details
     */
    Linear() : Polynomial<Coordinates>()
    {
        Polynomial<Coordinates>::setPolynomialType(PolynomialType::LINEAR);
    }
};

typedef CGFunctions::Linear<CGCoordinates::O> LinearO;
typedef CGFunctions::Linear<CGCoordinates::X> LinearX;
typedef CGFunctions::Linear<CGCoordinates::Y> LinearY;
typedef CGFunctions::Linear<CGCoordinates::Z> LinearZ;
typedef CGFunctions::Linear<CGCoordinates::YZ> LinearYZ;
typedef CGFunctions::Linear<CGCoordinates::XZ> LinearXZ;
typedef CGFunctions::Linear<CGCoordinates::XY> LinearXY;
typedef CGFunctions::Linear<CGCoordinates::XYZ> LinearXYZ;
    
} //namespace CGFunctions
#endif
