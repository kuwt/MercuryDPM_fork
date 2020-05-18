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
#ifndef Base_X_Y_Z_H
#define Base_X_Y_Z_H

#include <GeneralDefine.h>
#include <iostream>
#include "Math/Vector.h"
#include <vector>
#include "BaseCoordinates.h"

class BaseParticle;

class BaseInteraction;

class DPMBase;

namespace CGCoordinates
{

/*!
 * \brief Contains common member functions of the X, Y, and Z classes.
 * \details As X, Y, and Z share a lot of functionality, the shared functions 
 * are stored in this common base class.
 */
class Base_X_Y_Z : public BaseCoordinates
{
public:
    
    /*!
     * \brief Computes the prefactor of the Gauss CGFunction, which is dependent
     * on the number of non-averaged dimensions.
     */
    static Mdouble getGaussPrefactor(Mdouble width, Mdouble cutoff);
    
    /*!
     * \brief Computes the prefactor of the Gauss line integral, which is dependent
     * on the number of non-averaged dimensions.
     */
    static Mdouble getGaussIntegralPrefactor(Mdouble distance, Mdouble width, Mdouble cutoff);
    
    /*!
     * \brief Normalises the coefficients of Polynomial CGFunction such that
     * the integral over all non-averaged dimensions is unity.
     */
    static void normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff);
    
    /*!
     * returns the number of variables (in this case one)
     */
    static const unsigned countVariables();
    
};
    
}
#endif
