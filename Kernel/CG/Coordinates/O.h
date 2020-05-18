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
#ifndef O_H
#define O_H

#include <GeneralDefine.h>
#include <iostream>
#include "Math/Vector.h"
#include <vector>
#include <array>
#include "BaseCoordinates.h"

class BaseParticle;

class BaseInteraction;

class DPMBase;

namespace CGCoordinates
{

/*!
 * \brief Defines the non-averaged directions on which spatial coarse-graining 
 * is applied (none for O); all other directions (all for O) are averaged  
 * over homogeneously.
 * \details See XYZ for details.
 */
class O : public BaseCoordinates
{
public:
    
    /*!
     * \brief Writes the coordinate names in human-readable form to an ostream.
     */
    static void writeNames(std::ostream& os);
    
    /*!
     * \brief Writes the coordinates in human-readable form to an ostream.
     */
    void write(std::ostream& os) const;
    
    /*!
     * \brief returns the factor the CGFunction has to be divided by, due to
     * integrating the variables over the averaged dimensions, 1.0 for XYZ.
     */
    static Mdouble getVolumeOfAveragedDimensions(const Vec3D& min, const Vec3D& max);
    
    /*!
     * \brief Returns the square of the distance between the particle p and
     * the current CGPoint, in the non-averaged directions.
     */
    Mdouble getDistanceSquared(const Vec3D& p) const;
    
    /*!
     * \brief Normalises the coefficients of Polynomial CGFunction such that
     * the integral over all non-averaged dimensions is unity.
     */
    static void normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff);
    
    /*!
     * returns the number of variables (in this case three)
     */
    static const unsigned countVariables();
    
    static Mdouble getGaussPrefactor(Mdouble width, Mdouble cutoff)
    { return 1.0; }

    static bool isResolvedIn(unsigned dim) {return false;}
    
    static std::string getName();
    
};

/*!
 * See \ref spaceEvenly for details.
 */
template<typename T>
typename std::enable_if<std::is_base_of<CGCoordinates::O, typename T::CoordinatesType>::value, void>::type
spaceEvenly(Vec3D min, Vec3D max, std::vector<std::size_t> nAll, std::vector<T>& points)
{
    points.resize(1);
}
    
}
#endif
