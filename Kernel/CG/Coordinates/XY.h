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
#ifndef XY_H
#define XY_H

#include <GeneralDefine.h>
#include <iostream>
#include "Math/Vector.h"
#include "Base_XY_XZ_YZ.h"
#include <vector>
#include <array>

class BaseParticle;

class BaseInteraction;

class DPMBase;

namespace CGCoordinates
{

/*!
 * \brief Defines the non-averaged directions on which spatial coarse-graining 
 * is applied (the x- and y-direction for XY); all other directions are averaged  
 * over homogeneously.
 * \details See XYZ for details.
 */
class XY : public Base_XY_XZ_YZ
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
     * \brief Returns the length of the input vector in the non-averaged directions.
     */
    static Mdouble getLength(const Vec3D& p);
    
    /*!
     * \brief Returns the position of the current CGPoint, in the non-averaged
     * directions.
     */
    void setXY(Mdouble x, Mdouble y);
    
    /*!
     * \brief For the Interaction between particles/walls P and I, this function
     * returns the dot product between the normal vector of the interaction and
     * the branch vector from the current CGPoint towards I.
     */
    Mdouble getINormal(const BaseInteraction& c, const Vec3D& normal) const;
    
    /*!
     * \brief For the Interaction between particles/walls P and I, this function
     * returns the dot product between the normal vector of the interaction and
     * the branch vector from the current CGPoint towards P.
     */
    Mdouble getPNormal(const BaseInteraction& c, const Vec3D& normal) const;
    
    /*!
     * \brief For the Interaction between particles/walls P and I, this function
     * returns the dot product between the normal vector of the interaction and
     * the branch vector from the current CGPoint towards the contact point.
     */
    Mdouble getCNormal(const BaseInteraction& c, const Vec3D& normal) const;
    
    /*!
     * \brief For the Interaction between particles/walls P and I, this function
     * returns the square of the minimum distance between the the current
     * CGPoint and the branch vector between P and I.
     */
    Mdouble getTangentialSquared(const BaseInteraction& c, Mdouble pNormal) const;

    static bool isResolvedIn(unsigned i) {return i==2?false:true;}

    static std::string getName();

protected:
    
    /*!
     * The x-position of the current CGPoint.
     */
    Mdouble x_;
    /*!
     * The y-position of the current CGPoint.
     */
    Mdouble y_;
};

/*!
 * See \ref spaceEvenly for details.
 */
template<typename T>
typename std::enable_if<std::is_base_of<CGCoordinates::XY, typename T::CoordinatesType>::value, void>::type
spaceEvenly(Vec3D min, Vec3D max, std::vector<std::size_t> nAll, std::vector<T>& points)
{
    std::size_t n0 = nAll[0];
    std::size_t n1 = nAll[1];
    Mdouble delta0 = (max.X - min.X) / n0;
    Mdouble delta1 = (max.Y - min.Y) / n1;
    Mdouble start0 = min.X + 0.5 * delta0;
    Mdouble start1 = min.Y + 0.5 * delta1;
    points.resize(n0 * n1);
    for (std::size_t i = 0; i < n0; i++)
    {
        for (std::size_t j = 0; j < n1; j++)
        {
            points[i * n1 + j].coordinates.setXY(start0 + delta0 * i, start1 + delta1 * j);
        }
    }
}
    
}
#endif
