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
#ifndef R_H
#define R_H

#include <GeneralDefine.h>
#include <iostream>
#include "Math/Vector.h"
#include "Base_X_Y_Z.h"
#include <vector>
#include <array>
#include "Math/ExtendedMath.h"

class BaseParticle;

class BaseInteraction;

class DPMBase;

namespace CGCoordinates
{

/*!
 * \brief Defines the non-averaged directions on which spatial coarse-graining 
 * is applied (the x-direction for R); all other directions are averaged  
 * over homogeneously.
 * \details See XYZ for details.
 */
class R : public Base_X_Y_Z
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
    void setR(Mdouble r);
    
    const Mdouble getR() const
    { return r_; }
    
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
     * \brief Normalises the coefficients of Polynomial CGFunction such that
     * the integral over all non-averaged dimensions is unity.
     */
    static void normalisePolynomialCoefficients(std::vector<Mdouble>& coefficients, Mdouble cutoff);
    
    /*!
     * returns the weight function needed in getTotal (required as functions cannot be )
     * @return
     */
    Mdouble getWeight() override;
    
    static Mdouble getDomainVolume(const Vec3D& min, const Vec3D& max);

    static bool isResolvedIn(unsigned i) {return i==0?true:false;}

    static std::string getName();

protected:
    
    /*!
     * The x-position of the current CGPoint.
     */
    Mdouble r_;
};

/*!
 * See \ref spaceEvenly for details.
 * \todo limit to r>0
 */
template<typename T>
typename std::enable_if<std::is_base_of<CGCoordinates::R, typename T::CoordinatesType>::value, void>::type
spaceEvenly(Vec3D min, Vec3D max, std::vector<std::size_t> nAll, std::vector<T>& points)
{
    //note, the x-coordinate represents the r-coordinate here
    std::size_t n = nAll[0];
    Mdouble start = min.X;
    Mdouble delta = (max.X - start) / n;
//    Mdouble end = sqrt(mathsFunc::square(std::fmax(-min.X,max.X))+mathsFunc::square(std::fmax(-min.Y,max.Y)));
//    Mdouble start = 0;
//    Mdouble delta = (end - start) / n;
    if (start < 0.0)
    {
        if (max.X > 0.0)
        {
            start = 0.0;
        }
        else
        {
            start = -max.X;
            delta = (-min.X - start) / n;
        }
        std::cerr << "Warning: r should be positive, but x-limits are not;"
                     "grid is set within " << start << "<r<" << start + delta * (n + 1) << std::endl;
    }
    start += 0.5 * delta;
    points.resize(n);
    for (std::size_t i = 0; i < n; i++)
    {
        points[i].coordinates.setR(start + delta * i);
    }
}
    
}
#endif
