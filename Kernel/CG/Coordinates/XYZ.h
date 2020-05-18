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
#ifndef XYZ_H
#define XYZ_H

#include <GeneralDefine.h>
#include <iostream>
#include "Math/Vector.h"
#include <vector>
#include <array>
#include "BaseCoordinates.h"

class BaseParticle;

class BaseInteraction;

class DPMBase;

/*!
 * \namespace CGCoordinates
 * \brief The class in this namespace contain the position of a CGPoint, in the 
 * non-averaged directions, and functions that only depend on which non-averaged 
 * directions are used.
 * \details Contains base classes of CGFunctions, which in turn contains the 
 * base classes of CGPoint; CGPoint's and CGFunctions are always templated with 
 * one of these classes. The classes currently are: O, X, Y, Z, YZ, XZ, XY, XYZ. 
 * 
 * As the 1D classes X, Y, Z have some common functionality, they are derived 
 * from a common base class Base_X_Y_Z. Similarly, XY, XZ, YZ are derived from 
 * Base_XY_XZ_YZ.
 * 
 * \details See member CGCoordinates::XYZ for more details.
 */
namespace CGCoordinates
{

/*!
 * \brief Defines the position of the CGPoint, in the non-averaged directions, 
 * i.e. all directions on which spatial coarse-graining is applied (all 
 * directions for XYZ); all other directions are averaged over homogeneously.
 * \details In addition to defining the spatial variable, the classes in 
 * CGCoordinates contain all functions that only depend on how many coordinate 
 * directions are locally resolved (i.e. not averaged over); e.g. 
 * getDistanceSquared and getVolumeOfAveragedDimensions.
 * 
 * The CGCoordinates class should be chosen such that the system in homogeneous 
 * (symmetric) in the averaged directions. E.g., periodic chutes are homogeneous 
 * in x, y and t, but varying in z, so Z should be used.
 * 
 * See member CGCoordinates for more details.
 */
class XYZ : public BaseCoordinates
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
    void setXYZ(Vec3D p);
    
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
     * returns the number of variables (in this case three)
     */
    static const unsigned countVariables();

    static bool isResolvedIn(unsigned i) {return true;}

    static std::string getName();

protected:
    
    /*!
     * The position of the current CGPoint.
     */
    Vec3D p_;
};

/*!
 * \anchor spaceEvenly
 * \brief Creates a spatial mesh of CGPoints, i.e. the spatial positions at 
 * which the cg-variables are evaluated.
 * \details The mesh is created over the domain [min.X,max.X].[min.Y,max.Y].[min.Z,max.Z].
 * The mesh is defined by splitting the domain into nx.ny.nz smaller cubical 
 * subdomains; the centers of these cubes are the mesh points. This is done to 
 * keep the results comparable to a binning method.
 * \param[in] min the lower limit coordinate values of the meshed domain.
 * \param[in] max the upper limit coordinate values of the meshed domain.
 * \param[in] n the number of points in each spatial direction.
 * \param[out] points the vector of CGPoint's, which are now on a spatial mesh over the given domain.
 */
template<typename T>
typename std::enable_if<std::is_base_of<CGCoordinates::XYZ, typename T::CoordinatesType>::value, void>::type
spaceEvenly(Vec3D min, Vec3D max, std::vector<std::size_t> n, std::vector<T>& points)
{
    Vec3D delta = max - min;
    delta.X /= n[0];
    delta.Y /= n[1];
    delta.Z /= n[2];
    Vec3D start = min + 0.5 * delta;
    points.resize(n[0] * n[1] * n[2]);
    for (std::size_t i = 0; i < n[0]; i++)
    {
        for (std::size_t j = 0; j < n[1]; j++)
        {
            for (std::size_t k = 0; k < n[2]; k++)
            {
                points[(i * n[1] + j) * n[2] + k].coordinates.
                        setXYZ({start.X + delta.X * i,
                                start.Y + delta.Y * j,
                                start.Z + delta.Z * k});
            }
        }
    }
}
    
}
#endif
