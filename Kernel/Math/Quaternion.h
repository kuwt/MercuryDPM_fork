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

/** Implementation of a 3D quaternion (by Vitaliy).
 */
///Modifications
/// 21:9:2009 - Added the inclusion guards and some include objects
/// \todo Need to generalize this to n-dimensional quaternions of any type
#ifndef QUATERNION_H
#define QUATERNION_H

#include <cmath>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "Math/Vector.h"
#include "Math/MatrixSymmetric.h"

#include "GeneralDefine.h"
#include "SmallVector.h"
#include "SmallMatrix.h"

/*!
 * \class Quaternion
 * \brief This class contains the 4 components of a quaternion and the standard operators and functions needed for quaternion arithmetic.
 * \details
 * A quaternion is a four-dimensional vector q = (q0,q1,q2,q3) that satisfies |q|^2=q0^2+q1^2+q2^2+q3^2=1.
 * It can be used to describe the orientation of an object in space.
 * 
 * The properties (e.g. inertia tensor) of any interactable (particle or wall) in MercuryDPM are described in its
 * reference frame. The quaternion is used to rotate the object into the lab frame (the system geometry).
 * A few examples to help visualising quaternions can be found
 * <a href="http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/examples/index.htm"> here </a>.
 * 
 * The unit quaternion, q=(1,0,0,0) denotes the state where the lab frame and the reference frame is identical.
 * 
 * To see how to convert a quaternion to Euler angles or to compare it to a rotation of an object around a axis,
 * see <a href="https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles">Wikipedia</a> for details.
 */
class Quaternion
{
public:
    
    /*!
     * \brief the zeroth component of the quaternion q = (q0,q1,q2,q3)
     */
    Mdouble q0;
    /*!
     * \brief the first component of the quaternion q = (q0,q1,q2,q3)
     */
    Mdouble q1;
    /*!
     * \brief the second component of the quaternion q = (q0,q1,q2,q3)
     */
    Mdouble q2;
    /*!
     * \brief the third component of the quaternion q = (q0,q1,q2,q3)
     */
    Mdouble q3;
    
    /*!
     * \brief Constructor; sets quaternion value to (1,0,0,0)
     */
    Quaternion();
    
    /*!
     * \brief Alternative constructor. Sets quaternion value to (q0,q1,q2,q3)
     */
    Quaternion(Mdouble q0, Mdouble q1, Mdouble q2, Mdouble q3);
    
    /**
     * \todo should be explicit as teh conversion is not unique
     * @param normal
     */
    Quaternion(Vec3D normal)
    {
        setOrientationViaNormal(normal);
    }
    
    /*!
     * \brief Sets quaternion value to (1,0,0,0)
     */
    void setUnity();
    
    /*!
     * \brief Checks if the quaternion value is (1,0,0,0) 
     * \details Checks if ALL elements are zero
     * \return          TRUE if q0 equals one and ALL other elements are zero
     * \bug use isEqual instead of ==
    */
    bool isUnity() const
    {
        return q0 == 1.0 && q1 == 0.0 && q2 == 0.0 && q3 == 0.0;
    };
    
    /*!
     * \brief Adds another quaternion and returns the result.
     */
    Quaternion operator+(const Quaternion& a) const;
    
    /*!
     * \brief Subtracts another quaternion and returns the result.
     */
    Quaternion operator-(const Quaternion& a) const;
    
    /*!
     * \brief Multiplies by a scalar
     */
    Quaternion operator*(Mdouble a) const;
    
    /*!
     * \brief Divides by a scalar.
     */
    Quaternion operator/(Mdouble a) const;
    
    /*!
     * \brief Adds another quaternion
     */
    Quaternion& operator+=(const Quaternion& a);
    
    /*!
     * \brief Subtracts another quaternion.
     */
    Quaternion& operator-=(const Quaternion& a);
    
    /*!
     * \brief Multiplies *this by a scalar
     */
    Quaternion& operator*=(Mdouble a);
    
    /*!
     * \brief Divides by a scalar
     */
    Quaternion& operator/=(Mdouble a);
    
    /*!
     * \brief Makes this Quaternion unit length |q|=1
     */
    void normalise();
    
    /*!
     * \brief Make this Quaternion a certain length  |q|=length
     */
    void setLength(Mdouble length);
    
    /*!
     * \brief Calculates the distance between two Quaternion: \f$ \sqrt{\left(a-b\right) \cdot \left(a-b\right)} \f$
     */
    static Mdouble getDistance(const Quaternion& a, const Quaternion& b);
    
    /*!
     * \brief Calculates the squared distance between two Quaternion: \f$ \left(a-b\right) \cdot \left(a-b\right) \f$
     */
    static Mdouble getDistanceSquared(const Quaternion& a, const Quaternion& b);
    
    /*!
     * \brief Calculates the length of a Quaternion: \f$ \sqrt{a\cdot a} \f$
     */
    static Mdouble getLength(const Quaternion& a);
    
    /*!
     * \brief Calculates the squared length of a Quaternion: \f$ a\cdot a \f$
     */
    static Mdouble getLengthSquared(const Quaternion& a);
    
    /*!
     * \brief Calculates the length of this Quaternion: \f$ \sqrt{a\cdot a} \f$
     */
    Mdouble getLength() const;
    
    /*!
     * \brief Calculates the squared length of this Quaternion: \f$ a\cdot a \f$
     */
    Mdouble getLengthSquared() const;
    
    /*!
     * \brief Returns the requested component of this Quaternion
     */
    Mdouble getComponent(int index) const;
    
    ///\todo should the arguments be passed by reference?
    Quaternion angularVelocityBodyFixedFrameToAngularDisplacement(Vec3D v) const;
    
    /*!
     * \brief Converts an angular momentum v=omega into a 
     * quaternion rate of change, q(t+dt)-q(t)/dt
     */
    Quaternion angularDisplacementTimeDerivative(Vec3D v) const;
    
    void updateAngularDisplacement(Vec3D angularVelocityDt);
    
    /*!
     * \brief Converts quaternion rate of change into an angular momentum omega
     */
    Vec3D applyCInverse(Quaternion q) const;
    
    /*!
     * \brief Sets the requested component of this Quaternion to the requested value
     */
    void setComponent(int index, double val);
    
    /*!
     * \brief Checks if the length this Quaternion is equal the length of other with a certain tolerance
     */
    bool isEqualTo(const Quaternion& other, double tol) const;
    
    /*!
     * \brief Returns a unit Quaternion based on a.
     */
    static Quaternion getUnitQuaternion(const Quaternion& a);
    
    /*!
     * \brief Adds elements to an output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Quaternion& a);
    
    /*!
     * \brief Adds elements to an input stream
     */
    friend std::istream& operator>>(std::istream& is, Quaternion& a);
    
    /*!
     * \brief Adds a scalar to a quaternion
     */
    friend Quaternion operator+(Mdouble a, const Quaternion& b);
    
    /*!
     * \brief Subtracts the elements of a quaternion from a scalar
     */
    friend Quaternion operator-(Mdouble a, const Quaternion& b);
    
    /*!
     * \brief Subtracts a quaternion 
     */
    friend Quaternion operator-(const Quaternion& a);
    
    /*!
     * \brief Multiplies all elements by a scalar
     */
    friend Quaternion operator*(Mdouble a, const Quaternion& b);
    
    /*!
     * \brief Convert a quaternion to Euler angles.
     *  See <a href="https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles">Wikipedia</a> for details.
     */
    Vec3D getEuler() const;
    
    
    /*!
     * \brief Convert Euler angles to a quaternion.
     *  See <a href="https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles">Wikipedia</a> for details.
     */
    void setEuler(const Vec3D& e);
    
    /*!
     * \brief Converts a quaternion to the rotation angle in the XY plane (for Mercury2D).
     *  See <a href="https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles">Wikipedia</a> for details.
     */
    Mdouble getAngleZ() const;
    
    /*!
     * \brief Converts the rotation angle in the XY plane into a quaternion (for Mercury2D).
     *  See <a href="https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles">Wikipedia</a> for details.
     */
    void setAngleZ(Mdouble psi);
    
    /*!
     * \brief Converts the inverse inertia tensor from the reference frame to the lab frame; see 
     *  See <a href="QuaternionsWouter.pdf">QuaternionsWouter.pdf</a> for details, where this operation is denoted by \f$A*I^{-1}*A^T\f$.
     */
    MatrixSymmetric3D rotateInverseInertiaTensor(const MatrixSymmetric3D& invI) const;
    
    void rotateTensor(SmallMatrix<3, 3> I) const;
    
    /*!
     * \brief Converts the quaternions into a normal vector by rotating the vector x=(1,0,0); see
     *  See <a href="https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix">Wiki</a> for details.
     */
    Vec3D getAxis() const;
    
    /**
     * Used to the the normal of an InfiniteWall that has a normal into the x-direction by default.
     * This can be changed by resetting the wall orientation; thus the normal is the vector (1,0,0) rotated by this quaternion.
     */
    void setOrientationViaNormal(Vec3D normal);
    
    /**
     * Retrieves the rotation matrix
     * \param[out] A The rotation matrix
     */
    void getRotationMatrix(SmallMatrix<3, 3>& A) const; //note: A is output parameter
    
    /**
     * Applies the rotation to a position
     */
    void rotate(Vec3D& position) const;
    
    /**
     * Applies the rotation to a position
     */
    void rotate(SmallVector<3>& position) const;
    
    /**
     * Applies the inverse rotation to a position
     */
    void rotateBack(Vec3D& position) const;
    
    /**
     * Calculates the distance from a wall through p0 whose normal is the vector (1,0,0). Used for the calculation of the distance in InfiniteWalls
     */
    Mdouble getDistance(Vec3D p, Vec3D p0) const;
};

#endif
