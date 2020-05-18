/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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

//Note: This code is copied and adapted from hpGEM (see license above), version 22th of January 2016. It has been
//integrated into MercuryDPM at 16th of March 2017.

#ifndef MERCURY_SMALLVECTOR_H
#define MERCURY_SMALLVECTOR_H

#include "GeneralDefine.h"
#include "Logger.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <array>

template<unsigned int numberOfRows>
class SmallVector
{

public:
    
    SmallVector()
            : data_()
    {
    }
    
    SmallVector(const SmallVector& other)
            : data_(other.data_)
    {
    }
    
    SmallVector(SmallVector&& other)
            : data_(std::move(other.data_))
    {
    }
    
    SmallVector(const Mdouble array[])
            : data_()
    {
        std::copy(array, array + numberOfRows, data_.begin());
    }
    
    SmallVector(std::initializer_list<Mdouble> data)
            : data_()
    {
        logger.assert(data.size() == numberOfRows, "provided array has size %, but should have size %",
                      data.size(), numberOfRows);
        std::copy(data.begin(), data.end(), data_.begin());
    }
    
    SmallVector& operator=(const SmallVector& right)
    {
        std::copy(right.data_.begin(), right.data_.end(), data_.begin());
        return *this;
    }
    
    SmallVector& operator=(const std::array<Mdouble, numberOfRows> l)
    {
        std::copy(l.begin(), l.end(), data_.begin());
        return *this;
    }
    
    SmallVector operator+(const SmallVector& right) const
    {
        SmallVector result;
        std::transform(data_.begin(), data_.end(), right.data_.begin(), result.data_.begin(), std::plus<Mdouble>());
        return result;
    }
    
    SmallVector operator-(const SmallVector& right) const
    {
        SmallVector result;
        std::transform(data_.begin(), data_.end(), right.data_.begin(), result.data_.begin(), std::minus<Mdouble>());
        return result;
    }
    
    SmallVector operator*(const Mdouble& right) const
    {
        SmallVector result;
        std::transform(data_.begin(), data_.end(), result.data_.begin(), std::bind(std::multiplies<Mdouble>(),
                                                                                   std::placeholders::_1, right));
        return result;
    }
    
    ///Computes inner product between two vectors.
    Mdouble operator*(const SmallVector& right) const
    {
        return std::inner_product(right.data_.begin(), right.data_.end(), data_.begin(), 0.0);
    }
    
    SmallVector& operator/=(const Mdouble& right)
    {
        std::transform(data_.begin(), data_.end(), data_.begin(), std::bind(std::divides<Mdouble>(),
                                                                            std::placeholders::_1, right));
        return *this;
    }
    
    SmallVector operator/(const Mdouble& right) const
    {
        SmallVector result;
        std::transform(data_.begin(), data_.end(), result.data_.begin(), std::bind(std::divides<Mdouble>(),
                                                                                   std::placeholders::_1, right));
        return result;
    }
    
    void axpy(Mdouble a, const SmallVector& x)
    {
        for (unsigned int i = 0; i < numberOfRows; ++i)
        {
            data_[i] += a * x[i];
        }
    }
    
    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator==(const SmallVector& right) const
    {
        for (unsigned int i = 0; i < numberOfRows; ++i)
        {
            if (data_[i] != right[i])
            {
                return false;
            }
        }
        return true;
    }
    
    /// This function is dangerous to use, since it compares doubles without
    /// a tolerance interval to see if they are equal.
    bool operator<(const SmallVector& right) const
    {
        for (unsigned int i = 0; i < numberOfRows; ++i)
        {
            if (data_[i] < right[i])
            {
                return true;
            }
            if (data_[i] > right[i])
            {
                return false;
            }
        }
        return false;
    }
    
    SmallVector& operator+=(const SmallVector& right)
    {
        std::transform(data_.begin(), data_.end(), right.data_.begin(), data_.begin(), std::plus<Mdouble>());
        return *this;
    }
    
    SmallVector& operator-=(const SmallVector& right)
    {
        std::transform(data_.begin(), data_.end(), right.data_.begin(), data_.begin(), std::minus<Mdouble>());
        return *this;
    }
    
    SmallVector& operator*=(const double& right)
    {
        std::transform(data_.begin(), data_.end(), data_.begin(), std::bind(std::multiplies<Mdouble>(),
                                                                            std::placeholders::_1, right));
        return *this;
    }
    
    Mdouble& operator[](unsigned int n)
    {
        logger.assert(n < numberOfRows, "Requested entry %, but there are only % entries", n, numberOfRows);
        return data_[n];
    }
    
    const Mdouble& operator[](unsigned int n) const
    {
        logger.assert(n < numberOfRows, "Requested entry %, but there are only % entries", n, numberOfRows);
        return data_[n];
    }
    
    Mdouble& operator()(unsigned int n)
    {
        logger.assert(n < numberOfRows, "Requested entry %, but there are only % entries", n, numberOfRows);
        return data_[n];
    }
    
    const Mdouble& operator()(unsigned int n) const
    {
        logger.assert(n < numberOfRows, "Requested entry %, but there are only % entries", n, numberOfRows);
        return data_[n];
    }
    
    unsigned int size() const
    {
        return numberOfRows;
    }
    
    const Mdouble* data() const
    {
        return data_.data();
    }
    
    Mdouble* data()
    {
        return data_.data();
    }
    
    SmallVector operator-() const
    {
        return *this * -1.0;
    }
    
    Mdouble length() const
    {
        Mdouble sum = 0;
        for (Mdouble x : data_)
        {
            sum += x * x;
            logger(DEBUG, "x: %, sum: %", x, sum);
        }
        return std::sqrt(sum);
    }
    
    SmallVector<numberOfRows> getNormalised() const
    {
        return (*this) / length();
    }

private:
    std::array<Mdouble, numberOfRows> data_;
    
};

template<unsigned int numberOfRows>
SmallVector<numberOfRows> operator*(const Mdouble& left, const SmallVector<numberOfRows>& right)
{
    return right * left;
}

template<unsigned int numberOfRows>
std::ostream& operator<<(std::ostream& os, const SmallVector<numberOfRows>& A)
{
    os << "(";
    for (std::size_t i = 0; i < numberOfRows; ++i)
    {
        os << A[i] << " ";
    }
    os << ")";
    return os;
}

#endif //MERCURY_SMALLVECTOR_H
