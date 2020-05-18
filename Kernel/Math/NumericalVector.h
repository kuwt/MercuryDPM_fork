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

//Note: This code is copied and adapted from hpGEM (see license above), version 17th of September 2015. It has been
//integrated into MercuryDPM at 17th of September 2015.


#ifndef NumericalVector_H_
#define NumericalVector_H_

#include "Logger.h"
#include <vector>

/// \class NumericalVector<T>
/// \brief This is a vector of doubles
///
/// \details
/// This implements a vector of doubles and all the standard operators for it.
template<typename T = Mdouble>
class NumericalVector
{
public:
    
    NumericalVector<T>()
            : data_()
    {
    }
    
    explicit NumericalVector<T>(std::size_t m)
            : data_(m)
    {
    }
    
    NumericalVector<T>(std::initializer_list<T> l)
            : data_(l)
    {
    }
    
    NumericalVector<T>(const NumericalVector& other)
            : data_(other.data_)
    {
    }
    
    NumericalVector<T>(NumericalVector&& other)
            : data_(std::move(other.data_))
    {
    }
    
    NumericalVector<T>(const T array[], std::size_t size)
            : data_(array, array + size)
    {
    }
    
    void resize(std::size_t size)
    {
        if (size != data_.size())
        {
            data_.resize(size);
        }
    }
    
    NumericalVector<T>& operator=(const NumericalVector<T>& right)
    {
        data_ = right.data_;
        return *this;
    }
    
    NumericalVector<T>& operator=(const std::initializer_list<T> l)
    {
        data_ = l;
        return *this;
    }
    
    NumericalVector<T> operator+(const NumericalVector<T>& right) const
    {
        NumericalVector<T> result(*this);
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
        for (std::size_t i = 0; i < data_.size(); i++)
            result.data_[i] += right.data_[i];
        return result;
    }
    
    NumericalVector<T> operator-(const NumericalVector<T>& right) const
    {
        NumericalVector<T> result(*this);
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
        for (std::size_t i = 0; i < data_.size(); i++)
            result.data_[i] -= right.data_[i];
        return result;
    }
    
    NumericalVector<T> operator*(const T& right) const
    {
        NumericalVector<T> result(*this);
        for (T& d : result.data_)
            d *= right;
        return result;
    }
    
    
    T operator*(const NumericalVector& right) const
    {
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have equal length.");
        T sum = 0;
        for (std::size_t i = 0; i < data_.size(); i++)
            sum += data_[i] * right.data_[i];
        return sum;
    }
    
    NumericalVector<T>& operator/=(const T& right)
    {
        for (T& d : data_)
            d /= right;
        
        return *this;
    }
    
    NumericalVector<T> operator/(const T& right) const
    {
        NumericalVector<T> result(*this);
        return (result /= right);
        
    }
    
    NumericalVector<T>& operator+=(const NumericalVector<T>& right)
    {
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
        for (std::size_t i = 0; i < data_.size(); i++)
            data_[i] += right.data_[i];
        return *this;
    }
    
    
    NumericalVector<T>& operator-=(const NumericalVector<T>& right)
    {
        logger.assert(data_.size() == right.data_.size(), "Vectors don't have the same size");
        for (std::size_t i = 0; i < data_.size(); i++)
            data_[i] -= right.data_[i];
        return *this;
    }
    
    NumericalVector<T>& operator*=(const T& right)
    {
        for (T& d : data_)
            d *= right;
        return *this;
    }
    
    T& operator[](std::size_t n)
    {
        logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
        return data_[n];
    }
    
    const T& operator[](std::size_t n) const
    {
        logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
        return data_[n];
    }
    
    T& operator()(std::size_t n)
    {
        logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
        return data_[n];
    }
    
    const T& operator()(std::size_t n) const
    {
        logger.assert(n < data_.size(), "Requested entry %, but there are only % entries", n, data_.size());
        return data_[n];
    }
    
    std::size_t size() const
    {
        return data_.size();
    }
    
    const T* data() const
    {
        return data_.data();
    }
    
    T* data()
    {
        return data_.data();
    }

private:
    std::vector<T> data_;
    
    
};

template<typename T = Mdouble>
NumericalVector<T> operator*(const T& left, const NumericalVector<T>& right);

template<typename T = Mdouble>
NumericalVector<T> operator-(const NumericalVector<T>& right);

template<typename T = Mdouble>
std::ostream& operator<<(std::ostream& os, const NumericalVector<T>& A);


#endif /* NumericalVector_H_ */


