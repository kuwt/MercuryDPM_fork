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

#ifndef MERCURY_SMALLMATRIX_H
#define MERCURY_SMALLMATRIX_H

#include "GeneralDefine.h"
#include "Logger.h"
#include "SmallVector.h"

/// \class SmallMatrix
/// \brief Data type for small dense matrix.
///
/// \details Stores small dense matrix efficiently.
/// It only store Mdoubles as this is the main type linear algebra is done on in hpGEM
/// It stores the matrix in fortran style (column-major) to give quicker access to extern BLAS libraries.
/// For example,  the order they are stored in a 2x2 matrix is
///     0   2
///     1   3
template<unsigned int numberOfRows, unsigned int numberOfColumns>
class SmallMatrix
{
public:
    
    /// \brief Constructs a matrix of size n-rows by m-columns.
    SmallMatrix()
            : data_()
    {
    }
    
    SmallMatrix(const SmallVector<numberOfRows>& other)
            : data_()
    {
        logger.assert(numberOfColumns == 1, "Trying to construct a matrix with more than 1 columns from a vector");
        std::copy(other.data(), other.data() + numberOfRows, data_.begin());
    }
    
    /// \brief Constructs a matrix of size n-rows by m-columns and initialises all entry to a constant
    SmallMatrix(const Mdouble& c)
            : data_()
    {
        data_.fill(c);
    }
    
    SmallMatrix(const std::initializer_list<SmallVector<numberOfRows>>& entries)
            : data_()
    {
        logger.assert(entries.size() == numberOfColumns, "expected a matrix with % "
                                                         "columns, but got a matrix with % columns", numberOfColumns,
                      entries.size());
        unsigned int column = 0;
        for (const SmallVector<numberOfRows>& entry : entries)
        {
            for (unsigned int i = 0; i < numberOfRows; ++i)
            {
                (*this)(i, column) = entry[i];
            }
            ++column;
        }
    }
    
    /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B and A are both matrices
    SmallMatrix(const SmallMatrix& other)
            : data_()
    {
        std::copy(other.data_.begin(), other.data_.end(), data_.begin());
    }
    
    /// \brief Glues one or more vectors with the same number of rows together
    SmallMatrix(std::array<SmallVector<numberOfRows>, numberOfColumns> entries)
            : data_()
    {
        for (unsigned int i = 0; i < numberOfRows; ++i)
        {
            for (unsigned int j = 0; j < numberOfColumns; ++j)
            {
                (*this)(i, j) = entries[j][i];
            }
        }
    }
    
    /// \brief Move Matrix from another Matrix
    SmallMatrix(SmallMatrix&& other)
            : data_(std::move(other.data_))
    {
    }
    
    /// \brief defines the operator(n,m) to access the element on row n and column m
    Mdouble& operator()(unsigned int n, unsigned int m)
    {
        logger.assert(n < numberOfRows, "Requested row number % for a matrix with only % rows", n, numberOfRows);
        logger.assert(m < numberOfColumns, "Requested column number % for a matrix with only % columns", m,
                      numberOfColumns);
        return data_[n + m * numberOfRows];
    }
    
    /// \brief defines the operator(n,m) to access the element on row n and column m
    const Mdouble& operator()(unsigned int n, unsigned int m) const
    {
        logger.assert(n < numberOfRows, "Requested row number % for a matrix with only % rows", n, numberOfRows);
        logger.assert(m < numberOfColumns, "Requested column number % for a matrix with only % columns", m,
                      numberOfColumns);
        return data_[n + m * numberOfRows];
    }
    
    /// \brief Access the n linear element in the matrix.
    Mdouble& operator[](const unsigned int n)
    {
        logger.assert(n < numberOfRows * numberOfColumns, "Requested entry % for a matrix with only % entries", n,
                      numberOfRows * numberOfColumns);
        return data_[n];
    }
    
    const Mdouble& operator[](const unsigned int n) const
    {
        logger.assert(n < numberOfRows * numberOfColumns, "Requested entry % for a matrix with only % entries", n,
                      numberOfRows * numberOfColumns);
        return data_[n];
    }
    
    /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j= A_ij B_,j
    SmallVector<numberOfRows> operator*(SmallVector<numberOfColumns>& right);
    
    SmallVector<numberOfRows> operator*(SmallVector<numberOfColumns>& right) const;
    
    /// \brief Does matrix A_ij=scalar*B_ij
    SmallMatrix operator*(const Mdouble& right) const
    {
        SmallMatrix result;
        std::transform(data_.begin(), data_.end(), result.data_.begin(),
                       std::bind(std::multiplies<Mdouble>(), std::placeholders::_1, right));
        return result;
    }
    
    /// \brief Does matrix A_ij = B_ik * C_kj
    template<unsigned int K>
    SmallMatrix<numberOfRows, K> operator*(const SmallMatrix<numberOfColumns, K>& other);
    
    template<unsigned int K>
    SmallMatrix<numberOfRows, K> operator*(const SmallMatrix<numberOfColumns, K>& other) const;
    
    SmallMatrix& operator+=(const SmallMatrix& other)
    {
        std::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::plus<Mdouble>());
        return *this;
    }
    
    SmallMatrix& operator-=(const SmallMatrix& other)
    {
        std::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::minus<Mdouble>());
        return *this;
    }
    
    SmallMatrix operator+(const SmallMatrix& other) const
    {
        SmallMatrix result;
        std::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::plus<Mdouble>());
        return result;
    }
    
    SmallMatrix operator-(const SmallMatrix& other) const
    {
        SmallMatrix result;
        std::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::minus<Mdouble>());
        return result;
    }
    
    SmallMatrix operator-() const
    {
        return *this * -1.;
    }
    
    /// \brief Does matrix A_ij=scalar*A_ij
    SmallMatrix& operator*=(const Mdouble& scalar)
    {
        std::transform(data_.begin(), data_.end(), data_.begin(),
                       std::bind(std::multiplies<Mdouble>(), std::placeholders::_1, scalar));
        return *this;
    }
    
    /// \brief Does matrix A_ij = A_ik * B_kj
    /// note that other must be square because this is a fixed-size matrix
    SmallMatrix& operator*=(const SmallMatrix<numberOfColumns, numberOfColumns>& other);
    
    /// \brief Does matrix A_ij=scalar*A_ij
    SmallMatrix& operator/=(const Mdouble& scalar)
    {
        std::transform(data_.begin(), data_.end(), data_.begin(),
                       std::bind(std::divides<Mdouble>(), std::placeholders::_1, scalar));
        return *this;
    }
    
    /// \brief this does element by divided by a scalar
    SmallMatrix operator/(const Mdouble& scalar) const
    {
        SmallMatrix result;
        std::transform(data_.begin(), data_.end(), result.data_.begin(),
                       std::bind(std::divides<Mdouble>(), std::placeholders::_1, scalar));
        return result;
    }
    
    /// \brief Assigns one matrix to another.
    SmallMatrix& operator=(const SmallMatrix& right)
    {
        std::copy(right.data_.begin(), right.data_.end(), data_.begin());
        return *this;
    }
    
    /// \brief Assigns one matrix to another.
    SmallMatrix& operator=(SmallMatrix&& right)
    {
        std::move(right.data_.begin(), right.data_.end(), data_.begin());
        return *this;
    }
    
    /// \brief computeWedgeStuffVector.
    SmallVector<numberOfRows> computeWedgeStuffVector() const;
    
    /// \brief Applies the matrix y=ax + y, where x is another matrix and a is a scalar
    void axpy(Mdouble a, const SmallMatrix& x)
    {
        for (unsigned int i = 0; i < numberOfRows * numberOfColumns; ++i)
        {
            data_[i] += a * x[i];
        }
    }
    
    /// \brief Get total number of Matrix entries
    unsigned int size() const
    {
        return numberOfRows * numberOfColumns;
    }
    
    /// \brief Get the number of rows
    unsigned int getNumberOfRows() const
    {
        return numberOfRows;
    }
    
    ///\deprecated Does not conform naming convention, please use getNumberOfRows instead.
    unsigned int getNRows() const
    {
        return getNumberOfRows();
    }
    
    /// \brief Get the number of columns
    unsigned int getNumberOfColumns() const
    {
        return numberOfColumns;
    }
    
    unsigned int getNCols() const
    {
        return getNumberOfColumns();
    }
    
    /// \brief get the j^th column
    SmallVector<numberOfRows> getColumn(unsigned int j) const
    {
        logger.assert(j < numberOfColumns, "Asked for column %, but there are only % columns", j, numberOfColumns);
        return SmallVector<numberOfRows>(data() + j * numberOfRows);
    }
    
    /// \brief get the i^th row
    SmallVector<numberOfColumns> getRow(unsigned int i) const
    {
        logger.assert(i < numberOfRows, "Asked for row %, but there are only % rows", i, numberOfRows);
        SmallVector<numberOfColumns> result;
        for (unsigned int j = 0; j < numberOfColumns; ++j)
        {
            result[j] = (*this)(i, j);
        }
        return result;
    }
    
    /// \brief Return the LUfactorisation of the matrix
    SmallMatrix LUfactorisation() const;
    
    Mdouble determinant() const;
    
    /// \brief return the inverse in the vector result. The size of result matches the matrix.
    SmallMatrix inverse() const;
    
    SmallMatrix<numberOfColumns, numberOfRows> transpose() const
    {
        SmallMatrix<numberOfColumns, numberOfRows> result;
        for (unsigned int i = 0; i < numberOfRows; ++i)
        {
            for (unsigned int j = 0; j < numberOfColumns; ++j)
            {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }
    
    /// \brief solves Ax=B where A is the current matrix and B is passed in. The result is returned in B.
    template<unsigned int numberOfRightHandSideColumns>
    void solve(SmallMatrix<numberOfRows, numberOfRightHandSideColumns>& B) const;
    
    /// \brief solves Ax=b where A is the current matrix and NumericalVector b
    /// is the input parameter. The result is returned in b.
    void solve(SmallVector<numberOfRows>& b) const;
    
    Mdouble* data()
    {
        return data_.data();
    }
    
    const Mdouble* data() const
    {
        return data_.data();
    }

private:
    /// The actually data of the matrix class
    std::array<Mdouble, numberOfRows * numberOfColumns> data_;
};

/// Writes nicely formatted entries of the Matrix A to the stream os.
template<unsigned int numberOfRows, unsigned int numberOfColumns>
std::ostream& operator<<(std::ostream& os, const SmallMatrix<numberOfRows, numberOfColumns>& A)
{
    for (unsigned int i = 0; i < numberOfRows; ++i)
    {
        os << A.getRow(i) << std::endl;
    }
    return os;
}

///Multiplies a matrix with a Mdouble
template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns>
operator*(const Mdouble d, const SmallMatrix<numberOfRows, numberOfColumns>& mat)
{
    return mat * d;
}

///Multiplies a matrix with a vector
template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallVector<numberOfColumns> operator*(SmallVector<numberOfRows>& vec, SmallMatrix<numberOfRows, numberOfColumns>& mat);

#include "SmallMatrix_impl.h"


#endif //MERCURY_SMALLMATRIX_H
