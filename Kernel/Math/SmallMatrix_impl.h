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

//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


#include "SmallMatrix.h"


extern "C"
{

/// This does matrix times vector and is from blas level 2.
void dgemv_(const char* trans, int* m, int* n, double* alpha, double* A, int* LDA, double* x, int* incx, double* beta,
            double* y, int* incy);

///This is the gernal matrix multiplication from blas level 3
int
dgemm_(const char* transA, const char* transB, int* M, int* N, int* k, double* alpha, double* A, int* LDA, double* B,
       int* LDB, double* beta, double* C, int* LDC);

///This is the gerneral scalar times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
int daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);

/// This is LU factorisation of the matrix A. This has been taken from LAPACK
void dgetrf_(int* M, int* N, double* A, int* lda, int* IPIV, int* INFO);

/// This is the inverse calulation also from LAPACK. Calculates inverse if you pass it the LU factorisation.
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

/// This is used for solve Ax=B for x. Again this is from LAPACK.
void dgesv_(int* N, int* NRHS, double* A, int* lda, int* IPIV, double* B, int* LDB, int* INFO);
}


template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallVector<numberOfRows> SmallMatrix<numberOfRows, numberOfColumns>::operator*(SmallVector<numberOfColumns>& right)
{
    if (numberOfRows == 0)
    {
        logger(WARN, "Trying to multiply a vector with a matrix without any rows.");
        return SmallVector<numberOfRows>();
    }
    if (numberOfColumns == 0)
    {
        logger(WARN, "Trying to multiply a vector with a matrix without any columns.");
        return SmallVector<numberOfRows>();
    }
    int nr = numberOfRows;
    int nc = numberOfColumns;
    
    int i_one = 1;
    double d_one = 1.0;
    double d_zero = 0.0;
    
    SmallVector<numberOfRows> result;
    
    logger(DEBUG, "Matrix size: % x % \n Vector size: %", nr, nc, right.size());
    
    dgemv_("N", &nr, &nc, &d_one, this->data(), &nr, right.data(), &i_one, &d_zero, result.data(), &i_one);
    return result;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallVector<numberOfRows>
SmallMatrix<numberOfRows, numberOfColumns>::operator*(SmallVector<numberOfColumns>& right) const
{
    if (numberOfRows == 0)
    {
        logger(WARN, "Trying to multiply a vector with a matrix without any rows.");
        return SmallVector<numberOfRows>();
    }
    if (numberOfColumns == 0)
    {
        logger(WARN, "Trying to multiply a vector with a matrix without any columns.");
        return SmallVector<numberOfRows>();
    }
    int nr = numberOfRows;
    int nc = numberOfColumns;
    
    int i_one = 1;
    double d_one = 1.0;
    double d_zero = 0.0;
    
    SmallVector<numberOfRows> result;
    
    logger(DEBUG, "Matrix size: % x % \n Vector size: %", nr, nc, right.size());
    
    dgemv_("N", &nr, &nc, &d_one, (const_cast<double*>(this->data())), &nr, right.data(), &i_one, &d_zero,
           result.data(), &i_one);
    return result;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
template<unsigned int K>
SmallMatrix<numberOfRows, K>
SmallMatrix<numberOfRows, numberOfColumns>::operator*(const SmallMatrix<numberOfColumns, K>& other)
{
    int i = numberOfRows;
    int j = numberOfColumns;
    int k = K;
    
    if (numberOfColumns == 0)
    {
        logger(WARN, "Trying to multiply a matrix with a matrix without any columns.");
        return SmallMatrix<numberOfRows, K>();
    }
    //The result of the matrix is left.numberOfRows, right.numberOfColumns()
    SmallMatrix<numberOfRows, K> C;
    
    double d_one = 1.0;
    double d_zero = 0.0;
    
    //Let the actual multiplication be done by Fortran
    dgemm_("N", "N", &i, &k, &j, &d_one, this->data(), &i, const_cast<double*>(other.data()), &j, &d_zero, C.data(),
           &i);
    
    return C;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
template<unsigned int K>
SmallMatrix<numberOfRows, K>
SmallMatrix<numberOfRows, numberOfColumns>::operator*(const SmallMatrix<numberOfColumns, K>& other) const
{
    int i = numberOfRows;
    int j = numberOfColumns;
    int k = K;
    
    if (numberOfColumns == 0)
    {
        logger(WARN, "Trying to multiply a matrix with a matrix without any columns.");
        return SmallMatrix<numberOfRows, K>();
    }
    //The result of the matrix is left.Nrows, right.NCols()
    SmallMatrix<numberOfRows, K> C;
    
    double d_one = 1.0;
    double d_zero = 0.0;
    
    //Let the actual multiplication be done by Fortran
    dgemm_("N", "N", &i, &k, &j, &d_one, const_cast<double*>(this->data()), &i, const_cast<double*>(other.data()), &j,
           &d_zero, C.data(), &i);
    
    return C;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns>&
SmallMatrix<numberOfRows, numberOfColumns>::operator*=(const SmallMatrix<numberOfColumns, numberOfColumns>& other)
{
    //blas does not support in-place multiply
    return (*this) = (*this) * other;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallVector<numberOfRows> SmallMatrix<numberOfRows, numberOfColumns>::computeWedgeStuffVector() const
{
    //copied from MiddleSizeMatrix to prevent constructing a temporary MiddleSizeMatrix
    logger.assert_debug(numberOfColumns == numberOfRows - 1,
                  "Matrix has wrong dimensions to construct the wedge stuff vector");
    SmallVector<numberOfRows> result;
    
    switch (numberOfRows)
    {
        case 2:
            result[0] = -(*this)(1, 0);
            result[1] = +(*this)(0, 0);
            break;
        case 3:
            result[0] = (*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1);
            result[1] = (*this)(0, 1) * (*this)(2, 0) - (*this)(0, 0) * (*this)(2, 1); // includes minus sign already!
            result[2] = (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1);
            break;
        case 4:
            result[0] = (*this)(1, 0) * (-(*this)(2, 1) * (*this)(3, 2) + (*this)(3, 1) * (*this)(2, 2)) +
                        (*this)(2, 0) * ((*this)(1, 1) * (*this)(3, 2) - (*this)(3, 1) * (*this)(1, 2)) +
                        (*this)(3, 0) * (-(*this)(1, 1) * (*this)(2, 2) + (*this)(2, 1) * (*this)(1, 2));
            
            result[1] = (*this)(0, 0) * ((*this)(2, 1) * (*this)(3, 2) - (*this)(3, 1) * (*this)(2, 2)) +
                        (*this)(2, 0) * (-(*this)(0, 1) * (*this)(3, 2) + (*this)(3, 1) * (*this)(0, 2)) +
                        (*this)(3, 0) * ((*this)(0, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(0, 2));
            result[2] = (*this)(0, 0) * (-(*this)(1, 1) * (*this)(3, 2) + (*this)(3, 1) * (*this)(1, 2)) +
                        (*this)(1, 0) * ((*this)(0, 1) * (*this)(3, 2) - (*this)(3, 1) * (*this)(0, 2)) +
                        (*this)(3, 0) * (-(*this)(0, 1) * (*this)(1, 2) + (*this)(1, 1) * (*this)(0, 2));
            result[3] = (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) +
                        (*this)(1, 0) * (-(*this)(0, 1) * (*this)(2, 2) + (*this)(2, 1) * (*this)(0, 2)) +
                        (*this)(2, 0) * ((*this)(0, 1) * (*this)(1, 2) - (*this)(1, 1) * (*this)(0, 2));
            break;
        default:
            logger(ERROR, "Wedge product not implemented for this dimension");
    } //end switch
    
    return (result);
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns> SmallMatrix<numberOfRows, numberOfColumns>::LUfactorisation() const
{
    int nr = numberOfRows;
    int nc = numberOfColumns;
    int nPivot = std::min(numberOfRows, numberOfColumns);
    int iPivot[nPivot];
    
    SmallMatrix result(*this);
    
    int info;
    
    dgetrf_(&nr, &nc, result.data(), &nr, iPivot, &info);
    
    return result;
}

//class template specialization for this one function is a waste of code duplication
//just let the compiler figure out which case it needs
template<unsigned int numberOfRows, unsigned int numberOfColumns>
Mdouble SmallMatrix<numberOfRows, numberOfColumns>::determinant() const
{
    logger.assert_debug(numberOfRows == numberOfColumns, "Matrix should be square to have a determinant!");
    
    switch (numberOfRows)
    {
        case 0:
            return 1;
        case 1:
            return (*this)(0, 0);
        case 2:
            return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
        
        case 3:
            return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1)) -
                   (*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) - (*this)(2, 0) * (*this)(1, 2)) +
                   (*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1));
        
        case 4:
            return ((*this)(3, 0) * (*this)(2, 1) * (*this)(0, 3) - (*this)(2, 0) * (*this)(3, 1) * (*this)(0, 3)) *
                   (*this)(1, 2) +
                   (-(*this)(3, 0) * (*this)(0, 3) * (*this)(2, 2) + (*this)(2, 0) * (*this)(0, 3) * (*this)(3, 2)) *
                   (*this)(1, 1) +
                   ((*this)(3, 1) * (*this)(0, 3) * (*this)(2, 2) - (*this)(2, 1) * (*this)(0, 3) * (*this)(3, 2)) *
                   (*this)(1, 0) +
                   (-(*this)(3, 0) * (*this)(2, 1) * (*this)(1, 3) + (*this)(2, 0) * (*this)(3, 1) * (*this)(1, 3) +
                    (-(*this)(2, 0) * (*this)(3, 3) + (*this)(3, 0) * (*this)(2, 3)) * (*this)(1, 1) +
                    ((*this)(2, 1) * (*this)(3, 3) - (*this)(3, 1) * (*this)(2, 3)) * (*this)(1, 0)) * (*this)(0, 2) +
                   ((*this)(3, 0) * (*this)(1, 3) * (*this)(2, 2) - (*this)(2, 0) * (*this)(1, 3) * (*this)(3, 2) +
                    ((*this)(2, 0) * (*this)(3, 3) - (*this)(3, 0) * (*this)(2, 3)) * (*this)(1, 2) +
                    (-(*this)(2, 2) * (*this)(3, 3) + (*this)(2, 3) * (*this)(3, 2)) * (*this)(1, 0)) * (*this)(0, 1) +
                   (-(*this)(3, 1) * (*this)(1, 3) * (*this)(2, 2) + (*this)(2, 1) * (*this)(1, 3) * (*this)(3, 2) +
                    ((*this)(3, 1) * (*this)(2, 3) - (*this)(2, 1) * (*this)(3, 3)) * (*this)(1, 2) +
                    (*this)(1, 1) * ((*this)(2, 2) * (*this)(3, 3) - (*this)(2, 3) * (*this)(3, 2))) * (*this)(0, 0);
            // ... says Maple; this can possibly be done more efficiently,
            // maybe even with LU (with pivoting, though...)
        default:
            logger(ERROR, "Computing the Determinant for size % is not implemented", numberOfRows);
            break;
    }
    return 0;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns> SmallMatrix<numberOfRows, numberOfColumns>::inverse() const
{
    logger.assert_debug(numberOfRows == numberOfColumns, "Cannot invert a non-square matrix");
    SmallMatrix<numberOfRows, numberOfColumns> result = (*this);
    
    int nr = numberOfRows;
    int nc = numberOfColumns;
    
    int nPivot = numberOfRows;
    int iPivot[nPivot];
    
    int info = 0;
    
    dgetrf_(&nr, &nc, result.data(), &nr, iPivot, &info);
    
    int lwork = numberOfRows * numberOfColumns;
    SmallMatrix work;
    dgetri_(&nc, result.data(), &nc, iPivot, work.data(), &lwork, &info);
    
    return result;
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
template<unsigned int numberOfRightHandSideColumns>
void SmallMatrix<numberOfRows, numberOfColumns>::solve(SmallMatrix<numberOfRows, numberOfRightHandSideColumns>& B) const
{
    logger.assert_debug(numberOfRows == numberOfColumns, "can only solve for square matrixes");
    
    int n = numberOfRows;
    int nrhs = numberOfRightHandSideColumns;
    int info;
    
    int IPIV[numberOfRows];
    SmallMatrix<numberOfRows, numberOfColumns> matThis = *this;
    dgesv_(&n, &nrhs, matThis.data(), &n, IPIV, B.data(), &n, &info);
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
void SmallMatrix<numberOfRows, numberOfColumns>::solve(SmallVector<numberOfRows>& b) const
{
    logger.assert_debug(numberOfRows == numberOfColumns, "can only solve for square matrixes");
    
    int n = numberOfRows;
    int nrhs = 1;
    int info;
    
    int IPIV[numberOfRows];
    SmallMatrix matThis = *this;
    dgesv_(&n, &nrhs, matThis.data(), &n, IPIV, b.data(), &n, &info);
}

template<unsigned int numberOfRows, unsigned int numberOfColumns>
SmallVector<numberOfColumns> operator*(SmallVector<numberOfRows>& vec, SmallMatrix<numberOfRows, numberOfColumns>& mat)
{
    if (numberOfColumns == 0)
    {
        logger(WARN, "Trying to multiply a vector with a matrix without any columns.");
        return SmallVector<numberOfColumns>();
    }
    int nr = numberOfRows;
    int nc = numberOfColumns;
    
    int i_one = 1;
    double d_one = 1.0;
    double d_zero = 0.0;
    
    SmallVector<numberOfColumns> result;
    
    logger(DEBUG, "Matrix size: % x % \n Vector size: %", nr, nc, vec.size());
    
    dgemv_("T", &nr, &nc, &d_one, mat.data(), &nr, vec.data(), &i_one, &d_zero, result.data(), &i_one);
    return result;
}

