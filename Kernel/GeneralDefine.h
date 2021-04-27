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

#ifndef GENERALDEFINE_H
#define GENERALDEFINE_H

#include <limits>

#ifdef HIGH_PRECISION
typedef long double Mdouble;
#else
typedef double Mdouble;
#endif

#define MERCURY_DEPRECATED [[deprecated]]

#define UNUSED  __attribute__ ((__unused__))

namespace constants
{
    const Mdouble NaN = std::numeric_limits<Mdouble>::quiet_NaN();
    const Mdouble inf = std::numeric_limits<Mdouble>::infinity();
    const int intMax = std::numeric_limits<int>::max();
    const unsigned unsignedMax = std::numeric_limits<unsigned>::max();
}

//The number of digits are important for MAX_PROC, not the actual number. It is used to determine unique communication tags
//IMPORTANT: Always take a value 10^a where a is a natural number above 0
#define MAX_PROC 1000

/*!
 * \brief For the MPI communication routines this quantity is often required.
 * defining this macro makes the code a bit cleaner. Sadly it can't be defind
 * as a constant global variable as on compile time it is not known what it should be
 */
#ifdef MERCURY_USE_MPI
#define NUMBER_OF_PROCESSORS static_cast<unsigned>(MPIContainer::Instance().getNumberOfProcessors())
#define PROCESSOR_ID MPIContainer::Instance().getProcessorID()
#else
#define NUMBER_OF_PROCESSORS 1
#define PROCESSOR_ID 0
#endif

#ifdef MERCURY_USE_OMP
#define OMP_THREAD_NUM omp_get_thread_num()
#else
#define OMP_THREAD_NUM 0
#endif

/*!
 * \brief An enum that indicates the direction in Cartesian coordinates
 * \details This is heavily used in the parallel code to avoid confusion with numbers 
 */
enum Direction
{
    XAXIS = 0, YAXIS = 1, ZAXIS = 2
};


#endif
