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

#include "Helpers/MiscHelpers.h"
#include "Logger.h"

#include <chrono>

#ifdef __linux__
#include <unistd.h>
#endif

/*!
 * \note Not used function.
 */
void helpers::gnuplot(std::string command)
{
#ifdef __CYGWIN__
    logger(WARN, "[helpers::gnuplot] is not supported on Cygwin");
#elif _WIN32
    logger(WARN, "[helpers::gnuplot] is not supported on Windows");
#else
    FILE* pipe = popen("gnuplot -persist", "w");
    fprintf(pipe, "%s", command.c_str());
    fflush(pipe);
#endif
}

Mdouble helpers::getRealTime()
{
    // record start time
    static auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    return diff.count();
}

int helpers::qSortCompare(const void* x, const void* y)
{
    double xx = *(double*) x, yy = *(double*) y;
    if (xx < yy) return -1;
    if (xx > yy) return 1;
    return 0;
}

/*!
 * \brief Returns the 100*perc-th percentile of array.
 * \details array should be sorted, e.g.
 *  qsort(xs, n, sizeof(double), qSortCompare);
 *  and perc should be a number between 0 and 1.
 */
double helpers::getPercentile(const double* array, size_t nel, double perc)
{
    size_t lower_ind = floor(perc * (nel - 1));
    size_t upper_ind = ceil(perc * (nel - 1));
    double lambda = (perc * (nel - 1)) - lower_ind;
    double lower_x = array[lower_ind];
    double upper_x = array[upper_ind];
    double percentile = (1 - lambda) * lower_x + lambda * upper_x;
    return percentile;
}
