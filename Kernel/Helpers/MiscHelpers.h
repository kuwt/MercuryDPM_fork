<<<<<<<< HEAD:Kernel/Helpers/MiscHelpers.h
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

// Miscellaneous helpers

#ifndef MERCURYDPM_MISC_HELPERS_H
#define MERCURYDPM_MISC_HELPERS_H

#include "Math/ExtendedMath.h"

#include <string>

namespace helpers
{
    /*!
     * \brief Plots to a gnuplot window.
     */
    void gnuplot(std::string command);

    Mdouble getRealTime();

    /*!
     * \brief For use with qsort.
     */
    int qSortCompare(const void* x, const void* y);

    /*!
     * \brief Returns the 100*perc-th percentile of array.
     */
    double getPercentile(const double* array, size_t nel, double perc);
}

#endif // MERCURYDPM_MISC_HELPERS_H
========
//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
// For the list of developers, see <http://www.MercuryDPM.org/Team>.
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


// SimpleOpt feature, opt_main.cpp: Runs the python optimization tool with
// proper mercury dir structure passed via command line

#include"CMakeDefinitions.h"
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include "Walls/InfiniteWall.h"
#include <stdio.h>

// Helper function to implement extra OS commands after the driver code is done
std::string ExecCommand(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        std::string ret = "Could not run an extra OS command: "; ret.append(cmd);
        throw std::runtime_error(ret);
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

// Runs python optimization tool
int main(int argc, char** argv)
{
    std::string command;
    command = "python3 " + getMercurySourceDir() + "/Tools/SimpleOpt/BatchRun.py " + getMercurySourceDir() + " " + getMercuryBuildDir();
    std::cout<<ExecCommand(command.c_str())<<std::endl;
    return 0;
}
>>>>>>>> feature/Clumps:Drivers/Clump/Domino/DominoBatch.cpp
