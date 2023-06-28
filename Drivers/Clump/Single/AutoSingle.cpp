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

// Automatic script to clean working directories, recompile the code and post-process the output for Paraview


#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include<CMakeDefinitions.h>

// This function executes the OS command 'cmd' and returns the 'result'
std::string exec_command(const char* cmd) {
    std::array<char, 256> buffer;
    std::string result;

    // Prepare a pipe to write a to execute a systm command and read the result
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}


int main(int argc, char* argv[])
{
    // Automatic script to clean working directories, recompile the code and post-process the output for Paraview
    std::string command;
    std::string name = "Single";

    // Remove data for stl sequence (Blender) visualizations
    //command = "rm ClumpSeq.txt";
    //exec_command(command.c_str());

    // Make
    command = "make " + name;
    exec_command(command.c_str());

    // Run
    command = "./" + name;
    exec_command(command.c_str());

    // Clean the old paraview output directory
    command = "rm -rf paraview_" + name;
    exec_command(command.c_str());

    // Create new paraview output directory
    command = "mkdir paraview_" + name;
    exec_command(command.c_str());

    // Data2pvd tool run
    command = "../../../Tools/data2pvd " + name + ".data paraview_" + name + "/" + name;
    exec_command(command.c_str());

    // Paraview energy data postprocessing tool
    command = "python " + getMercuryDPMSourceDir() + "/Tools/MClump/PlotEnergies.py " +
            getMercuryDPMBuildDir() + "/Drivers/Clump/" + name + "/ " + name;
    exec_command(command.c_str());
    return 0;
}
