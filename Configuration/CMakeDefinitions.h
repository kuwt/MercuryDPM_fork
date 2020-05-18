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

/// This file is used for generating defitions that give access to CMakeVaribles from within a cpp file (defintions have only been added as required
#ifndef CMAKEDEFINTIONS_H
#define CMAKEDEFINTIONS_H

#include <string>

const std::string getMercurySourceDir();

const std::string getMercuryBuildDir();

const int getSVNRevision();

const std::string getSVNURL();

const std::string getVersion();

/*
* This maps the CMake variables to defines used in the Logger.
* Because CMake uses ON and OFF opposed to 0 and 1, we have to
* use some witchery to make it work too.
* We'll clean up the mess afterwards by undefining ON and OFF...
*/
#define ON 1
#define OFF 0
#if @Mercury_BACKTRACE_ENABLE@ == ON
// This symbol is only defined as true, when the stacktrace code should be compiled in.
#define MERCURY_STACKTRACE_SHOW 1
#else
#define MERCURY_STACKTRACE_SHOW 0
#endif

#if @Mercury_BACKTRACE_DEMANGLE@ == ON
// This symbol is only defined as true, when every system required for demangling is present.
#define MERCURY_STACKTRACE_DEMANGLE 1
#else
#define MERCURY_STACKTRACE_DEMANGLE 0
#endif
// Cleaning up our symbols.
#undef ON
#undef OFF

#endif
