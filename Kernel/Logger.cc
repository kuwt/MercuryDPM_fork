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

#include "Logger.h"

#ifdef MERCURY_STACKTRACE_SHOW
//To create stacktraces, we need the
//  backtrace(3)
//function calls.. (really, we don't want to do this by hand)
#include <execinfo.h>

#ifdef MERCURY_STACKTRACE_DEMANGLE
//However, we would end up with mangled function names...
//So, instead we'll use the abi::__cxa_demangle function call..
#include <cxxabi.h>
//Oh, and we really want to get the function names as well..
#include <dlfcn.h>
#endif

#endif

#include <cstdlib>
#include <iostream>
#include <csignal>

/*!
 *  \brief Definition of the different loglevels by its wrapper class LL. These are used as tags in template
 *  metaprogramming for the Logger class.
 */
LL<Log::FATAL> FATAL;
LL<Log::ERROR> ERROR;
LL<Log::WARN> WARN;
LL<Log::INFO> INFO;
LL<Log::DEFAULT> DEFAULT;
LL<Log::VERBOSE> VERBOSE;
LL<Log::DEBUG> DEBUG;


/*!
 *  \brief Definition of different loggers with certain modules. A user can define its own custom logger here.
 */
/* Actual definition of the default logger. */
Logger<MERCURY_LOGLEVEL> logger("MercuryKernel");
/* Actual definition of the default logger. */
Logger<CG_LOGLEVEL> cgLogger("MercuryCG");

/*!
 * \brief Prints messages of loglevel INFO.
 *
 * \param[in] module        The module name of the current logger invocation.
 * \param[in] msg           formatted message to be printed.
 * \param[in] doFlush       Flusher enum class object to enable/disable flushing of the output.
 */
static void printInfo(std::string module, std::string msg, Flusher doFlush)
{
#ifdef MERCURY_USE_MPI
    //Check if MPI is initialised
    initialiseMPI();
    MPIContainer& communicator = MPIContainer::Instance();
    std::cout << "[Process: " << communicator.getProcessorID() << "]: " << msg;
    if (doFlush == Flusher::FLUSH)
    {
        std::cout << std::endl;
    }
#else
    std::cout << msg;
    if (doFlush == Flusher::FLUSH)
    {
        std::cout << std::endl;
    }
#endif
}

/*!
 * \brief Prints messages of loglevel WARN.
 *
 * \param[in] module        The module name of the current logger invocation.
 * \param[in] msg           formatted message to be printed.
 * \param[in] doFlush       Flusher enum class object to enable/disable flushing of the output.
 */
static void printMessage(std::string module, std::string msg, Flusher doFlush)
{
#ifdef MERCURY_USE_MPI
    //Check if MPI is initialised
    initialiseMPI();
    MPIContainer& communicator = MPIContainer::Instance();
    std::cout << "\033[1;33mModule " << module << ":\033[0m\n" << "[Processor: " << communicator.getProcessorID() << "]" << msg;
    if (doFlush == Flusher::FLUSH)
    {
        std::cout << std::endl;
    }
#else
    std::cout << "\033[1;33mMessage " << module << ":\033[0m\n" << msg;
    if (doFlush == Flusher::FLUSH)
    {
        std::cout << std::endl;
    }
#endif
}

/*!
 * \brief Prints messages of loglevel ERROR.
 *
 * \param[in] module        The module name of the current logger invocation.
 * \param[in] msg           formatted message to be printed.
 * \param[in] doFlush       Flusher enum class object to enable/disable flushing of the output.
 */
[[noreturn]] static void printError(std::string module, std::string msg, Flusher doFlush)
{
#ifdef MERCURY_USE_MPI
    //Check if MPI is initialised
    initialiseMPI();
    MPIContainer& communicator = MPIContainer::Instance();
    std::cout << "\033[1;33mError " << module << ":\033[0m\n" << "[Processor: " << communicator.getProcessorID() << "]" << msg << std::endl;
#else
    std::cout << "\033[1;31mAn error has occured"
              << "\n\033[1;31mModule  :" << module
              << "\n\033[1;31mMessage :" << msg << std::endl;
#endif
#ifdef MERCURY_STACKTRACE_SHOW
    std::cerr << "\n-----------------[Stack Trace]-----------------\n";
    
    void* stackBuffer[64]; //This should be enough for all purposes..
    //First, we retrieve the addresses of the entire stack...
    int nStackFrames = backtrace(stackBuffer, 64);
#ifndef MERCURY_STACKTRACE_DEMANGLE
    //We don't have the demangling infra, so just use backtrace_symbols.
    char** functionNames = backtrace_symbols(stackBuffer, nStackFrames);
    for( int i = 0; i < nStackFrames; i++ )
    {   
        std::cerr << '\t' << functionNames[i] << '\n';
    }
    std::cerr << "Exiting.\n" << std::endl;

    //DO NOT USE DELETE HERE. THIS SHOULD BE free()'d!
    // -- dducks
    free(functionNames);
#else
    //We request the symbol information ourselves, in order to be able to demangle it.
    //And request the function names using dladdr.
    Dl_info infoStruct;
    for (int i = 4; i < nStackFrames; i++)
    {
        if (dladdr(stackBuffer[i], &infoStruct))
        { // We succesfully loaded the address...
            int demangleStatus;
            char* fnDemangled = abi::__cxa_demangle(infoStruct.dli_sname, NULL, NULL, &demangleStatus);
            if (infoStruct.dli_sname == nullptr)
                continue;
            
            //We even succesfully demangled the symbol...
            if (demangleStatus == 0)
            {
                std::cerr << fnDemangled << " +" << (void*) ((char*) stackBuffer[i] - (char*) infoStruct.dli_saddr) << "\t(" << infoStruct.dli_fname << ")\n";
                free(fnDemangled);
            }
            else
            { //Well, we tried. Lets output at least our raw symbol name.
                std::cerr << infoStruct.dli_sname << " +" << (void*) ((char*) stackBuffer[i] - (char*) infoStruct.dli_saddr) << "\t(" << infoStruct.dli_fname << ")\n";
            }
        }
        else
        { //Name lookup failed.
            std::cerr << stackBuffer[i] << ": ?????" << std::endl;
        }
    }
#endif
#endif
    //send a signal first, in case a debugger can catch it
    std::raise(SIGTERM);
    //call exit() for the specific loglevel
    std::exit(static_cast<int>(Log::ERROR));
}

/*!
 * \brief Prints messages of loglevel FATAL.
 *
 * \param[in] module        The module name of the current logger invocation.
 * \param[in] msg           formatted message to be printed.
 * \param[in] doFlush       Flusher enum class object to enable/disable flushing of the output.
 */
// Default implementation for logging errors / fatals
// [[noreturn]] indicates this function may not return
[[noreturn]] static void printFatalError(const std::string& module, const std::string& msg, Flusher doFlush)
{
#ifdef MERCURY_USE_MPI
    //Check if MPI is initialised
    initialiseMPI();
    MPIContainer& communicator = MPIContainer::Instance();
    std::cout << "\033[1;33mError " << module << ":\033[0m\n" << "[Processor: " << communicator.getProcessorID() << "]" << msg << std::endl;
#else
    std::cout << "\033[1;31mA fatal error has occured"
              << "\n\033[1;31mModule  :" << module
              << "\n\033[1;31mMessage :" << msg << std::endl;
#endif
#ifdef MERCURY_STACKTRACE_SHOW
    std::cerr << "\n-----------------[Stack Trace]-----------------\n";
    
    void* stackBuffer[64]; //This should be enough for all purposes..
    //First, we retrieve the addresses of the entire stack...
    int nStackFrames = backtrace(stackBuffer, 64);
#ifndef MERCURY_STACKTRACE_DEMANGLE
    //We don't have the demangling infra, so just use backtrace_symbols.
    char** functionNames = backtrace_symbols(stackBuffer, nStackFrames);
    for( int i = 0; i < nStackFrames; i++ )
    {
        std::cerr << '\t' << functionNames[i] << '\n';
    }
    std::cerr << "Exiting.\n" << std::endl;

    //DO NOT USE DELETE HERE. THIS SHOULD BE free()'d!
    // -- dducks
    free(functionNames);
#else
    //We request the symbol information ourselves, in order to be able to demangle it.
    //And request the function names using dladdr.
    Dl_info infoStruct;
    for (int i = 4; i < nStackFrames; i++)
    {
        if (dladdr(stackBuffer[i], &infoStruct))
        { // We succesfully loaded the address...
            int demangleStatus;
            char* fnDemangled = abi::__cxa_demangle(infoStruct.dli_sname, NULL, NULL, &demangleStatus);
            if (infoStruct.dli_sname == nullptr)
                continue;
            
            //We even succesfully demangled the symbol...
            if (demangleStatus == 0)
            {
                std::cerr << fnDemangled << " +" << (void*) ((char*) stackBuffer[i] - (char*) infoStruct.dli_saddr) << "\t(" << infoStruct.dli_fname << ")\n";
                free(fnDemangled);
            }
            else
            { //Well, we tried. Lets output at least our raw symbol name.
                std::cerr << infoStruct.dli_sname << " +" << (void*) ((char*) stackBuffer[i] - (char*) infoStruct.dli_saddr) << "\t(" << infoStruct.dli_fname << ")\n";
            }
        }
        else
        { //Name lookup failed.
            std::cerr << stackBuffer[i] << ": ?????" << std::endl;
        }
    }
#endif
#endif
    //send a signal first, in case a debugger can catch it
    std::raise(SIGTERM);
    //call exit for the specific loglevel
    std::exit(static_cast<int>(Log::FATAL));
}

// Default output methods.
LoggerOutput loggerOutputDefaultImpl = {printFatalError, //onFatal
                                        printError, //onError
                                        printMessage, //onWarn
                                        printInfo, //onInfo
                                        printInfo, //onVerbose
                                        printMessage //onDebug
};

//And we assign them.
LoggerOutput* loggerOutput = &loggerOutputDefaultImpl;
