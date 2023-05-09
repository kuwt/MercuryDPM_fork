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

#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <sstream>
#include <functional>
#include <type_traits>
#include <iomanip>
#include "GeneralDefine.h"
#include <iostream>

#ifndef MERCURYDPM_LOGLEVEL
#define MERCURYDPM_LOGLEVEL Log::DEFAULT
#endif

#ifndef CG_LOGLEVEL
#define CG_LOGLEVEL Log::DEFAULT
#endif

#ifdef assert
#undef assert
//#error You included assert before the logger. Please use logger.assert_debug() instead.
#endif

#ifdef MERCURYDPM_FORCE_ASSERTS
#define MERCURYDPM_ASSERTS true
#else
#ifdef MERCURYDPM_NO_ASSERTS
#define MERCURYDPM_ASSERTS false
#else
#ifdef NDEBUG
#define MERCURYDPM_ASSERTS false
#else
#define MERCURYDPM_ASSERTS true
#endif
#endif
#endif

/*!
 * \brief The Logger class provides ability to write log messages in your own customized format
 *
 * \details logger objects can be globally called on different loglevels whose usage is briefly explained where they
 * have been declared down below.
 *
 * The previous version of the Logger used the syntax
 * Logger<Log::LEVEL> logger;
 * which is now predefined in logger.cc.
 *
 * The operator() function is utilized in this logger to keep log messages in a single statement of format:
 *
 *      logger(Log::LEVEL, "Message", args...);
 *
 * Where args... is a user defined pack of parameters. Each parameter is denoted by a percentage character %, see
 * example below.
 *
 * logger(INFO, "The current timestep is %.", getTimeStep());
 * OUTPUT: The current timestep is 0.0104124.
 *
 * The output of arguments an also be customized by a certain precision or width. A number following the percentage
 * character is defined as the precision of the parameter (e.g. %12), whereas a full stop followed by a number is defined as
 * the output width of a parameter (e.g. %.10). Here std::left was chosen as default positioning. See examples below.
 *
 * logger(INFO, "The current timestep is %12.", getTimestep());
 * OUTPUT: The current timestep is 0.01020304917.
 *
 * logger(INFO, "The current timestep is %.12.", getTimestep());
 * OUTPUT: The current timestep is 0.0102051    .
 *
 * logger(INFO, "The current timestep is %12.14.", getTimestep());
 * OUTPUT: The current timestep is 0.01020304917  .
 *
 * After each invocation of the logger the output is flushed. To speed-up the code it is key to avoid flushing of the
 * output wherever possible. Therefore, escape characters such as "\n" are preferred over std::endl. Furthermore if a
 * parameter pack of the logger contains the argument Flusher::NO_FLUSH it will skip the use of std::endl after logger
 * invocation.
 *
 * Default loglevel is Log::DEFAULT = 0.
 */


/*!

 * \brief Enum class which enables/disables flushing of output.
 * \details the Flusher class is only used if the loglevel is below VERBOSE and DEBUG and if the CMAKE_BUILD_TYPE is
 * not "Debug".
 * If Flusher::NO_FLUSH is added as an argument to a logger invocation it will prevent the logger to call std::endl
 * and flush the ouput at the end of the call. This is mainly utilized to speed up the code when logging messages.
 */
enum class Flusher
{
    FLUSH,
    NO_FLUSH
};

/*!
 * \brief The different loglevels.
 *
 * \details The different loglevels, represented as signed characters,
 * in descending order of severeness. Worst is FATAL, best is DEBUG.
 *
 * Please, use the tags FATAL/ERROR/etc without class/enum/namespace.
 */
enum class Log
        : signed char
{
    FATAL = -20, ERROR = -15, WARN = -10, INFO = -5, DEFAULT = 0, VERBOSE = 5, DEBUG = 10
};

/*!
 * \brief Internally used to filter on loglevel.
 * Do not edit, as this is required for an optimised logger.
 * \details used for the implementation of different loglevels; the operator compares current loglevel to a fixed one
 * and thus checks if a message should be logged.
 */
constexpr bool operator<=(const Log rhs, const Log lhs)
{
    return ((static_cast<signed char>(rhs)) <= (static_cast<signed char>(lhs)));
}

/*!
 * \brief Default functions for output generation
 *
 *  \details These handlers will be called on generation of the message.
 * The functions are of signature
 *    void (std::string module Name, std::string message, Flusher doFlush_);
 *
 * These functions may not return but call std::exit() instead.
 * They may also throw any exception to allow code to gracefully
 * recover.
 */
class LoggerOutput
{
public:
    std::function<void(std::string, std::string, Flusher)> onFatal;
    std::function<void(std::string, std::string, Flusher)> onError;
    std::function<void(std::string, std::string, Flusher)> onWarn;
    std::function<void(std::string, std::string, Flusher)> onInfo;
    std::function<void(std::string, std::string, Flusher)> onVerbose;
    std::function<void(std::string, std::string, Flusher)> onDebug;
};

/*!
 * \brief Declaration of the output functions.
 * \details If the output needs to be redirected, please
 * swap the loggerOutput pointer to your preferred
 * LoggerOutput instance, and make sure this exists
 * until _AFTER_ an std::exit() invocation.
 * (e.g. clean up with std::atexit())
 */
extern LoggerOutput* loggerOutput;

// Forward declaration..
template<Log L = Log::DEFAULT, bool ASSERTS = MERCURYDPM_ASSERTS>
class Logger;

/*!
 * \brief Tag for template metaprogramming
 *
 * This tag class serves as a way to correctly
 * resolve the implementation (if any) of the
 * operator() or .log() method of the Logger.
 * Please, don't change it at all nor give it
 * any members.
 *
 * \details The LL class is a wrapper for the Log::LogLevel class. This wrapper is only used to enable the loglevels
 * to be written as single words.
 */
template<Log Level>
class LL
{
public:
};

/*!
 * The following are the loglevels which should be used to control the logger. They are declared here but defined in
 * Logger.cc.
 */

/*!
 * \brief Fatal log level
 * 
 * Fatal, as in, the program has suffered from the worst possible failure and there is no
 * way it can gracefully recover. The difference to ERROR is mainly that this type of failure is most often caused by
 * non-human error
 *
 * Example: No memory allocations possible
 *
 * Default behaviour: log to std::cerr, followed by std::exit().
 */
extern LL<Log::FATAL> FATAL;

/*!
 * \brief Error log level
 *
 * Error, as in, the program has found a severe problem which it cannot resolve
 * any further. It does not know how to recover in a sane way. The difference to FATAL is mainly that this type of
 * failure is most often caused by human error.
 *
 * Example: Negative time step, Infinite end time and no override of the
 * continuation function.
 *
 * Default behaviour: log to std::cerr, followed by std::exit().
 */
extern LL<Log::ERROR> ERROR;

/*!
 * \brief Warning log level
 *
 * Warning, as in, the program has detected a problem but does know a solution.
 * The simulation can continue with this fix, but the user should look at
 * fixing his / her simulation so this won't occur in the future.
 *
 * Example: Setting a smaller Xmax than Xmin.
 *
 * Default behaviour: log to std::cerr, returns afterwards.
 */
extern LL<Log::WARN> WARN;

/*!
 * \brief Info log level
 *
 * Useful information, small oddities and statistics which should be of no real effect
 * to the user, but still give useful information about the current state and progress of the program.
 *
 * Example: Finished inserting 381 particles.
 *
 * Default behaviour: log to std::cout, returns afterwards.
 */
extern LL<Log::INFO> INFO;

/*!
 * \brief Default log level
 *
 * Only useful for defining the loglevel of the logger itself. Should not actually be used.
 */
extern LL<Log::DEFAULT> DEFAULT;

/*!
 * \brief Verbose information
 *
 * Information which is not useful to anybody except those looking for weird behaviour.
 * These should however still be clear in meaning.
 *
 * Example: Finished creating a particle.
 *
 * Default behaviour: ignore.
 */
extern LL<Log::VERBOSE> VERBOSE;

/*!
 * \brief Debug information
 *
 * Only used for internal development. Can be very cryptic, as it is only meant for finding
 * bugs / oddities by the internal development team.
 * 
 * Example: Collision found between Particle #38201 and Wall #5
 *
 * Default behaviour: ignore.
 */
extern LL<Log::DEBUG> DEBUG;

/*!
 * \brief the Logger class is the main class of the logger implementation. It holds all the functions which invoke
 * certain methods to create messages based on input parameter deductions.
 *
 * \tparam L The log level defined in cMake configuration. Messages of higher level than L are ignored.
 * 
 * Usage: logger(FATAL, "Error in (here) because % < %!\n", var1, var2)
 * OUTPUT: Error in (here) because 2 < 1!
 *
 *
 * Define custom loggers by:
 * #ifndef HG_LOGLEVEL_CUSTOMMOD
 * #define HG_LOGLEVEL_CUSTOMMOD Log::Debug
 * #endif
 * Logger<HG_LOGLEVEL_CUSTOMMOD> customLogger; 
 */
template<Log L, bool ASSERTS>
class Logger
{

private:
    /*!
     * \brief The module name of this actual logger.
     */
    const std::string module;
    /*!
     * \brief Can prevent the logger from flushing the buffer via std::endl. doFlush_ is set automatically based on
     * build and loglevel settings.
     */
    Flusher doFlush_ = Flusher::FLUSH;

public:
    
    /*!
     * \brief constructor
     *
     * \param[in] name   The name in this module used in output messages.
     */
    explicit Logger(const std::string name)
            : module(name)
    {
    }
    
    /*!
     * \brief destructor
     */
    ~Logger()
    = default;
    
    /*!
     *
     * \brief Log implementation of this function
     *
     * Actual implementation of the log function.
     * If the user defined loglevel L is lower than the called LOGLEVEL it will evaluate to an empty body
     * function below. If L is greater than the called LOGLEVEL it will invoke this function.
     *
     * \param[in] log        Loglevel, either FATAL, ERROR, WARN, INFO, VERBOSE, DEBUG
     * \param[in] format     Message format, where % can be used as a placeholder for arguments.
     * \param[in] arg     Any arguments which replace all the % characters.
     */
    
    template<Log LOGLEVEL, typename ... Args>
    typename std::enable_if<!((L < LOGLEVEL) && (MERCURYDPM_LOGLEVEL < LOGLEVEL)), void>::type
    operator()(const LL<LOGLEVEL> log, const char* format UNUSED, Args&& ... arg UNUSED)
    {
        std::stringstream msgstream;
        createMessage(msgstream, format, arg...);
        if (LOGLEVEL <= Log::FATAL)
        {
            loggerOutput->onFatal(module, msgstream.str(), doFlush_);
        }
        else if (LOGLEVEL <= Log::ERROR)
        {
            loggerOutput->onError(module, msgstream.str(), doFlush_);
        }
        else if (LOGLEVEL <= Log::WARN)
        {
            loggerOutput->onWarn(module, msgstream.str(), doFlush_);
            doFlush_ = Flusher::FLUSH;
        }
        else if (LOGLEVEL <= Log::INFO)
        {
            loggerOutput->onInfo(module, msgstream.str(), doFlush_);
            doFlush_ = Flusher::FLUSH;
        }
        else if (LOGLEVEL <= Log::VERBOSE)
        {
            loggerOutput->onVerbose(module, msgstream.str(), doFlush_);
            doFlush_ = Flusher::FLUSH;
        }
        else
        {
            loggerOutput->onDebug(module, msgstream.str(), doFlush_ = Flusher::FLUSH);
            doFlush_ = Flusher::FLUSH;
        }
    }
    
    /*!
     * \brief Empty body function utilized to suppress logger messages above a certain user defined loglevel L.
     */
    template<Log LOGLEVEL, typename... Args>
    typename std::enable_if<L < LOGLEVEL && MERCURYDPM_LOGLEVEL < LOGLEVEL, void>::type
    operator()(const LL<LOGLEVEL> log, const char* format UNUSED, Args&& ... arg UNUSED)
    {
    }
    
    /*!
     * \brief Converts a std::string message into a char array terminated by a null character (C-string).
     */
    //std::string is sometimes convenient, but always slow, so where possible, don't convert the const char* to a string
    //before converting it back
    template<Log LOGLEVEL, typename... Args>
    void operator()(const LL<LOGLEVEL> log, const std::string& format UNUSED, Args&& ... arg
    UNUSED)
    {
        (*this)(log, format.c_str(), arg...);
    }
    
    /*!
     *
     * \brief Asserts on this logger
     *
     * If ASSERTS are activated, evaluates an assertion, prints an error message and aborts
     * in case of a failure. This message can be redirected and will be send to loggerOuptput->onFatal.
     * If ASSERTS are not activated it will redirect to the empty body function below.
     *
     * \param[in] assertion       An assertion, which must be true.
     * \param[in] format          Message format, where % can be used as a placeholder for arguments.
     * \param[in] arg             Any arguments which replaces all the % characters.
     */
    
    
    template<typename... Args>
    typename std::enable_if<(ASSERTS) && (sizeof...(Args) >= 0), void>::type
    assert_debug(bool assertion, const char* format, Args&& ... arg)
    {
        assert_always(assertion, format, arg...);
    }
    
    template<typename... Args>
    typename std::enable_if<!((ASSERTS) && sizeof...(Args) >= 0), void>::type
    assert_debug(bool assertion, const char* format, Args&& ... arg)
    {
    }
    
    /*!
     * \brief Converts a std::string message into a char array terminated by a null character (C-string).
     */
    //the conversion from "" to a std::string is so slow, it takes 50% of the total run time for a release build...
    template<typename... Args>
    void assert_debug(bool assertion, const std::string format, Args&& ... arg)
    {
        assert_debug(assertion, format.c_str(), arg...);
    }
    
    template<typename... Args>
    void assert_always(bool assertion, const char* format, Args&& ... arg)
    {
        if (!assertion)
        {
            std::stringstream msgstream;
            createMessage(msgstream, format, arg...);
            loggerOutput->onFatal(module, msgstream.str(), doFlush_);
        }
        
    }
    
    /*!
     * \brief Converts a std::string message into a char array terminated by a null character (C-string).
     */
    template<typename... Args>
    void assert_always(bool assertion, const std::string format, Args&& ... arg)
    {
        assert_always(assertion, format.c_str(), arg...);
    }
    
    /*!
     * \brief Oldschool log method.
     * \deprecated Use operator() instead.
     */
    template<typename... Args>
    MERCURYDPM_DEPRECATED
    void log(const Log loglevel, const std::string& format, Args&& ... arg)
    {
        if (loglevel <= L || loglevel <= MERCURYDPM_LOGLEVEL)
        {
            std::stringstream msgstream;
            createMessage(msgstream, format.c_str(), arg...);
            if (loglevel <= Log::FATAL)
            {
                loggerOutput->onFatal(module, msgstream.str(), doFlush_);
            }
            else if (loglevel <= Log::ERROR)
            {
                loggerOutput->onError(module, msgstream.str(), doFlush_);
            }
            else if (loglevel <= Log::WARN)
            {
                loggerOutput->onWarn(module, msgstream.str(), doFlush_);
                doFlush_ = Flusher::FLUSH;
            }
            else if (loglevel <= Log::INFO)
            {
                loggerOutput->onInfo(module, msgstream.str(), doFlush_);
                doFlush_ = Flusher::FLUSH;
            }
            else if (loglevel <= Log::VERBOSE)
            {
                loggerOutput->onVerbose(module, msgstream.str(), doFlush_);
                doFlush_ = Flusher::FLUSH;
            }
            else
            {
                loggerOutput->onDebug(module, msgstream.str(), doFlush_ = Flusher::FLUSH);
                doFlush_ = Flusher::FLUSH;
            }
        }
    }

private:
    
    /*!
     *
     * \brief Edits the message to a certain format and writes it to a stringstream by recursively replacing all %
     * characters with the arguments values.
     * \details The creation of messages is divided into three different overloaded functions. the function
     * createMessage is recursively called and each of the functions below is called for a certain case dependent on
     * the amount and type of parameters.
     *
     * \param[in] msg       stringstream which represents the output message.
     * \param[in] fmt       char array of the yet unformatted message.
     * \param[in] arg       argument to replace the next % character.
     * \param[in] args      parameter pack of the remaining arguments.
     *
     */
    template<typename Arg1, typename... Args>
    void createMessage(std::stringstream& msg, const char* fmt,
                       Arg1&& arg, Args&& ... args)
    {
        bool doSkipNext = false;
        while (*fmt != '%' || doSkipNext)
        {
            //Make sure we're not running past the end of our formatting string.
            if (*fmt == '\0')
                return;
            
            if (*fmt == '\\' && !doSkipNext)
            { //Escape for the % character
                doSkipNext = true;
                fmt++;
            }
            else
            {
                msg << *fmt;
                fmt++;
                doSkipNext = false;
            }
        }
        fmt++; //Consume the % character
        int precision = 0;
        int width = 0;
        // if precision and width or only precision is defined
        if (isdigit(*fmt))
        {
            precision = std::atoi(fmt);
            while (isdigit(*fmt))
            {
                fmt++;
            }
            if (std::ispunct(*fmt))
            {
                fmt++;
                if (std::isdigit(*fmt))
                {
                    width = std::atoi(fmt);
                    while (isdigit(*fmt))
                    {
                        fmt++;
                    }
                }
                    // else the char is a real full stop so set the pointer back to full stop.
                else
                {
                    fmt--;
                }
            }
        }
            // if only a width and no precision defined
        else if (std::ispunct(*fmt))
        {
            fmt++;
            if (std::isdigit(*fmt))
            {
                width = std::atoi(fmt);
                while (isdigit(*fmt))
                {
                    fmt++;
                }
            }
                // else the char is a real full stop so set the pointer back to full stop.
            else
            {
                fmt--;
            }
        }
        if (width != 0 && precision != 0)
        {
            msg << std::setprecision(precision) << std::left << std::setw(width) << arg;
        }
        else if (precision != 0)
        {
            msg << std::setprecision(precision) << arg;
        }
        else if (width != 0)
        {
            msg << std::left << std::setw(width) << arg;
        }
        else
        {
            msg << arg;
        } //include args somehow..
        createMessage(msg, fmt, args...);//and recursively call ourselve / the method below.
    }
    
    
    /*!
     * \brief Overloaded version of createMessage to catch arguments of Flusher and suppress input flushing via
     * std::endl. If there is an argument which should be catched from the logger, overloading the function is the
     * way to go.
     *
     * \param[in] msg       stringstream which represents the output message.
     * \param[in] fmt       char array of the yet unformatted message.
     * \param[in] arg       argument of type Flusher which will be skipped and does not replace the next % character.
     * \param[in] args      parameter pack of the remaining parameters.
     */
    //terminating case for Flusher not needed. This function is also called when the parameter pack Args&& args is
    // empty
    template<typename... Args>
    void createMessage(std::stringstream& msg, const char* fmt,
                       Flusher arg, Args&& ... args)
    {
        // only suppress flushing if Mercury is not in CMAKE_BUILD_TYPE "Debug" and if the user defined loglevel from
        // cMake is below VERBOSE/DEBUG (<=5)
#ifndef MERCURYDPM_DEBUG
        if (arg != Flusher::FLUSH && MERCURYDPM_LOGLEVEL <= Log::VERBOSE)
        {
            doFlush_ = Flusher::NO_FLUSH;
        }
#endif
        // skip this argument by recursively calling this function again
        createMessage(msg, fmt, args...);
    }
    
    
    /*!
     * \brief Terminating case / Argument call.
     * Overloaded function for a logger message with only one argument or where only one argument is left.
     *
     * \param[in] msg       stringstream which represents the output message.
     * \param[in] fmt       char array of the yet unformatted message.
     * \param[in] arg       argument to replace the next % character.
     */
    // faster than above function and non recursive that is why it is good to have it
    template<typename Arg1>
    void createMessage(std::stringstream& msg, const char* fmt, Arg1&& arg)
    {
        bool doSkipNext = false;
        while (*fmt != '%' || doSkipNext)
        {
            if (*fmt == '\0') // End of string
                return;
            
            if (*fmt == '\\' && !doSkipNext)
            { //Escape for the % character and the \ character
                doSkipNext = true;
                fmt++;
            }
            else
            { //invoke the replacement
                msg << *fmt;
                fmt++;
                doSkipNext = false;
            }
        }
        fmt++; //Consume the % character
        int precision = 0;
        int width = 0;
        // if precision and width or only precision is defined
        if (isdigit(*fmt))
        {
            precision = std::atoi(fmt);
            while (isdigit(*fmt))
            {
                fmt++;
            }
            if (std::ispunct(*fmt))
            {
                fmt++;
                if (std::isdigit(*fmt))
                {
                    width = std::atoi(fmt);
                    while (isdigit(*fmt))
                    {
                        fmt++;
                    }
                }
                    // else the char is a real full stop so set the pointer back to full stop.
                else
                {
                    fmt--;
                }
            }
        }
            // if only a width and no precision defined
        else if (std::ispunct(*fmt))
        {
            fmt++;
            if (std::isdigit(*fmt))
            {
                width = std::atoi(fmt);
                while (isdigit(*fmt))
                {
                    fmt++;
                }
            }
                // else the char is a real full stop so set the pointer back to full stop.
            else
            {
                fmt--;
            }
        }
        if (width != 0 && precision != 0)
        {
            msg << std::setprecision(precision) << std::left << std::setw(width) << arg << fmt;
        }
        else if (precision != 0)
        {
            msg << std::setprecision(precision) << arg << fmt;
        }
        else if (width != 0)
        {
            msg << std::left << std::setw(width) << arg << fmt;
        }
        else
        {
            msg << arg << fmt;
        }
    }
    
    
    /*!
     * \brief Terminating case / no argument call
     * Overloaded function for a logger message without arguments.
     *
     * \param[in] msg       stringstream which represents the output message.
     * \param[in] message   char array of the message.
     */
    void createMessage(std::stringstream& msg, const char* message)
    {
        msg << message;
    }
};

/*! Default logger.
 * Use this for general logging.
 *
 * For very specific modules, define your own logger.
 * If you want to make extensive use of Debug messages,
 * please use a custom logger as well, to prevent polluting
 * the output.
 */
extern Logger<MERCURYDPM_LOGLEVEL> logger;

extern Logger<CG_LOGLEVEL> cgLogger;

//just emptying the functions is not sufficiently aggressive in disabling the actual (costly) comparison
#if !MERCURYDPM_ASSERTS
#define assert(e,...) assert(true,"")
#endif

#ifdef MERCURYDPM_USE_MPI
#include "MpiContainer.h"
#endif

#endif
