#Here is the check for CX14 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)

#Add the flag to compiler options if the current compiler is found
if (COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. We require a compiler that has :support for CXX14, for example GCC 5.0 or later")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    #Extra check for GCC if it is correct version or not
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
		message(FATAL_ERROR "Insufficient gcc version, you need at least version 5.0 of GCC")
    endif()

endif()
