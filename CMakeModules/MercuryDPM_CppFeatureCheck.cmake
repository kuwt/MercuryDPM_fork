#Here is the check for C++ support of different versions
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

#Add the flag to compiler options if the current compiler is found
if (COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
elseif (COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif (COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no CXX11, CXX14 or CXX17 support. We require a compiler
	that has support for CXX11, CXX14 or CXX17 for example GCC 5.1 or later")
endif ()

if (CMAKE_COMPILER_IS_GNUCC)
    #check for GCC if it is correct version or not. The version is the minimum version supporting CXX17
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
        message(FATAL_ERROR "Insufficient gcc version, you need at least version 5.1 of GCC")
    endif ()

endif ()