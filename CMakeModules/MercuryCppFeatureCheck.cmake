#Here is the check for CX11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

#Add the flag to compiler options if the current compiler is found
if (COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif (COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()	
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. We require a compiler that has :support for CXX11, for example GCC 4.8 or later")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message("GCC 4.x does not support the DEPRECATED attribute; therefore, compiler warnings about setting the attribute are ignored")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes")
    endif()

    #Extra check for GCC if it is correct version or not
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
		message(FATAL_ERROR "Insufficient gcc version, you need at least version 4.8 of GCC")
    endif()

endif()