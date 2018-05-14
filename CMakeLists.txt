cmake_minimum_required (VERSION 2.8)
project (Mercury C CXX Fortran)

set(Mercury_VERSION_MAJOR 0)
set(Mercury_VERSION_MINOR 10)

set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")


enable_testing()

if (UNIX)
    add_definitions( -DUNIX )
endif(UNIX)

if (WIN32)
    message("Windows operating system")
    add_definitions( -DWINDOWS )
endif(WIN32)

if (CYGWIN)
	message("You are using CYGWIN; which is no longer supported. However the code should still work")
	add_definitions( -DCYGWIN )
endif(CYGWIN)


include(CMakeModules/MercuryMPI.cmake)

include(CMakeModules/MercuryLogger.cmake)

include(CMakeModules/MercuryBacktrace.cmake)

include(CMakeModules/MercuryCppFeatureCheck.cmake)

include(CMakeModules/MercuryDocumentation.cmake)

include(CMakeModules/MercuryXballs.cmake)



#If user has not given build type set it to Release
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE 
	STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(Drivers)
add_subdirectory(Kernel)
add_subdirectory(Tools)
add_subdirectory(XBalls)

#Add the command fullTest, which rebuilds the cache, cleans the project and then compiles and tests everything.
add_custom_target(fullTest
    COMMAND make rebuild_cache
	COMMAND make clean
	COMMAND make
	COMMAND make test
	COMMENT "This will clean the project, rebuild and test." )

#This section generates a configuration file that contains the CMake variables#
#############################################

configure_file(Configuration/CMakeDefinitions.cc
	${PROJECT_SOURCE_DIR}/Kernel/CMakeDefinitions.cc)
configure_file(Configuration/CMakeDefinitions.h
	${PROJECT_SOURCE_DIR}/Kernel/CMakeDefinitions.h)


#This section is about the build of documentation#
##################################################

option(Mercury_BUILD_USER_DIR "Compile USER DIR in testing" OFF)



FIND_PACKAGE(Python)
configure_file(Configuration/fpdiff.py
	${PROJECT_SOURCE_DIR}/Scripts/fpdiff.py)

#find_program(PYTHON_EXEC "python" 
#	DOC "Location of python executable to use")






