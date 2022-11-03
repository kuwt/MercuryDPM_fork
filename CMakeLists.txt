cmake_minimum_required (VERSION 2.8.12)
project (Mercury C CXX Fortran)
#set(CMAKE_CXX_STANDARD 17)

# CMAKE_CXX_STANDARD is a cmake 3.1 feature
# CMAKE_CXX_STANDARD 17 is supported in 3.18, but not in 3.5
# set(CMAKE_CXX_STANDARD 17)

# Trunk will become 0.14.alpha
set(Mercury_VERSION_MAJOR 0)
set(Mercury_VERSION_MINOR 14)

set(REVISION, "0")
set(REPOSITORY_URL "")

find_package(Git)
if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
	if(Git_FOUND)
		execute_process(
			COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
			WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
			OUTPUT_VARIABLE REVISION
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
		# For Git it doesn't make sense to specify a URL since the
		# commit might not have been pushed to a server. (Using
		# `git remote -v` doesn't help if the repository has more than
		# one remote set up.)
		set(REPOSITORY_URL "https://bitbucket.org/mercurydpm/mercurydpm.git")

		message("Git: commit is ${REVISION}")

	else()
		message("The project is in a Git repository, but the git executable is not available")
	endif()
else()
	message("This Mercury version is not in Git")
	find_package(Subversion)
	if(EXISTS ${PROJECT_SOURCE_DIR}/.svn)
		if(Subversion_FOUND)
			Subversion_WC_INFO(${PROJECT_SOURCE_DIR} SVN)
			set(REPOSITORY_URL ${SVN_WC_URL})
			set(REVISION ${Project_WC_REVISION})
		endif()
	else()
		message("This Mercury version in not in svn")
	endif()
endif()

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
    add_definitions(-DCYGWIN)
endif (CYGWIN)

#If user has not given build type set it to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE
            STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#If cMake builds in debug mode add the MERCURY_DEBUG flag
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DMERCURY_DEBUG)
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

option(Mercury_TRIANGLE_WALL_CORRECTION "Smooth transition between triangle walls" ON)
if(Mercury_TRIANGLE_WALL_CORRECTION)
	add_definitions( -DMercury_TRIANGLE_WALL_CORRECTION )
endif()

include(CMakeModules/MercuryMPI.cmake)

include(CMakeModules/MercuryOMP.cmake)

include(CMakeModules/MercuryLogger.cmake)

include(CMakeModules/MercuryBacktrace.cmake)

include(CMakeModules/MercuryCppFeatureCheck.cmake)

include(CMakeModules/MercuryDocumentation.cmake)

include(CMakeModules/MercuryXballs.cmake)

#This section generates a configuration file that contains the CMake variables#
#############################################
configure_file(Configuration/CMakeDefinitions.cc
		Kernel/CMakeDefinitions.cc)
configure_file(Configuration/CMakeDefinitions.h
		Kernel/CMakeDefinitions.h)

add_subdirectory(Kernel)
add_subdirectory(Tools)
add_subdirectory(XBalls)

# Determine whether to couple MercuryDPM with oomph-lib or not
option(Mercury_OOMPH_COUPLING "Couple MercuryDPM with oomph-lib" OFF)
# Couple with oomph-lib if OOMPH_COUPLING is on
include(CMakeModules/OomphCoupling.cmake)
# this is done after oomph is included since it needs the includes defined there
add_subdirectory(Drivers)

#Add the command fullTest, which rebuilds the cache, cleans the project and then compiles and tests everything.
add_custom_target(fullTest
    	COMMAND make rebuild_cache
	COMMAND make clean
	COMMAND make -j 4
	COMMAND make test -j 4
	COMMENT "This will clean the project, rebuild and test." )

#Add the command update, which rebuilds the cache and compiles the codes (no test).
add_custom_target(update
		COMMAND make rebuild_cache
		COMMAND make -j 4
		COMMENT "Quick rebuild" )

#This section is about the build of documentation#
##################################################

option(Mercury_BUILD_USER_DIR "Compile USER DIR in testing" ON)



###### These needs to be fixed as the new scripts line above will always create this, so the Python problem is back
################################################

FIND_PACKAGE(Python)
# this line causes problems when you have multiple builds
FILE(COPY Configuration/fpdiff.py DESTINATION ${PROJECT_SOURCE_DIR}/Scripts/)

#find_program(PYTHON_EXEC "python" 
#	DOC "Location of python executable to use")
