# the extern keyword is treated compiler-specific
#\todo add keywords for other fortran compilers
if(CMAKE_COMPILER_IS_GNUG77)
    add_definitions( -DgFortran )
endif()

# the preprocessor directives on a std. linux system
#\todo make generic, i.e. imitate the autogen behaviour
add_definitions( -DOOMPH_HAS_STACKTRACE )
add_definitions( -DOOMPH_HAS_UNISTDH -DOOMPH_HAS_MALLOCH)
add_definitions( -DOOMPH_HAS_TRIANGLE_LIB )
add_definitions( -DOOMPH_HAS_TETGEN_LIB -DTETLIBRARY )
add_definitions( -DUSING_OOMPH_SUPERLU )

#add CMakeCache options
option(OOMPH_HAS_HYPRE "Whether HYPRE is installed" OFF)
if (OOMPH_HAS_HYPRE)
    add_definitions( -DOOMPH_HAS_HYPRE )
endif()

#add CMakeCache options
option(OOMPH_HAS_TRILINOS "Whether TRILINOS is installed" OFF)
if (OOMPH_HAS_TRILINOS)
    add_definitions( -DOOMPH_HAS_TRILINOS )
endif()

#add CMakeCache option
option(OOMPH_HAS_MPI "Whether MPI is installed" OFF)
if (OOMPH_HAS_MPI)
    add_definitions( -DOOMPH_HAS_MPI )
    FIND_PACKAGE(MPI REQUIRED)
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "The option you have chosen requires mpi and you do not have this installed. Please install")
    endif()
    include_directories(${MPI_CXX_INCLUDE_PATH})
    add_definitions( -DUSING_OOMPH_SUPERLU_DIST )
endif()

# the #include statements in oomph-lib assume these paths to be set
include_directories(${OOMPH_DIR})
include_directories(${OOMPH_DIR}/src/generic)
include_directories(${OOMPH_DIR}/external_src)

# for #include "oomph_crbond_bessel.h,oomph_tetgen.h,fpu_control.h" to work
include_directories(${OOMPH_DIR}/external_src/oomph_crbond_bessel)
include_directories(${OOMPH_DIR}/external_src/oomph_tetgen)
include_directories(${OOMPH_DIR}/external_src/oomph_triangle)
include_directories(${OOMPH_DIR}/external_src/oomph_metis_from_parmetis_3.1.1)
if (OOMPH_HAS_MPI)
    include_directories(${OOMPH_DIR}/external_src/oomph_superlu_dist_3.0)
endif()

# for #include "../../../src/*" to work
include_directories(${OOMPH_DIR}/demo_drivers/solid/static_fish)

# first compile external_src, then src, then demo_drivers
add_subdirectory(external_src)
add_subdirectory(src)

option(OOMPH_COMPILE_DEMO_DRIVERS "Compile demo_drivers directory" OFF)
if (OOMPH_COMPILE_DEMO_DRIVERS)
    add_subdirectory(demo_drivers)
endif()