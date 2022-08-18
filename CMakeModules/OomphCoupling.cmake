#Do nothing if coupling is turned off.
if(NOT OOMPH_COUPLING)
    #message(STATUS "Coupling with oomph-lib is disabled")
    return()
endif()

# Clone oomph-lib if has not been cloned before
set(OOMPH_DIR ${PROJECT_SOURCE_DIR}/oomph-lib)
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/oomph-lib/src)
    message(STATUS "Cloning https://github.com/oomph-lib/oomph-lib.git")
    execute_process(COMMAND git clone https://github.com/oomph-lib/oomph-lib.git ${OOMPH_DIR})
#else()
    #message(STATUS "Oomph-lib is residing here: ${OOMPH_DIR}")
endif()

option(OOMPH_CMAKE "Use cmake to couple oomph-lib (if coupling turned on)" ON)
if (OOMPH_CMAKE)
    message(STATUS "Compiling oomph with CMake")
    # copy CMakeFiles
    FILE(COPY CMakeModules/oomph/CMakeLists.txt DESTINATION ${OOMPH_DIR})
    FILE(COPY CMakeModules/oomph/src/CMakeLists.txt DESTINATION ${OOMPH_DIR}/src)
    FILE(COPY CMakeModules/oomph/src/generic/CMakeLists.txt DESTINATION ${OOMPH_DIR}/src/generic)
    FILE(COPY CMakeModules/oomph/src/solid/CMakeLists.txt DESTINATION ${OOMPH_DIR}/src/solid)
    FILE(COPY CMakeModules/oomph/external_src/CMakeLists.txt DESTINATION ${OOMPH_DIR}/external_src)
    FILE(COPY CMakeModules/oomph/demo_drivers/CMakeLists.txt DESTINATION ${OOMPH_DIR}/demo_drivers)
    FILE(COPY CMakeModules/oomph/demo_drivers/solid/CMakeLists.txt DESTINATION ${OOMPH_DIR}/demo_drivers/solid)
    FILE(COPY CMakeModules/oomph/demo_drivers/solid/CMakeLists.txt DESTINATION ${OOMPH_DIR}/demo_drivers/solid)
    if (${APPLE})
        message (STATUS "MAC OS X: making dummy fpu_control.h")
        IF(NOT EXISTS ${OOMPH_DIR}/external_src/oomph_triangle/fpu_control.h)
            FILE(COPY ${OOMPH_DIR}/external_src/oomph_triangle/dummy_fpu_control.h DESTINATION ${OOMPH_DIR})
            FILE(RENAME ${OOMPH_DIR}/dummy_fpu_control.h ${OOMPH_DIR}/external_src/oomph_triangle/fpu_control.h)
        ENDIF()
        #-DCMAKE_C_FLAGS=-Wno-implicit-function-declaration -DCMAKE_CXX_FLAGS=-Wno-undefined-var-template -DOOMPH_CMAKE=ON -DCMAKE_OSX_DEPLOYMENT_TARGET=11.3
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-function-declaration")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
    endif ()
    add_subdirectory(${OOMPH_DIR})
else()
    # Otherwise, we assume it is already compiled
    message(STATUS "Coupling to external oomph directory" ${OOMPH_DIR})
    #check if the oomph-lib src directory exists
    if (NOT IS_DIRECTORY ${OOMPH_DIR}/src)
        message(FATAL_ERROR "${OOMPH_DIR}/src does not exist.\n $ENV{HOME} Set OOMPH_DIR to the directory where the src/ and build/lib/ folders of oomphlib reside")
    endif()
    #check if the oomph-lib standard libraries have been built
    if (NOT IS_DIRECTORY ${OOMPH_DIR}/build/lib)
        message(FATAL_ERROR "${OOMPH_DIR}/build/lib does not exist.\n Please compile oomph-lib first")
    endif()
    # include the folder where the oomph-lib libraries reside
    # (so we can link against the libraries, see e.g. Drivers/OomphlibCoupling/CMakeLists.txt)
    link_directories(${OOMPH_DIR}/build/lib)
endif()

# include the folders where the oomph-lib src files reside
# (so you can write e.g. #include "mesh.h")
include_directories(${OOMPH_DIR}/src ${OOMPH_DIR}/src/poisson ${OOMPH_DIR}/src/generic ${OOMPH_DIR}/src/solid ${OOMPH_DIR}/src/constitutive ${OOMPH_DIR}/external_src)

add_library(oomph STATIC ${CMAKE_SOURCE_DIR}/Kernel/Logger.cc)

if (OOMPH_HAS_MPI)
    set(MPILIBS oomph_superlu_dist_3.0)
    add_definitions( -DOOMPH_HAS_MPI )
endif()

target_link_libraries(oomph steady_axisym_advection_diffusion young_laplace advection_diffusion advection_diffusion_reaction axisym_advection_diffusion axisym_foeppl_von_karman axisym_linear_elasticity axisym_navier_stokes axisym_poroelasticity axisym_spherical_solid beam biharmonic constitutive darcy fluid_interface flux_transport foeppl_von_karman fourier_decomposed_helmholtz generalised_newtonian_axisym_navier_stokes generalised_newtonian_navier_stokes helmholtz linear_elasticity linear_wave linearised_navier_stokes linearised_axisym_navier_stokes mesh_smoothing meshes multi_physics navier_stokes ode poisson polar_navier_stokes poroelasticity rigid_body shell solid spherical_advection_diffusion spherical_navier_stokes  time_harmonic_fourier_decomposed_linear_elasticity time_harmonic_linear_elasticity unsteady_heat womersley generic oomph_hsl oomph_arpack oomph_crbond_bessel oomph_triangle oomph_tetgen oomph_superlu_4.3 ${MPILIBS} oomph_lapack oomph_flapack oomph_blas  oomph_metis_from_parmetis_3.1.1)
# reynolds_averaged_navier_stokes