# Src/solid contains a default CMakeList.txt; copy it into all other folders
# exception: generic has a special CMakeList.txt
FILE(GLOB ALL_DIRS *)
file(GLOB EXCLUDE_DIRS solid generic)
list(REMOVE_ITEM ALL_DIRS ${EXCLUDE_DIRS})
FOREACH (DIR ${ALL_DIRS})
    IF (EXISTS "${DIR}/Makefile.am")
        FILE(COPY solid/CMakeLists.txt DESTINATION ${DIR})
    ENDIF ()
ENDFOREACH ()

# add all subdirectories with a CMakeLists file
FILE(GLOB ALL_DIRS *)
FOREACH (DIR ${ALL_DIRS})
    IF (EXISTS "${DIR}/CMakeLists.txt")
        add_subdirectory(${DIR})
        get_filename_component(NAME ${DIR} NAME)
    ENDIF ()
ENDFOREACH ()

## create an empty library named oomph (I used one dummy source file as I can't make an empty library)
#add_library(oomph STATIC generic/oomph_definitions.cc)
#
## link all Src libraries to oomph
## exception: axisym_displ_based_foeppl_von_karman pml_fourier_decomposed_helmholtz pml_helmholtz pml_time_harmonic_linear_elasticity are in conflict with other parts of oomph (same class names)
## note: generic is already linked to all the external_src libraries
#target_link_libraries(oomph steady_axisym_advection_diffusion young_laplace advection_diffusion advection_diffusion_reaction axisym_advection_diffusion axisym_foeppl_von_karman axisym_linear_elasticity axisym_navier_stokes axisym_poroelasticity axisym_spherical_solid beam biharmonic constitutive darcy fluid_interface flux_transport foeppl_von_karman fourier_decomposed_helmholtz generalised_newtonian_axisym_navier_stokes generalised_newtonian_navier_stokes helmholtz linear_elasticity linear_wave linearised_axisym_navier_stokes mesh_smoothing meshes multi_physics navier_stokes ode  poisson polar_navier_stokes poroelasticity rigid_body shell solid spherical_advection_diffusion spherical_navier_stokes  time_harmonic_fourier_decomposed_linear_elasticity time_harmonic_linear_elasticity unsteady_heat womersley generic)
#
