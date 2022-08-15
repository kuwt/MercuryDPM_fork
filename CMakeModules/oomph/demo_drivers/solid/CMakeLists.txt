# This file is a copy of the CMakeLists.txt template in demo_drivers/solid; only change the original file in that folder
include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

# extract name of current directory
get_filename_component(folder ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# add necessary include statements for driver codes
FILE(GLOB ALL_FILES ${OOMPH_DIR}/src/*)
FOREACH(DIR ${ALL_FILES})
	IF(EXISTS "${DIR}/CMakeLists.txt")
		include_directories(${DIR})
	ENDIF()
ENDFOREACH()
include_directories(${OOMPH_DIR}/src/)
include_directories(${OOMPH_DIR}/external_src/)

#All source files will be made to exec files
file(GLOB CPPFILES *.cc)
foreach(CPPFILE ${CPPFILES})
	#extract the actually file name
	get_filename_component(FILENAME ${CPPFILE} NAME)
	#extract the filename minus the cpp. This will be the name of exe file
	get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
	#Make the executable; first the exceptions, then the default behaviour

	# executables with several cc files
	if((EXECNAME STREQUAL "refineable_simple_shear")
		OR (EXECNAME STREQUAL "simple_shear"))
		add_executable(${EXECNAME} ${FILENAME} precompiled_mesh.cc)
	elseif(EXECNAME STREQUAL "two_d_poisson2")
		#exception: refineable_simple_shear uses a precompiled mesh
		add_executable(${EXECNAME} ${FILENAME} two_d_poisson2_mesh.cc)
	elseif(EXECNAME STREQUAL "two_dim_hex")
		#exception: refineable_simple_shear uses a precompiled mesh
		add_executable(${EXECNAME} ${FILENAME} homo_lin_elasticity_elements.cc Thomo_lin_elasticity_elements.cc)
	elseif(EXECNAME STREQUAL "two_dim")
		#exception: refineable_simple_shear uses a precompiled mesh
		add_executable(${EXECNAME} ${FILENAME} homo_lin_elasticity_elements.cc Thomo_lin_elasticity_elements.cc)

	#executables with conditional compilation
	elseif((NOT OOMPH_HAS_HYPRE) AND (EXECNAME STREQUAL "HypreSolver_test"))
		#message("skipping ${EXECNAME} b/c hypre not installed")
		continue()
	elseif((NOT OOMPH_HAS_TRILINOS) AND ((EXECNAME STREQUAL "TrilinosSolver_test") OR (EXECNAME STREQUAL "two_d_multi_poisson") OR (EXECNAME STREQUAL "harmonic")))
		#message("skipping ${EXECNAME} b/c trilinos not installed")
		continue()
	elseif(((NOT OOMPH_HAS_HYPRE) OR (NOT OOMPH_HAS_TRILINOS)) AND (EXECNAME STREQUAL "pseudo_solid_collapsible_tube"))
		#message("skipping ${EXECNAME} b/c hypre or trilinos not installed")
		continue()

	#executables excluded in teh original oomph-lib
	elseif((EXECNAME STREQUAL "precompiled_mesh")
			OR (EXECNAME STREQUAL "two_d_poisson2_mesh")
			OR (EXECNAME STREQUAL "refineable_linearised_axisym_navier_stokes_elements")
			OR (EXECNAME STREQUAL "linearised_axisym_navier_stokes_elements")
			OR (EXECNAME STREQUAL "homo_lin_elasticity_elements")
			OR (EXECNAME STREQUAL "Thomo_lin_elasticity_elements")
			OR (EXECNAME STREQUAL "two_dim_hex")
			)
		#message("skipping ${EXECNAME} b/c no main")
		continue()
	elseif((EXECNAME STREQUAL "driven_cavity_with_simple_lsc_preconditioner")
			OR (EXECNAME STREQUAL "linearised_axisym_navier_stokes_elements")
			OR (EXECNAME STREQUAL "linearised_axisymmetric_fluid_interface_elements")
			OR (EXECNAME STREQUAL "perturbed_spines")
			OR (EXECNAME STREQUAL "two_layer_perturbed_spine_mesh")
			)
		#message("skipping ${EXECNAME} b/c in extra_dist")
		continue()
	elseif(EXECNAME STREQUAL "marangoni_convection")
		#message("skipping ${EXECNAME} b/c [no idea whyit's excluded in oomph-lib]")
		continue()
	elseif((EXECNAME STREQUAL "steady_axisym_advection_diffusion")
			OR (EXECNAME STREQUAL "young_laplace")
			OR (EXECNAME STREQUAL "inclined_plane")
			OR (EXECNAME STREQUAL "counter_rotating_disks")
			)
		#message("executable for ${CPPFILE} set to ${EXECNAME}_2 b/c name conflict")
		set(EXECNAME "${EXECNAME}_2")
		add_executable(${EXECNAME} ${FILENAME})
	elseif(${CMAKE_CURRENT_SOURCE_DIR} MATCHES "^.*/demo_drivers/mpi/*")
		#message("executable for ${CPPFILE} set to ${EXECNAME}_2 b/c name conflict")
		set(EXECNAME "${EXECNAME}_MPI")
		add_executable(${EXECNAME} ${FILENAME})
	#default behaviour
	else()
		add_executable(${EXECNAME} ${FILENAME})
	endif()

    #link executable against all src and external_src libraries (except axisym_displ_based_foeppl_von_karman pml_fourier_decomposed_helmholtz pml_helmholtz pml_time_harmonic_linear_elasticity)
	#exceptions: folders that are in conflict with the standard oomph libraries
	if (${folder} STREQUAL "axisym_displ_based_foeppl_von_karman")
		target_link_libraries(${EXECNAME} axisym_displ_based_foeppl_von_karman generic)
	elseif(${CMAKE_CURRENT_SOURCE_DIR} MATCHES "^.*/pml_fourier_decomposed_helmholtz/oscillating_sphere")
		target_link_libraries(${EXECNAME} pml_fourier_decomposed_helmholtz generic)
	elseif(${CMAKE_CURRENT_SOURCE_DIR} MATCHES "^.*/pml_helmholtz/scattering")
		target_link_libraries(${EXECNAME} pml_helmholtz meshes generic)
	elseif(${CMAKE_CURRENT_SOURCE_DIR} MATCHES "^.*/pml_time_harmonic_linear_elasticity")
		target_link_libraries(${EXECNAME} pml_time_harmonic_linear_elasticity  generic)
	else()
		target_link_libraries(${EXECNAME} oomph)
	endif()
endforeach()

#add all subdirectories with a CMakeLists.txt file
FILE(GLOB ALL_FILES *)
FOREACH(DIR ${ALL_FILES})
	IF(EXISTS "${DIR}/Makefile.am")
		FILE(COPY CMakeLists.txt DESTINATION ${DIR})
		add_subdirectory(${DIR})
	ENDIF()
ENDFOREACH()

#copy certain files to the build
FILE(GLOB INPUT_FILES  "*.ele" "*.node" "*.face" "*.poly")
FILE(COPY ${INPUT_FILES} DESTINATION .)

