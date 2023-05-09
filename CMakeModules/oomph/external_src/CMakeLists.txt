# src/solid contains the default CMakeList.txt; copy it into all external_src folders
# exceptions: oomph_metis_from_parmetis_3.1.1 oomph_parmetis_3.1.1 oomph_superlu_dist_3.0
FOREACH(DIR oomph_arpack oomph_blas oomph_crbond_bessel oomph_hsl oomph_lapack oomph_flapack oomph_superlu_4.3 oomph_tetgen oomph_triangle oomph_metis_from_parmetis_3.1.1)
    FILE(COPY ../src/solid/CMakeLists.txt DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/${DIR})
ENDFOREACH()
if (OOMPH_HAS_MPI)
    FILE(COPY ../src/solid/CMakeLists.txt DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/oomph_superlu_dist_3.0)
else()
    FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/oomph_superlu_dist_3.0/CMakeLists.txt)
endif()

# add all subdirectories with a CMakeLists file
FILE(GLOB ALL_FILES *)
FOREACH(DIR ${ALL_FILES})
    IF(EXISTS "${DIR}/CMakeLists.txt")
        add_subdirectory(${DIR})
    ENDIF()
ENDFOREACH()
