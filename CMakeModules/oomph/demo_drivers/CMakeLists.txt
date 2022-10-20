# subfolder solid contains a default CMakeList.txt; copy it other subfolders
FILE(GLOB ALL_DIRS *)
if (OOMPH_HAS_MPI)
    file(GLOB EXCLUDE_DIRS solid)
else()
    file(GLOB EXCLUDE_DIRS solid mpi)
endif()
list(REMOVE_ITEM ALL_DIRS ${EXCLUDE_DIRS})
FOREACH (DIR ${ALL_DIRS})
    IF (EXISTS "${DIR}/Makefile.am")
        FILE(COPY solid/CMakeLists.txt DESTINATION ${DIR})
    ENDIF ()
ENDFOREACH ()

#allow multiple targets with the same name
cmake_policy(SET CMP0002 OLD)

# add all subdirectories with a CMakeLists.txt file
# exception: the MPI directory
FILE(GLOB ALL_FILES *)
if (NOT OOMPH_HAS_MPI)
    file(GLOB EXCLUDE mpi)
    list(REMOVE_ITEM ALL_FILES ${EXCLUDE})
endif ()
FOREACH(DIR ${ALL_FILES})
    IF(EXISTS "${DIR}/CMakeLists.txt")
        add_subdirectory(${DIR})
    ENDIF()
ENDFOREACH()
