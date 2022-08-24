if(NOT OOMPH_COUPLING)
  #message(STATUS "Coupling with oomph-lib is disabled")
  return()
endif()

include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)
include_directories(${Mercury_SOURCE_DIR}/oomph-lib/src/generic)
include_directories(${Mercury_SOURCE_DIR}/oomph-lib/src/solid)
include_directories(${Mercury_SOURCE_DIR}/oomph-lib/src/constitutive)
include_directories(${Mercury_SOURCE_DIR}/oomph-lib/src/)
include_directories(${Mercury_SOURCE_DIR}/oomph-lib/external_src/)

FILE(GLOB ALL_FILES *)
FOREACH(FILE ${ALL_FILES})
  IF(EXISTS "${FILE}/CMakeLists.txt")
    #MESSAGE(STATUS "add_sub(${FILE})")
    add_subdirectory(${FILE})
  ENDIF()
ENDFOREACH()

#All cpp files will be made to exec files
file(GLOB CPPFILES "*.cpp")
#for every cpp found
foreach (CPPFILE ${CPPFILES})
  #extract the actually file name
  get_filename_component(FILENAME ${CPPFILE} NAME)
  #extract the filename minus the cpp. This will be the name of exe file
  get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
  #Make the exe
  add_executable(${EXECNAME} ${FILENAME})
  #All cpp folder and linked agaist MercuryBase
  target_link_libraries(${EXECNAME} oomphSolid MercuryBase)
endforeach ()
