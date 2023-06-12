file(GLOB ALL_FILES *)

foreach(FILE ${ALL_FILES})
  if(EXISTS "${FILE}/CMakeLists.txt")
    add_subdirectory(${FILE})
  endif()
endforeach()
