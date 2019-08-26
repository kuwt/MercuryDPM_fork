# add all subdirectories
FILE(GLOB ALL_FILES *)
FOREACH(FILE ${ALL_FILES})
    IF(EXISTS "${FILE}/CMakeLists.txt")
        add_subdirectory(${FILE})
    ENDIF()
ENDFOREACH()
