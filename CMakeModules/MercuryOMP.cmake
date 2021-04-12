# Creates a new option in cmake which includes the OpenMP libraries when enabled
OPTION (Mercury_USE_OpenMP "Use OpenMP" OFF)

IF(Mercury_USE_OpenMP)
    # Tries to find OpenMP; but does not require it.
    find_package(OpenMP)
    # Adds a precompiler flag, so you can use #ifdef MERCURY_USE_OMP
    add_definitions( -DMERCURY_USE_OMP)
    # If OpenMP cannot be found, we look for a local install.
    # If none is found, an error message is thrown.
    IF(NOT OPENMP_FOUND)
        if (EXISTS /usr/local/opt/libomp/)
            message(STATUS "Found OpenMP in /usr/local/opt/libomp/.")
            SET(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -lomp -I/usr/local/opt/libomp/include -L/usr/local/opt/libomp/lib")
        else()
            message(FATAL_ERROR "The option you have chosen requires openmp and you do not have this installed. Please install")
        endif()
    ENDIF()
    # Set compiler flags
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()
