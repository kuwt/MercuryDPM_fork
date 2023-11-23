# Creates a new option in cmake which includes the OpenMP libraries when enabled
OPTION (MercuryDPM_USE_OpenMP "Use OpenMP" OFF)

IF(MercuryDPM_USE_OpenMP)
    # Note: I tried to make OMP work on a mac using clang using the sources below, but I was unsuccessful.
    # - https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
    # - https://stackoverflow.com/questions/60005176/how-to-deal-with-clang-error-unsupported-option-fopenmp-on-travis
    # It works, however, with the gcc compiler. I used the following flags:
    # -DCMAKE_CXX_COMPILER=/opt/homebrew/Cellar/gcc/13.2.0/bin/g++-13 -DCMAKE_C_COMPILER=/opt/homebrew/Cellar/gcc/13.2.0/bin/gcc-13 -DMercuryDPM_USE_OpenMP=ON

    # Tries to find OpenMP; but does not require it.
    find_package(OpenMP)
    # Adds a precompiler flag, so you can use #ifdef MERCURYDPM_USE_OMP
    add_definitions( -DMERCURYDPM_USE_OMP)
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
