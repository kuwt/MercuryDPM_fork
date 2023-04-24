#This part searches automatically for MPI files when enabled
option(MercuryDPM_USE_MPI "Include mpi this will enable parallel computation" OFF)

if(MercuryDPM_USE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    add_definitions( -DMERCURYDPM_USE_MPI )
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "The option you have chosen requires mpi and you do not have this installed. Please install")
    endif()
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()
