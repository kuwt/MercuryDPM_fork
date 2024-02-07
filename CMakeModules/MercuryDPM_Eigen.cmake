#Do nothing if coupling is turned off.
if(NOT MercuryDPM_EIGEN)
    #message(STATUS "Eigen is disabled")
    return()
endif()

# Clone Eigen if has not been cloned before
set(EIGEN_DIR ${PROJECT_SOURCE_DIR}/Eigen)
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND git submodule init ${EIGEN_DIR})
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMAND git submodule update ${EIGEN_DIR})
if(EXISTS ${EIGEN_DIR}/demos)
    message(STATUS "Eigen was sucessfully initialised and cloned")
else()
    message(FATAL_ERROR "git clone eigen failed. If this problem persists you can manually clone Eigen by running: \n   git submodule init\n   git submodule update")
endif()

include_directories(${PROJECT_SOURCE_DIR}/Eigen/)
