cmake_minimum_required(VERSION 2.8)

include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

file(GLOB CPPFILES  "*.cpp")
#for every cpp found
foreach(CPPFILE ${CPPFILES})
    #extract the actually file name
    get_filename_component(FILENAME ${CPPFILE} NAME)
    #extract the filename minus the cpp. This will be the name of exe file
    get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
    #Make the exe
    add_executable(${EXECNAME} ${FILENAME})
    #All cpp folder and linked against DPMBase
    target_link_libraries(${EXECNAME}  MercuryBase)


file(COPY ${CMAKE_SOURCE_DIR}/Tools/CombineMultiTimeStepAllDataFiles.sh DESTINATION .)
endforeach()

