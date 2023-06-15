include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

#Required for Chute drivers
add_library(DieFilling STATIC
        DieFilling.cc
        )


#Part 1 : All cpp files will be made to exec files
#########################################
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
    target_link_libraries(${EXECNAME} oomphFluid MercuryBase DieFilling)
endforeach ()
