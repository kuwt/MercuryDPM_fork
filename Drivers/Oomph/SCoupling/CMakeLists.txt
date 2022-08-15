include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)
include_directories(${Mercury_SOURCE_DIR}/Oomph-lib/src/generic)
include_directories(${Mercury_SOURCE_DIR}/Oomph-lib/src/solid)
include_directories(${Mercury_SOURCE_DIR}/Oomph-lib/src/constitutive)
include_directories(${Mercury_SOURCE_DIR}/Oomph-lib/src/)
include_directories(${Mercury_SOURCE_DIR}/Oomph-lib/external_src/)

#Part 1 : All cpp files will be made to exec files
#########################################
#Collect all the names of the cpp, note at the moment limited to demos, but will be fixed later
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
    target_link_libraries(${EXECNAME} oomph MercuryBase)
endforeach ()

##Adds extra command to the build of SolidBagSelfTest
#add_custom_command(TARGET SolidBagSelfTest
#        COMMAND make -j SolidBag1_3SS_0 SolidBag1_3SS_1 SolidBag1_3SS_2 SolidBag1 SolidBag2 SolidBag2b SolidBag3 SolidBag4
#        COMMENT "SolidBagSelfTest depends on multiple SolidBag* executables, so needs to compile them as well" )

#FILE(GLOB INPUT_FILES  "*.ele" "*.node" "*.face" "*.poly" "*.txt" "Stage*" "*.in")
#FILE(COPY ${INPUT_FILES} DESTINATION .)
