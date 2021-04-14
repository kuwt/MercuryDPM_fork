add_subdirectory(Math)

#Part 1 : All cpp files will be made to exec files
#########################################
#Collect all the names of the cpp, note at the moment limited to demos, but will be fixed later
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
    target_link_libraries(${EXECNAME} MercuryBase)
endforeach()

# JMFT: This line needs to come at the bottom, after Parts 1 and 2.
# This is the other way round to the CMakeLists.txt file in other directories,
# where it has to come at the top;
# see e.g. ../Tutorials/CMakeLists.txt.
# This is because the tests in this directory will produce executables whose
# names end in 'demo', which need to be run _before_ outputs from them are
# compared.
# TODO: This naming should be changed to avoid this kludge.
include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)


