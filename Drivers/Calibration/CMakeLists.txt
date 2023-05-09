include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)
set(CMAKE_CXX_STANDARD 17)

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
	#All cpp folder and linked against Chute
	target_link_libraries(${EXECNAME}  MercuryBase)
endforeach()



#Calibration directory is in the source; this adds a link to the build, so the executables can be built
if(NOT EXISTS build)
  configure_file(py/build.in ${PROJECT_SOURCE_DIR}/Drivers/Calibration/py/build)
endif()