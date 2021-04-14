include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

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
	target_link_libraries(${EXECNAME}  MercuryBase)
endforeach()

#add_custom_target(POSTPROCESS ALL
#	COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/PostProcess ${CMAKE_CURRENT_BINARY_DIR}
#	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/Matlab/thomas/loadstatistics.m ${CMAKE_CURRENT_BINARY_DIR}
#COMMENT "Copying post-processing files to build directory"
#	COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/Drivers/MercuryCG/fstatistics ${CMAKE_CURRENT_BINARY_DIR}
#)

