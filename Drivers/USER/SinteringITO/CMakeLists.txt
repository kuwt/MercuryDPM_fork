include_directories(${MercuryDPM_SOURCE_DIR}/Kernel
	                    ${MercuryDPM_BINARY_DIR}/Kernel)
add_subdirectory(ARS2017)


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

