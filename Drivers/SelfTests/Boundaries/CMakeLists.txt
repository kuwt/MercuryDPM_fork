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
    #All cpp folder and linked against MercuryBase
	target_link_libraries(${EXECNAME} Chute)
endforeach()

#Copy specific files to your build directory
#file(COPY ${CMAKE_SOURCE_DIR}/Drivers/MercurySimpleDemos/TriangulatedWall.vtk DESTINATION .)

#copy restart input for self tests into build, deprecated because we don't want to copy the data over
#file(COPY ${CMAKE_SOURCE_DIR}/Drivers/SelfTests/Boundaries/InputData/StressStrainControlBoundarySelfTestInput.restart DESTINATION .)
