include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

add_library(Screw STATIC
		ScrewRectangularSectionClean.cc
		ScrewRoundEdgesClean.cc
		ScrewTriangularSectionClean.cc
		../Helicoid05.cc
)


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
	target_link_libraries(Screw MercuryBase)
	target_link_libraries(${EXECNAME}  Screw)
endforeach()
