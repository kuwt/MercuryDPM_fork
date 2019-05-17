include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

find_package(Boost 1.60.0)
if(NOT Boost_FOUND)
	message(FATAL_ERROR "Could not find boost!")
endif()
#set(Boost_INCLUDE_DIR c:/mingw/include/)
#set(Boost_LIBRARY_DIR c:/mingw/lib/)
#find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

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
	#All cpp folder and linked agaist DPMBase
	target_link_libraries(${EXECNAME}  MercuryBase)
endforeach()

