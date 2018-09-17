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
    #All cpp folder and linked agaist MercuryBase
	target_link_libraries(${EXECNAME}  MercuryBase)
endforeach()

#Part 2 : Copy example config file
#########################################
file(COPY test.config DESTINATION .)

#Part 3 : Copy all the tools
#########################################
#Collect all the bash and python scripts
file(GLOB TOOLS RELATIVE *.sh *.py)
#Copy and change permissions of the scripts
file(COPY ${TOOLS}
  DESTINATION .
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
#Copy gnuplot file and summary page
file(COPY visualize.gnuplot Summary.html DESTINATION .)
