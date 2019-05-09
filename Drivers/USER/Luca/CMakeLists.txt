include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)
add_subdirectory(cleanedAlgorithms)

file(GLOB CCFILES  "*.cc")


add_library(Luca STATIC
	ScrewAuger.cc
	NewInteractions/LinearPlasticViscoelasticInteractionCutoff.cc
	NewInteractions/LinearPlasticViscoelasticNormalSpeciesExtended.cc
	SphericalEnvelope.cc
	DualFeederCasing_BOT.cc
	DualFeederCasing_TOP.cc
	ScrewCC.cc
	CompressionPistonSurface.cc
	CompressionPiston.cc
	Helicoid.cc
	Helicoid02.cc
	Helicoid03.cc
	Helicoid04.cc
	Helicoid05.cc
	Helicoid06.cc
	OscillatingLattice1D.cc
	ScrewRectangularSection.cc
	ScrewRoundEdges.cc
	ScrewTriangularSection.cc
	HopperAndCasing.cc
	ScrewBottom.cc
	RotatingIntersectionOfWalls.cc
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
	#All cpp folder and linked agaist DPMBase
	target_link_libraries(Luca MercuryBase)
	target_link_libraries(${EXECNAME}  Luca)
endforeach()

