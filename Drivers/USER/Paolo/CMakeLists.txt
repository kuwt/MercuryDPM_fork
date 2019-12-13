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
	add_executable(${EXECNAME} ${FILENAME} ../../../Kernel/ClusterGenerator.h ../../../Kernel/ClusterGenerator.cc ../../../Kernel/BaseCluster.h ../../../Kernel/BaseCluster.cc
			BayesianCompression.h BayesianCompression.cc SingleClusterStressStrain.h SingleClusterStressStrain.cc SingleClusterCompression.h SingleClusterCompression.cc
			DoublePorosityStressStrain.h DoublePorosityStressStrain.cc DoublePorosityStressStrainA.h DoublePorosityStressStrainA.cc
			DoublePorosityStressStrainB.h DoublePorosityStressStrainB.cc DoublePorosityStressStrainC.h DoublePorosityStressStrainC.cc
			DoublePorosityStressStrainD.h DoublePorosityStressStrainD.cc BaseDoublePorosityStressStrain.h BaseDoublePorosityStressStrain.cc
			BaseDoublePorosityStressStrainD.h BaseDoublePorosityStressStrainD.cc)
	#All cpp folder and linked agaist DPMBase
	target_link_libraries(${EXECNAME}  MercuryBase)
endforeach()

