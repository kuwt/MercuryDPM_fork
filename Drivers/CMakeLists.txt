#add_subdirectory(Chute)
#add_subdirectory(Contraction)
add_subdirectory(SerializedDriver)
add_subdirectory(Leidenfrost)
add_subdirectory(GranularJet)
add_subdirectory(Avalanching)
add_subdirectory(CartesianShearCell)
add_subdirectory(Segregation)
add_subdirectory(UnitTests)
add_subdirectory(Clump)
add_subdirectory(MercurySimpleDemos)
add_subdirectory(MercuryCG)
add_subdirectory(MercuryMPITests)
add_subdirectory(ChuteDemos)
add_subdirectory(SuperQuadricDemos)
add_subdirectory(ImportTools)
add_subdirectory(Tutorials)
add_subdirectory(AdvancedTutorials)
add_subdirectory(CompareToStefan)
add_subdirectory(SelfTests)
add_subdirectory(Papers)
add_subdirectory(Sinter)
add_subdirectory(SpeedTest)
add_subdirectory(Demos)
add_subdirectory(Template)
add_subdirectory(ParallelDrum)
add_subdirectory(Clusters)
add_subdirectory(SimpleDrum)
add_subdirectory(Calibration)
add_subdirectory(TriangleWalls)
add_subdirectory(Oomph)
add_subdirectory(SimpleOpt)
add_subdirectory(RotatingDrum)




#IFCD: If you want (a folder of) the USER directory to be built, mark it so in your ccmake or cmake-gui

if (MercuryDPM_BUILD_USER_DIR)
    add_subdirectory(USER)
endif()
