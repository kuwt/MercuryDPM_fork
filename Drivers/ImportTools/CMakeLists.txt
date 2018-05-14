include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

add_executable(BinaryReaderSTLOneTriangleUnitTest BinaryReaderSTLOneTriangleUnitTest.cpp ${Mercury_SOURCE_DIR}/Kernel/BinaryReader.cc ${Mercury_SOURCE_DIR}/Kernel/Math/Vector.cc ${Mercury_SOURCE_DIR}/Kernel/Logger.cc ${Mercury_SOURCE_DIR}/Kernel/CMakeDefinitions.cc)

add_executable(BinaryReaderSTL1by1by1bySquareUnitTest BinaryReaderSTL1by1by1bySquareUnitTest.cpp ${Mercury_SOURCE_DIR}/Kernel/BinaryReader.cc ${Mercury_SOURCE_DIR}/Kernel/Math/Vector.cc ${Mercury_SOURCE_DIR}/Kernel/Logger.cc ${Mercury_SOURCE_DIR}/Kernel/CMakeDefinitions.cc)

