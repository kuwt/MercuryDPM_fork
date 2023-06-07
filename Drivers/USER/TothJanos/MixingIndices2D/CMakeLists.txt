include(${CMAKE_SOURCE_DIR}/Configuration/MasterMake.cmake)

file(GLOB CPPFILES  "*.cpp")
foreach(CPPFILE ${CPPFILES})
    get_filename_component(FILENAME ${CPPFILE} NAME)
    get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
    add_executable(${EXECNAME} ${FILENAME})
    target_link_libraries(${EXECNAME} MercuryBase)
endforeach()
