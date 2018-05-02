include_directories(${Mercury_SOURCE_DIR}/Kernel
                    ${Mercury_BINARY_DIR}/Kernel)

#Define a function to extract the number of required processors
function(get_number_of_cores EXECNAME NUMCORES)
  string(FIND ${EXECNAME} "MPI" POS1)
  string(FIND ${EXECNAME} "Test" POS2)
  math(EXPR START "${POS1} + 3")
  math(EXPR LENGTH "${POS2} - ${START}") 
  if (${LENGTH} STREQUAL 0)
    message(FATAL_ERROR "No number of cores specified for ${EXECNAME}. Format is *MPI<number_of_cores>Test.cpp")
  endif()
  string(SUBSTRING ${EXECNAME} ${START} ${LENGTH} NUMCORES )
  set(NUMCORES ${NUMCORES} PARENT_SCOPE)
endfunction()

#Part 2 : Make run test for each of the demo files
##################################################

file(GLOB SELFTESTS "*SelfTest.cpp")
file(GLOB UNITTESTS "*UnitTest.cpp")
file(GLOB MPITESTS  "*MPI*Test.cpp")
#for each demo add a test with the same name
if (Mercury_USE_MPI)
	foreach (TEST ${UNITTESTS} ${SELFTESTS})
        	get_filename_component(EXECNAME ${TEST} NAME_WE)
        	add_test(${EXECNAME} ${EXECNAME})
	endforeach()
	foreach (TEST ${MPITESTS})
        	get_filename_component(EXECNAME ${TEST} NAME_WE)
            get_number_of_cores(${EXECNAME} NUMCORES)
        	add_test(${EXECNAME} mpiexec -n ${NUMCORES} ./${EXECNAME})
	endforeach()

else()
	foreach (TEST ${UNITTESTS} ${SELFTESTS})
		get_filename_component(EXECNAME ${TEST} NAME_WE)
		add_test(${EXECNAME} ${EXECNAME} )
	endforeach()
endif()

#Part 3 : Make tests for each of the selftest_data files
########################################################


file(GLOB TESTDATAFILES "${CMAKE_CURRENT_SOURCE_DIR}/SelfTestData/*.*")
#for each file in the selftest_data folder create a test. Which checks the data against this old data. The actually testing is done my the script self_test.
foreach(TESTFILE ${TESTDATAFILES})
        get_filename_component(TESTNAME ${TESTFILE} NAME)
	add_test(${TESTNAME} ${CMAKE_SOURCE_DIR}/Scripts/self_test ${TESTFILE} ${TESTNAME})
        #Add the newly created files to the clean target
        set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${TESTNAME}")
endforeach()
