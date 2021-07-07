#Setup scripts directory
########################################################

file(GLOB SCRIPTFILES "*.sh")
foreach(SCRIPTFILE ${SCRIPTFILES})
    get_filename_component(FILENAME ${SCRIPTFILE} NAME)
    get_filename_component(EXECNAME ${SCRIPTFILE} NAME_WE)
    add_custom_target(${FILENAME}
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
            COMMAND chmod 700 ${FILENAME}
            #COMMAND ./${FILENAME}
            )
endforeach()

#First create the directory if does not exist this stops the dir does not exist error in the next line
#add_custom_command(TARGET ${EXECNAME} PRE_BUILD
#        COMMAND ${CMAKE_COMMAND} -E make_directory  ${CMAKE_CURRENT_BINARY_DIR}/Scripts/)


#Copy and copy back
#add_custom_command(TARGET ${EXECNAME} POST_BUILD
        #COMMAND ${CMAKE_COMMAND} -E copy_directory
        #${CMAKE_CURRENT_BINARY_DIR}/Scripts/ ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/
#        COMMAND ${CMAKE_COMMAND} -E copy_directory
#        ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/ ${CMAKE_CURRENT_BINARY_DIR}/Scripts/)