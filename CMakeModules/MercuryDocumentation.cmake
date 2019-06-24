cmake_policy(SET CMP0057 NEW)
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
    option(Mercury_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" ON)
else()
    message(STATUS "Doxygen found: ${DOXYGEN_FOUND}, dot found: ${DOXYGEN_DOT_FOUND}: Local documentation will not be built")
    option(Mercury_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
endif()


if(Mercury_BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen REQUIRED dot)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR
                "Doxygen is needed to build the documentation. Please install it correctly or turn off Mercury_BUILD_DOCUMENTATION")
    else()
        #This is the configure file for normal doxygen builds
        configure_file(Configuration/doxygen.conf
                ${PROJECT_BINARY_DIR}/Configuration/doxygen.conf  @ONLY IMMEDIATE)

        #The next four and for website doxygen builds. The should be hinded from the public cmake at some point
        configure_file(Configuration/web_doxygen.conf
                ${PROJECT_BINARY_DIR}/Configuration/web_doxygen.conf @ONLY IMMEDIATE)

        configure_file(Configuration/hpg.css
                ${PROJECT_BINARY_DIR}/Configuration/hpg.css @ONLY IMMEDIATE)

        configure_file(Configuration/new_footer.html
                ${PROJECT_BINARY_DIR}/Configuration/new_footer.html @ONLY IMMEDIATE)

        configure_file(Configuration/new_header.html
                ${PROJECT_BINARY_DIR}/Configuration/new_header.html @ONLY IMMEDIATE)

        #-- Add custom targets to both make and clean (delete) the documentation
        add_custom_target 	(doc
                COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Configuration/doxygen.conf
                SOURCES ${PROJECT_BINARY_DIR}/Configuration/doxygen.conf)

        add_custom_target	(docClean
                COMMAND mkdir -p ${PROJECT_BINARY_DIR}/Documentation
                COMMAND rm -r ${PROJECT_BINARY_DIR}/Documentation/*
                COMMENT "Cleaning (deleting) the documentation"	)

        add_custom_target	(docWeb
                COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Configuration/web_doxygen.conf
                SOURCES ${PROJECT_BINARY_DIR}/Configuration/web_doxygen.conf)
    endif()
endif()
