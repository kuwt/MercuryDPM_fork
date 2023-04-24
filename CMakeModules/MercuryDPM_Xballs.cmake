if (WIN32)
    message("Windows operating system: Turning off xballs support")
    option(MercuryDPM_INCLUDE_XBALLS_SUPPORT "Use the xballs of Stefan Luding for visualisation" OFF)
else(WIN32)
    option(MercuryDPM_INCLUDE_XBALLS_SUPPORT "Use the xballs of Stefan Luding for visualisation" OFF)
endif(WIN32)

if (MercuryDPM_INCLUDE_XBALLS_SUPPORT)
	FIND_PACKAGE(X11)
        configure_file(${PROJECT_SOURCE_DIR}/XBalls/xballs.txt ${PROJECT_BINARY_DIR}/XBalls/xballs.txt)
	 configure_file(${PROJECT_SOURCE_DIR}/XBalls/MakeMovie
	    		${PROJECT_BINARY_DIR}/XBalls/MakeMovie @ONLY IMMEDIATE)

    if (NOT X11_FOUND)
        message(FATAL_ERROR "X11 is needed, please install it")
    endif()
endif()
