#Configure the default loglevel using CMake..
set(Mercury_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of Mercury. DEFAULT is recommended.")
set_property(CACHE Mercury_LOGLEVEL PROPERTY STRINGS FATAL ERROR WARN INFO DEFAULT VERBOSE DEBUG)
mark_as_advanced(FORCE Mercury_LOGLEVEL)
add_definitions(-DMERCURY_LOGLEVEL=Log::${Mercury_LOGLEVEL})

#DEPRECATED debated about this and now added a pre-compiler flag. But this part of code is useful to show how to add a user defined flag to Mercury.
#Configure a flag which enables/disables the manipulation of flushing behaviour in the logger.
# if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
# set( Mercury_ENABLE_FLUSHER OFF CACHE BOOL "Enable/disable manipulation of flushing behaviour in the
# logger.")
# else()
# set( Mercury_ENABLE_FLUSHER ON CACHE BOOL "Enable/disable manipulation of flushing behaviour in the logger
# ." )
# endif()
# mark_as_advanced( FORCE Mercury_ENABLE_FLUSHER )
# add_definitions( -DMERCURY_ENABLE_FLUSHER=$<BOOL:${Mercury_ENABLE_FLUSHER}> )

