#Configure the default loglevel using CMake..
set( Mercury_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of Mercury. DEFAULT is recommended.")
set_property( CACHE Mercury_LOGLEVEL PROPERTY STRINGS FATAL ERROR WARN INFO DEFAULT VERBOSE DEBUG )
mark_as_advanced( FORCE Mercury_LOGLEVEL )
add_definitions( -DMERCURY_LOGLEVEL=Log::${Mercury_LOGLEVEL} )