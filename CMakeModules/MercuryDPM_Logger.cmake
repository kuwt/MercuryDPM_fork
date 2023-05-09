#Configure the default loglevel using CMake..
set( MercuryDPM_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of MercuryDPM. DEFAULT is recommended.")
set_property( CACHE MercuryDPM_LOGLEVEL PROPERTY STRINGS FATAL ERROR WARN INFO DEFAULT VERBOSE DEBUG )
mark_as_advanced( FORCE MercuryDPM_LOGLEVEL )
add_definitions( -DMERCURYDPM_LOGLEVEL=Log::${MercuryDPM_LOGLEVEL} )