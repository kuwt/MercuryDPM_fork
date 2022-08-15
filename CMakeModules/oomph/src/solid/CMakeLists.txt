# This file is a copy of the CMakeLists.txt template in src/solid;
# only change the original file in that folder

# extract name of current directory
get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# write all header file names into a generic header $dir.h
# get all header files
file(GLOB headers *.h)
# remove common header from header list
file(GLOB commonHeader ${dir}.h)
if (commonHeader)
    list(REMOVE_ITEM headers ${commonHeader})
endif()

# get all source files
file(GLOB sources *.cc *.f *.c *.cxx)

# exceptions for oomph_superlu_dist_3.0
if (${CMAKE_CURRENT_SOURCE_DIR} MATCHES "^.*/oomph_superlu_dist_3.0")
    set(headers superlu_defs.h  superlu_ddefs.h Cnames.h machines.h psymbfact.h supermatrix.h old_colamd.h   util_dist.h superlu_enum_consts.h)
    set(sources sp_ienv.c etree.c sp_colorder.c get_perm_c.c mmd.c comm.c memory.c util.c superlu_grid.c pxerbla.c superlu_timer.c GetDiagU.c mc64ad.c symbfact.c psymbfact.c psymbfact_util.c get_perm_c_parmetis.c dlangs.c dgsequ.c dlaqgs.c dutil.c dmemory.c dmyblas2.c dsp_blas2.c dsp_blas3.c dldperm.c ddistribute.c pdgstrf.c pdgstrs_Bglobal.c pdgstrs1.c pdgssvx_ABglobal.c pdgsrfs_ABXglobal.c pdgsmv_AXglobal.c pdgssvx.c pdgstrs.c pddistribute.c pdlangs.c pdutil.c pdgsequ.c pdlaqgs.c pdgsrfs.c pdgsmv.c pdgstrs_lsum.c pdsymbfact_distdata.c dreadhb.c dreadrb.c dreadtriple.c)
endif()

# create common header (unless no header files are found)
if (headers)
    foreach(header ${headers})
        get_filename_component(name ${header} NAME)
        set(content "${content}#include \"${name}\"\n")
    endforeach()
    file(WRITE ${dir}.h ${content})
endif()

# exclude some source files
file(GLOB exclude_sources Taxisym_linear_elasticity_elements.cc)
if (NOT (exclude_sources STREQUAL ""))
    list(REMOVE_ITEM sources ${exclude_sources})
endif()

# write header for all template files
file(GLOB templates *.template.cc)
foreach(template ${templates})
    get_filename_component(template_WE ${template} NAME_WE)
    file(WRITE ${template_WE}.h "#include \"${template_WE}.template.h\"\n#include \"${template_WE}.template.cc\"")
endforeach()

# compile to a library named $dir
# message(STATUS "create ${dir}")

if (NOT "${sources}" STREQUAL "")
    add_library(${dir} STATIC ${sources})
endif()

