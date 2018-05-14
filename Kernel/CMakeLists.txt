
#Include directories where the source files and the libraries are stored
include_directories(
        ${Mercury_SOURCE_DIR}/Kernel
        ${Mercury_BINARY_DIR}/Kernel
)

#For the Fortran-based linear algebra
add_subdirectory(Math)

if(Mercury_USE_MPI)
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

#Compile all .cc in the directories below, and add them to the library DPMBase
file(GLOB Boundaries_SRC "Boundaries/*.cc")
file(GLOB_RECURSE CG_SRC "CG/*.cc")
file(GLOB FMM_SRC "FMM/*.cc")
file(GLOB Math_SRC "Math/*.cc")
file(GLOB Math_Quad_SRC "Math/Lebedev_Quadrature/*.cc")
file(GLOB Particles_SRC "Particles/*.cc")
file(GLOB Strings_SRC "Strings/*.cc")
file(GLOB VTKWriter_SRC "VTKWriter/*.cc")
file(GLOB Walls_SRC "Walls/*.cc")
file(GLOB Interactions_SRC "Interactions/*.cc")
file(GLOB InteractionsN_SRC "Interactions/NormalForceInteractions/*.cc")
file(GLOB InteractionsF_SRC "Interactions/FrictionForceInteractions/*.cc")
file(GLOB InteractionsA_SRC "Interactions/AdhesiveForceInteractions/*.cc")
file(GLOB Species_SRC "Species/*.cc")
file(GLOB SpeciesN_SRC "Species/NormalForceSpecies/*.cc")
file(GLOB SpeciesF_SRC "Species/FrictionForceSpecies/*.cc")
file(GLOB SpeciesA_SRC "Species/AdhesiveForceSpecies/*.cc")
file(GLOB Nurbs_SRC "Nurbs/*.cc")


#Basic library without hGrid
add_library(DPMBase STATIC
        ${Boundaries_SRC}
        ${CG_SRC}
        ${FMM_SRC}
        ${Math_SRC}
        ${Math_Quad_SRC}
        ${Particles_SRC}
        ${Strings_SRC}
        ${VTKWriter_SRC}
        ${Walls_SRC}

        BaseObject.cc
        BaseInteractable.cc
        Logger.cc
        File.cc
        DPMBase.cc
        MpiContainer.cc
        MpiDataClass.cc



        BoundaryHandler.cc
        PeriodicBoundaryHandler.cc
        InteractionHandler.cc
        ParticleHandler.cc
        SpeciesHandler.cc
        WallHandler.cc
        DomainHandler.cc
        Domain.cc

        ${Interactions_SRC}
        ${InteractionsN_SRC}
        ${InteractionsF_SRC}
        ${InteractionsA_SRC}

        ${Species_SRC}
        ${SpeciesN_SRC}
        ${SpeciesF_SRC}
        ${SpeciesA_SRC}

        ${Nurbs_SRC}

        ${Mercury_BINARY_DIR}/Kernel/CMakeDefinitions.cc
)

target_link_libraries(DPMBase LinearAlgebra)

if(Mercury_BACKTRACE_DEMANGLE)
    if (Mercury_BACKTRACE_LINK_DL)
        target_link_libraries(DPMBase dl)
    endif()
    if (Mercury_BACKTRACE_LINK_EXECINFO)
        target_link_libraries(DPMBase execinfo)
    endif()
endif()

if(Mercury_USE_MPI)
    set( DEP_COMPILE_FLAGS "${DEP_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
    set( DEP_LINKER_FLAGS "${DEP_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
    target_link_libraries(DPMBase "${MPI_LIBRARIES}")
endif()



#Could later create a 2D and 3D library seperatly, but sure how many we want to create
add_library(MercuryBase STATIC
        HGrid.cc
        MercuryBase.cc
        Mercury2D.cc
        Mercury3D.cc
        HGridOptimiser.cc)

target_link_libraries(MercuryBase DPMBase)

#Required for Chute drivers
add_library(Chute STATIC
        ChuteBottom.cc
        ChuteWithHopper.cc
        Chute.cc)

target_link_libraries(Chute MercuryBase)
