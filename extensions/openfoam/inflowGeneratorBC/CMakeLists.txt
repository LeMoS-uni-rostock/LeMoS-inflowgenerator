set(PRJ inflowGeneratorBC)

set(SRC 
 boostRandomGen.C 
# inflowGeneratorBaseFvPatchVectorField.C
 inflowGeneratorFvPatchVectorFieldFwd.C
# modalTurbulenceFvPatchVectorField.C
# turbulentSpotVolumeInserters.C
# turbulentForcingFvOption.C
# turbulentScalarFvPatchField.C
 vortonqueue.cpp
 inflowgeneratorbase.cpp
 structurebasedinflowgenerator.cpp
 turbulentStructures/turbulentstructure.cpp
 turbulentStructures/anisotropicvorton.cpp
)

set(OF_INCLUDE_DIRS
)

set(OF_LIBS 
)

set(INCLUDE_DIRS
)

set(LIBS 
uniof
)

set(IS_OF_TARGS
 globalPatch
 fieldDataProvider
)

#set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}") #hack
#add_subdirectory(turbulentStructures/turbulentStructure)
#add_subdirectory(turbulentStructures/hatSpot)
#add_subdirectory(turbulentStructures/gaussianSpot)
#add_subdirectory(turbulentStructures/decayingTurbulenceVorton)
#add_subdirectory(turbulentStructures/anisotropicVorton)
# add_subdirectory(turbulentStructures/anisotropicVorton2)
#add_subdirectory(turbulentStructures/combinedVorton)
## add_subdirectory(turbulentStructures/homogeneousTurbulenceVorton)

setup_lib_target_OF(${PRJ} "${SRC}" "${OF_INCLUDE_DIRS}" "${OF_LIBS}" "${INCLUDE_DIRS}" "${LIBS}" "${IS_OF_TARGS}")
linkToolkit_OF_ShLib( ${PRJ} )
