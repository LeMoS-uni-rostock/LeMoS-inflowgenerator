set(PRJ superposeIsotropicTurbulence)

set(SRC superposeIsotropicTurbulence.C)

set(OF_INCLUDE_DIRS
#   OpenFOAM finiteVolume OSspecific/POSIX lduSolvers 
  turbulenceModels/incompressible/LES 
#   turbulenceModels/LES/LESfilters
#   meshTools sampling surfMesh
)

set(OF_LIBS 
# finiteVolume LESfilters sampling
)

set(INCLUDE_DIRS
  ${inflowGeneratorBCOF23x_SOURCE_DIR}
  ${inflowGeneratorBCOF23x_SOURCE_DIR}/turbulentStructures/turbulentStructure
  ${fieldDataProviderOF23x_SOURCE_DIR}
)

set(LIBS
 testcases
 refdata
 uniof
)
set(IS_OF_TARGS 
 inflowGeneratorBC
 globalPatch
)

setup_exe_target_OF(${PRJ} "${SRC}" "${OF_INCLUDE_DIRS}" "${OF_LIBS}" "${INCLUDE_DIRS}" "${LIBS}" "${IS_OF_TARGS}")
linkToolkit_OF_Exe(${PRJ})
