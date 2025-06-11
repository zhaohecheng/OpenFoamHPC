add_definitions(-DWM_LABEL_SIZE=64)
add_definitions(-DWM_DP)

add_executable(simpleFoamHPC
        solvers/incompressible/simpleFoamHPC/simpleFoamHPC.C
        ${CMAKE_SOURCE_DIR}/src/dynamicFvMesh/dynamicFvMesh/dynamicFvMesh.C
        ${CMAKE_SOURCE_DIR}/src/finiteVolume/finiteVolume/fvm/fvmDiv.C
        ${CMAKE_SOURCE_DIR}/src/finiteVolume/finiteVolume/convectionSchemes/convectionScheme/convectionScheme.C
        ${CMAKE_SOURCE_DIR}/src/finiteVolume/cfdTools/general/MRF/MRFZoneList.C
        ${CMAKE_SOURCE_DIR}/src/OpenFOAM/matrices/lduMatrix/lduMatrix/lduMatrix.C)

target_include_directories(simpleFoamHPC PRIVATE
        ${CMAKE_SOURCE_DIR}/src/finiteVolume/lnInclude
        ${CMAKE_SOURCE_DIR}/src/meshTools/lnInclude
        ${CMAKE_SOURCE_DIR}/src/sampling/lnInclude
        ${CMAKE_SOURCE_DIR}/src/TurbulenceModels/turbulenceModels/lnInclude
        ${CMAKE_SOURCE_DIR}/src/TurbulenceModels/incompressible/lnInclude
        ${CMAKE_SOURCE_DIR}/src/transportModels
        ${CMAKE_SOURCE_DIR}/src/transportModels/incompressible/singlePhaseTransportModel
        ${CMAKE_SOURCE_DIR}/src/dynamicMesh/lnInclude
        ${CMAKE_SOURCE_DIR}/src/dynamicFvMesh/lnInclude
        ${CMAKE_SOURCE_DIR}/src/OpenFOAM/lnInclude
        ${CMAKE_SOURCE_DIR}/src/TurbulenceModels/incompressible/lnInclude/)