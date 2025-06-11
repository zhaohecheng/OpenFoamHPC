add_definitions(-DWM_LABEL_SIZE=64)
add_definitions(-DWM_DP)

add_executable(simpleFoamHPC solvers/incompressible/simpleFoamHPC/simpleFoamHPC.C)

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