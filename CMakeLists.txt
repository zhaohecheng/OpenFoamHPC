cmake_minimum_required(VERSION 4.0)

add_definitions(-DWM_LABEL_SIZE=64)
add_definitions(-DWM_DP)

add_executable(
        buoyantSimpleFoam
        applications/solvers/heatTransfer/buoyantSimpleFoam/buoyantSimpleFoam.C
        src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C
        src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C)

target_include_directories(
        buoyantSimpleFoam PUBLIC
        src/finiteVolume/lnInclude
        src/meshTools/lnInclude
        src/sampling/lnInclude
        src/transportModels/compressible/lnInclude
        src/thermophysicalModels/basic/lnInclude
        src/thermophysicalModels/radiation/lnInclude
        src/TurbulenceModels/turbulenceModels/lnInclude
        src/TurbulenceModels/compressible/lnInclude
        src/regionFaModels/lnInclude
        src/OpenFOAM/lnInclude
)

add_executable(
        simpleFoam
        applications/solvers/incompressible/simpleFoam/simpleFoam.C
        src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C
        src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C)

target_include_directories(
        simpleFoam PUBLIC
        src/finiteVolume/lnInclude
        src/meshTools/lnInclude
        src/sampling/lnInclude
        src/TurbulenceModels/turbulenceModels/lnInclude
        src/TurbulenceModels/incompressible/lnInclude
        src/transportModels
        src/transportModels/incompressible/singlePhaseTransportModel
        src/dynamicMesh/lnInclude
        src/dynamicFvMesh/lnInclude
        src/OpenFOAM/lnInclude
        src/TurbulenceModels/incompressible/lnInclude/
)