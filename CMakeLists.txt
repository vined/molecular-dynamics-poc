cmake_minimum_required(VERSION 3.5)
project(molecular_dynamics_poc)

set(CMAKE_CXX_STANDARD 11)

add_executable(molecular_dynamics_poc src/main.cc
        src/Vectors.cc src/Vectors.h
        src/Parameters.cc src/Parameters.h
        src/Atoms.cc src/Atoms.h
        src/Energies.cc src/Energies.h
        src/utils/OutputUtils.cc src/utils/OutputUtils.h
        src/StepResult.h
        src/Quaternions.cc src/Quaternions.h
        src/LeapFrogResult.h
        src/utils/FileUtils.cc src/utils/FileUtils.h
        src/Sites.cc src/Sites.h
        src/Matrix.cc src/Matrix.h
        src/Molecules.cc src/Molecules.h
        src/PredictorCorrector.cc src/PredictorCorrector.h)
