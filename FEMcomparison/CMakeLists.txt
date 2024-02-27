cmake_minimum_required (VERSION 3.13)

project (FEMcomparison)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(cmake/StandardProjectSettings.cmake)

## Sanitizer options if supported by compiler
#include(cmake/Sanitizers.cmake)
#add_library(project_options INTERFACE)
#enable_sanitizers(project_options)

set(TOOLS_LIB_TYPE STATIC)
add_library(Tools ${TOOLS_LIB_TYPE} "")

#enabling MKL library
option(FEMCOMPARISON_USING_MKL "Whether the MKL library will be linked in" OFF)
if(FEMCOMPARISON_USING_MKL)
    include(cmake/EnableMKL.cmake)
    enable_mkl(Tools)
endif(FEMCOMPARISON_USING_MKL)

option(OPTMIZE_RUN_TIME "Whether the simulation is accelarated." OFF)
if(OPTMIZE_RUN_TIME)
    add_definitions(-DOPTMIZE_RUN_TIME)
endif(OPTMIZE_RUN_TIME)

option(TIMER_CONTRIBUTE "Whether timer measures is required." OFF)
if(TIMER_CONTRIBUTE)
    add_definitions(-DTIMER_CONTRIBUTE)
endif(TIMER_CONTRIBUTE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DFEMCOMPARISON_DEBUG)
endif()

option(FEMCOMPARISON_TIMER "Timer measurement" OFF)
if(FEMCOMPARISON_TIMER)
    add_definitions(-DFEMCOMPARISON_TIMER)
endif()

# Allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

# Enables NeoPZ library
include(cmake/EnableNeoPZ.cmake)
enable_pz()

# Enables log4cxx library
#TODO: Change all LOG4CXX calls to PZ_LOG or correctly implement the LOG4CXX call for this project;
#TODO: The current call for definition LOG4CXX is a dummy call that always return false.
#include(cmake/EnableLog4CXX.cmake)
#enable_log4cxx()

find_package(Threads)

find_package(TBB)

# This option enables a lot of warnings and treat them as errors, to ensure
# good programming practices are used. Since its behaviour is extreme, it
# should be turned off by default.
option(ENABLE_ALL_WARNINGS "Enable a lot of compiler warnings and treat them as errors" OFF)
if (ENABLE_ALL_WARNINGS)
    add_library(project_warnings INTERFACE FEM/MeshInit.cpp FEM/MeshInit.h FEM/Solver.cpp FEM/Solver.h FEM/Output.cpp FEM/Output.h FEM/InputTreatment.cpp FEM/InputTreatment.h)
    include(cmake/CompilerWarnings.cmake)
    set_project_warnings(project_warnings)
    set(CMAKE_CXX_FLAGS "")
endif ()

# Add subdirectories
add_subdirectory(CompEl)
add_subdirectory(Tools)
add_subdirectory(Material)
add_subdirectory(FEM)
add_subdirectory(CreateApproxSpace)

