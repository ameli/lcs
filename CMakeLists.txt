# =====================================================================================
#
#       Filename:  CMakeLists.txt
#
#    Description:  Lagrangian Coherence Structures
#
#        Version:  1.0
#        Created:  02/12/2014 05:16:02 PM
#       Revision:  none
#       Compiler:  gcc
#
#         Author:  Siavash Ameli
#   Organization:  University of California, Berkeley
#
# =====================================================================================

cmake_minimum_required(VERSION 2.8)
project(LagrangianCoherenceStrutures CXX)

# ======================
# VTK/ParaView Inclusion
# ======================

# ParaView plugin option
set(BUILD_PARAVIEW_PLUGIN ON CACHE BOOL "Build ParaView plugin?")

# Avoid confliction of Paraview and VTK inclusion
if(BUILD_PARAVIEW_PLUGIN)
    # include ParaView instead of VTK
    find_package(ParaView REQUIRED)
    if(ParaView_FOUND)
        include(${PARAVIEW_USE_FILE})
    else(ParaView_FOUND)
        message(FATAL_ERROR "ParaView not found.")
    endif(ParaView_FOUND)
else(BUILD_PARAVIEW_PLUGIN)
    # include VTK
    find_package(VTK REQUIRED)
    if(VTK_FOUND)
        include(${VTK_USE_FILE})
    else(VTK_FOUND)
        message(FATAL_ERROR "VTK not found.")
    endif(VTK_FOUND)
endif(BUILD_PARAVIEW_PLUGIN)

# ==============
# Compiler Flags
# ==============

# Initialize CXX FLAGS (Not in WINDOWS)
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS                "-Wall -std=c++0x ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
    set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
endif()

# Compiler-specific activations
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.4 OR GCC_VERSION VERSION_EQUAL 4.4))
        message(FATAL_ERROR "${PROJECT_NAME} Requires g++ 4.4 or greater.")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
# else("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    # Windows
endif()

# ============
# OpenMP Flags
# ============

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FALGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

# ================
# Source Inclusion
# ================

# Headers
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${PROJECT_INCLUDE_DIR})

# Main Source
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(MAIN_EXECUTABLE_FILENAME "lcs")
add_executable(${MAIN_EXECUTABLE_FILENAME} 
    ${PROJECT_SOURCE_DIR}/main.cxx
    ${PROJECT_SOURCE_DIR}/Deformation.cxx
    ${PROJECT_SOURCE_DIR}/IdentifyStructures.cxx
    ${PROJECT_SOURCE_DIR}/SmoothStructuredPoints.cxx
    ${PROJECT_SOURCE_DIR}/ImageDataDerivatives.cxx
    ${PROJECT_SOURCE_DIR}/Matrix.cxx
    ${PROJECT_SOURCE_DIR}/Tensor.cxx
    ${PROJECT_SOURCE_DIR}/Vector.cxx
    ${PROJECT_SOURCE_DIR}/CubeCell.cxx
    ${PROJECT_SOURCE_DIR}/PointDataType.cxx)
target_link_libraries(${MAIN_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

# Test Source
set(BUILD_TEST ON CACHE BOOL "Build Test?")
if(BUILD_TEST)
    # Test Vector
    set(TEST_EXECUTABLE_FILENAME "TestVector")
    set(PROJECT_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_executable(${TEST_EXECUTABLE_FILENAME}
        ${PROJECT_SOURCE_DIR}/Vector.cxx
        ${PROJECT_TEST_DIR}/TestVector.cxx)
    target_link_libraries(${TEST_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

    # Test Matrix
    set(TEST_EXECUTABLE_FILENAME "TestMatrix")
    set(PROJECT_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_executable(${TEST_EXECUTABLE_FILENAME}
        ${PROJECT_SOURCE_DIR}/Matrix.cxx
        ${PROJECT_TEST_DIR}/TestMatrix.cxx)
    target_link_libraries(${TEST_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

    # Generate Sample Data
    set(TEST_EXECUTABLE_FILENAME "GenerateSampleData")
    set(PROJECT_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_executable(${TEST_EXECUTABLE_FILENAME}
        ${PROJECT_TEST_DIR}/GenerateSampleData.cxx)
    target_link_libraries(${TEST_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

    # Test Identify Structures
    set(TEST_EXECUTABLE_FILENAME "TestIdentifyStructures")
    set(PROJECT_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_executable(${TEST_EXECUTABLE_FILENAME}
        ${PROJECT_SOURCE_DIR}/IdentifyStructures.cxx
        ${PROJECT_SOURCE_DIR}/SmoothStructuredPoints.cxx
        ${PROJECT_SOURCE_DIR}/ImageDataDerivatives.cxx
        ${PROJECT_SOURCE_DIR}/CubeCell.cxx
        ${PROJECT_SOURCE_DIR}/Matrix.cxx
        ${PROJECT_SOURCE_DIR}/Vector.cxx
        ${PROJECT_SOURCE_DIR}/Tensor.cxx
        ${PROJECT_SOURCE_DIR}/PointDataType.cxx
        ${PROJECT_TEST_DIR}/TestIdentifyStructures.cxx)
    target_link_libraries(${TEST_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

    # Test Smooth Structured Points
    set(TEST_EXECUTABLE_FILENAME "TestSmoothStructuredPoints")
    set(PROJECT_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_executable(${TEST_EXECUTABLE_FILENAME}
        ${PROJECT_SOURCE_DIR}/SmoothStructuredPoints.cxx
        ${PROJECT_SOURCE_DIR}/Deformation.cxx
        ${PROJECT_TEST_DIR}/TestSmoothStructuredPoints.cxx)
    target_link_libraries(${TEST_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})
endif(BUILD_TEST)

# ==================
# Output Directories
# ==================

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Directory for all executables.")
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Directory for all libraries.")

# ===============
# ParaView Plugin
# ===============

if(BUILD_PARAVIEW_PLUGIN)
    set(PROJECT_PLUGIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/plugin)
    set(PROJECT_PLUGIN_DOCUMENTATION_DIR ${CMKAE_CURRENT_SOURCE_DIR}/doc)
    set(PluginVersion "1.0")

    # Deformation Filter
    ADD_PARAVIEW_PLUGIN(
        DeformationPlugin
        ${PluginVersion}
        DOCUMENTATION_DIR       ${CMAKE_CURRENT_SOURCE_DIR}/doc
        SERVER_MANAGER_XML      ${PROJECT_PLUGIN_DIR}/Deformation.xml
        GUI_RESOURCE_FILES      ${PROJECT_PLUGIN_DIR}/DeformationGUI.xml
        GUI_RESOURCES           ${PROJECT_PLUGIN_DIR}/DeformationICN.qrc
        SERVER_MANAGER_SOURCES  ${PROJECT_SOURCE_DIR}/Deformation.cxx)

    # Identify Structures Filter
    ADD_PARAVIEW_PLUGIN(
        IdentifyStructuresPlugin
        ${PluginVersion}
        DOCUMENTATION_DIR       ${CMAKE_CURRENT_SOURCE_DIR}/doc
        SERVER_MANAGER_XML      ${PROJECT_PLUGIN_DIR}/IdentifyStructures.xml
        GUI_RESOURCE_FILES      ${PROJECT_PLUGIN_DIR}/IdentifyStructuresGUI.xml
        GUI_RESOURCES           ${PROJECT_PLUGIN_DIR}/IdentifyStructuresICN.qrc
        SERVER_MANAGER_SOURCES  ${PROJECT_SOURCE_DIR}/IdentifyStructures.cxx
        SOURCES                 ${PROJECT_SOURCE_DIR}/SmoothStructuredPoints.cxx
                                ${PROJECT_SOURCE_DIR}/ImageDataDerivatives.cxx
                                ${PROJECT_SOURCE_DIR}/CubeCell.cxx
                                ${PROJECT_SOURCE_DIR}/PointDataType.cxx
                                ${PROJECT_SOURCE_DIR}/Vector.cxx
                                ${PROJECT_SOURCE_DIR}/Matrix.cxx
                                ${PROJECT_SOURCE_DIR}/Tensor.cxx)
endif(BUILD_PARAVIEW_PLUGIN)

# ============
# Installation
# ============

option(INSTALL_PROJECT "Install project on host system" ON)

if(${INSTALL_PROJECT})
    set(PROJECT_INSTALL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/install)
    install(TARGETS ${MAIN_EXECUTABLE_FILENAME} DESTINATION ${PROJECT_INSTALL_PATH})
endif(${INSTALL_PROJECT})

# ================
# Generate Package
# ================

option(GENERATE_PACKAGE "Generate system dependent package" ON)
if(${GENERATE_PACKAGE})
    # include(InstallRequiredSystemLibraries)  # not work in WIN

    # Custom variables
    set(AUTHOR "Siavash Ameli")
    set(PACKAGE_NAME "LCS")
    set(PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

    # cpack variables
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.md")
    set(CPACK_PACKAGE_NAME ${PACKAGE_NAME})
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER ${AUTHOR})
    set(CPACK_PACKAGE_CONTACT "sameli@berkeley.edu")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")

    # Package directory
    set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_DIR}/${PACKAGE_NAME}")

    # Package path
    # set(CPACK_TOPLEVEL_TAG ${PACKAGE_DIR})
    # set(CPACK_PACKAGING_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CPACK_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    
    # Package install path
    if(NOT DEFINED PROJECT_INSTALL_PATH)
        set(CPACK_SET_DESTDIR true)
        set(CPACK_INSTALL_PREFIX ${CMAKE_CURRENT_PROJECT}/bin)
    endif(NOT DEFINED PROJECT_INSTALL_PATH)

    # Binary option
    option(BINARY_PACKAGE "Generate binary package" ON)

    # OS
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        if(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "DEB")
        else(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "SOURCE_TGZ")
        endif(${BINARY_PACKAGE})
    elseif()
        # TODO: other OS
    endif()

    include(CPack)
endif(${GENERATE_PACKAGE})
