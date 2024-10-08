# Declare minimum version of CMake required
cmake_minimum_required(VERSION 3.15)

#----------------OPTIONS------------------------------------------------------------------

message(" # Setting CMake Tool Mesh Partitioner")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set project name, version and languages
project(tool_meshConversorPar LANGUAGES C CXX Fortran)

# Configure MPI module
include(mpi)

# Configure HDF5 module
include(hdf5)

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#----------------COMPILE------------------------------------------------------------------
# Compile
if(SHARED_LIBS)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
    file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
    include_directories(${CMAKE_BINARY_DIR}/modules)
else()
    file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_constants.f90 
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_numerical_params.f90 
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi.f90 
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_utils.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_ijk_indices.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_custom_types.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_maths.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/quadrature_rules.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/elem_qua.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/elem_hex.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi_mesh.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_comms.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_comms_boundaries.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_hdf5.f90
                        ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_nvtx.f90
                        ${CMAKE_CURRENT_LIST_DIR}/src/mod_meshConversorTool.f90
                        ${CMAKE_CURRENT_LIST_DIR}/src/tool_meshConversorPar.f90)
endif()

add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE Fortran)

add_subdirectory(lib_gempaInterface)

set_mpi()
set_gempa()
if(SHARED_LIBS)
    target_link_libraries(${PROJECT_NAME} lib_gempaInterface lib_sod2d)
else()
    set_hdf5()
    target_link_libraries(${PROJECT_NAME} lib_gempaInterface)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
if(SHARED_LIBS)
    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include)
endif()

# FFLAGS depend on the compiler
# TODO: check if can safely be removed later
# get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)