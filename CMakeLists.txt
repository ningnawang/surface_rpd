# ###############################################################################
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
set(CMAKE_CUDA_ARCHITECTURES 61 75)
project(CXX_TEST LANGUAGES CXX)

# ###############################################################################
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_STANDARD 17)

# ###############################################################################
find_package(OpenMP)
find_package(CGAL REQUIRED)
set(CMAKE_BUILD_TYPE "Release")

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# add dependencies
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(AUTO_DOWNLOAD "Auto download dependencies" ON)

include(rpdDependencies)

add_subdirectory(extern/mshloader)

# ###############################################################################
# For VolumeVoronoiGPU binary
# ###############################################################################
set(CXX_HEADER_LIST
    "src/matfp/geogram/predicates/powerside1.h"
    "src/matfp/geogram/predicates/powerside2.h"
    "src/matfp/geogram/predicates/powerside3.h"
    "src/matfp/geogram/predicates/powerside4.h"
    "src/matfp/geogram/predicates.h"
)

set(CXX_SOURCE_LIST
    "src/main_gui_cxx.cxx"
    "src/main.cxx"

    "src/io.cxx"
    "src/medial_sphere.cxx"
    "src/triangulation.cxx"

    "src/matfp/geogram/predicates.cpp"

    # "src/matfp/geogram/RPD.cpp"
    # "src/matfp/geogram/RPD.h"
    # "src/matfp/geogram/generic_RPD.h"
    # "src/matfp/geogram/generic_RPD_vertex.h"
    # "src/matfp/geogram/generic_RPD_polygon.cpp"
    # "src/matfp/geogram/generic_RPD_polygon.h"
    # "src/matfp/geogram/generic_RPD_cell.cpp"
    # "src/matfp/geogram/generic_RPD_cell.h"
    # "src/matfp/geogram/generic_RPD_utils.h"
    # "src/matfp/geogram/RPD_mesh_builder.cpp"
    # "src/matfp/geogram/RPD_mesh_builder.h"
    # "src/matfp/geogram/RPD_callback.cpp"
    # "src/matfp/geogram/RPD_callback.h"
)

add_executable(${PROJECT_NAME} ${CXX_HEADER_LIST} ${CXX_SOURCE_LIST})
target_compile_definitions(${PROJECT_NAME} PUBLIC -Dgeogram_EXPORTS)
target_include_directories(${PROJECT_NAME} PRIVATE src)

if(${MSVC})
    if(${CMAKE_VERSION} VERSION_GREATER 3.6.0)
        set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "${PROJECT_NAME}")
    endif()
endif()

target_link_directories(${PROJECT_NAME} PUBLIC "${CMAKE_BINARY_DIR}/lib")
target_link_libraries(${PROJECT_NAME}
    optimized mshloader
    optimized polyscope
    optimized geogram
    optimized CGAL::CGAL
)

# ###############################################################################
# For VORO_GUI binary
# ###############################################################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -pthread -D_DEBUG")

add_executable(VORO_GUI src/test_gui.cpp)
target_link_libraries(VORO_GUI PUBLIC polyscope mshloader)