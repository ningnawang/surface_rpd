# ###############################################################################
# CMake download helpers
# ###############################################################################

# download external dependencies
include(rpdDownloadExternal)

# ###############################################################################
# Required dependencies
# ###############################################################################

# geogram
if(NOT TARGET geogram)
    tetwild_download_geogram()
    include(geogram) # check geogram.cmake
endif()

# # Boost
# if(TETWILD_WITH_HUNTER)
# hunter_add_package(Boost COMPONENTS thread system)
# endif()

# # fmt
# if(NOT TARGET fmt::fmt)
# rpd_download_fmt()
# add_subdirectory(${EXTERNAL_DIR}/fmt)
# endif()

# # spdlog
# if(NOT TARGET spdlog::spdlog)
# rpd_download_spdlog()
# add_library(spdlog INTERFACE)
# add_library(spdlog::spdlog ALIAS spdlog)
# target_include_directories(spdlog INTERFACE ${EXTERNAL_DIR}/spdlog/include)
# target_compile_definitions(spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)
# target_link_libraries(spdlog INTERFACE fmt::fmt)
# endif()

# # libigl
# if(NOT TARGET igl::core)
# rpd_download_libigl()
# endif()

# # pymesh loaders
# add_subdirectory(${EXTERNAL_DIR}/pymesh)

# # CL11
# if(NOT TARGET CLI11::CLI11)
# tetwild_download_cli11()
# add_subdirectory(${EXTERNAL_DIR}/cli11)
# endif()

# polyscope
if(NOT TARGET polyscope)
    rpd_download_polyscope()
    add_subdirectory(${EXTERNAL_DIR}/polyscope)
endif()

# json
if(NOT TARGET nlohmann_json::nlohmann_json)
    rpd_download_json()
endif()