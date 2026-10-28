cmake_minimum_required(VERSION 3.21)

project(zspace_toolsets VERSION 3.0.0 LANGUAGES CXX)

option(ZSPACE_TOOLSETS_BUILD_SHARED "Build zSpace toolsets as a DLL" ON)
option(ZSPACE_TOOLSETS_BUILD_TESTS "Build zSpace toolsets smoke tests" OFF)
option(ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE "Build against a sibling zspace_core checkout" ON)
set(ZSPACE_CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../zspace_core" CACHE PATH "Path to local zspace_core checkout")
set(ZSPACE_SDK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sdk/zspace" CACHE PATH "Path to installed zSpace SDK")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(ZSPACE_TOOLSETS_WITH_TETGEN "Enable harmonic slicing with an external TetGen library" OFF)
option(ZSPACE_TOOLSETS_SLICING_WITH_OPENUSD "Enable tagged slicer USD IO with an external OpenUSD SDK" OFF)
option(ZSPACE_TOOLSETS_SLICING_OPENMP "Enable OpenMP in the slicing engine" OFF)
set(ZSPACE_TETGEN_ROOT "" CACHE PATH "TetGen installation root (include/tetgen.h and lib/tetgen)")

set(ZSPACE_TOOLSETS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(ZSPACE_TOOLSETS_BUILD_SHARED)
  set(ZSPACE_TOOLSETS_LIBRARY_TYPE SHARED)
else()
  set(ZSPACE_TOOLSETS_LIBRARY_TYPE STATIC)
endif()

if(TARGET zspace::core AND TARGET zspace::interface AND TARGET zspace::io)
  message(STATUS "zspace_toolsets: using existing zSpace targets from parent project")
elseif(ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE)
  if(NOT EXISTS "${ZSPACE_CORE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE=ON requires ZSPACE_CORE_DIR to point to zspace_core")
  endif()

  set(ZSPACE_BUILD_INTERFACE ON CACHE BOOL "Build zSpace interface DLL" FORCE)
  set(ZSPACE_BUILD_IO ON CACHE BOOL "Build zSpace IO DLL" FORCE)
  set(ZSPACE_BUILD_DISPLAY OFF CACHE BOOL "Build zSpace display DLL" FORCE)
  set(ZSPACE_BUILD_INTEROP OFF CACHE BOOL "Build zSpace interop DLL" FORCE)
  set(ZSPACE_BUILD_TESTS OFF CACHE BOOL "Build zSpace tests" FORCE)
  set(ZSPACE_BUILD_SHARED ON CACHE BOOL "Build zSpace DLLs" FORCE)
  set(ZSPACE_WITH_OPENGL OFF CACHE BOOL "Disable zSpace display OpenGL" FORCE)
  set(ZSPACE_IO_WITH_TINYUSDZ OFF CACHE BOOL "Disable TinyUSDZ download for toolsets scaffold" FORCE)

  add_subdirectory("${ZSPACE_CORE_DIR}" "${CMAKE_BINARY_DIR}/zspace_core_build")
else()
  list(PREPEND CMAKE_PREFIX_PATH "${ZSPACE_SDK_DIR}")
  find_package(zspace 3 CONFIG REQUIRED)
endif()

option(ZSPACE_TOOLSETS_SLICING_ONLY "Build only the slicing implementation" OFF)
include("${ZSPACE_TOOLSETS_ROOT}/cmake/SlicingSources.cmake")
zspace_slicing_sources(ZSPACE_TOOLSETS_SOURCES "${ZSPACE_TOOLSETS_ROOT}" "${ZSPACE_TOOLSETS_WITH_TETGEN}")
if(NOT ZSPACE_TOOLSETS_SLICING_ONLY)
  file(GLOB_RECURSE other_sources CONFIGURE_DEPENDS
    "${ZSPACE_TOOLSETS_ROOT}/src/zToolsets/*.cpp")
  list(FILTER other_sources EXCLUDE REGEX "/fabrication/slicing/")
  list(APPEND ZSPACE_TOOLSETS_SOURCES ${other_sources})
endif()
file(GLOB_RECURSE ZSPACE_TOOLSETS_HEADERS CONFIGURE_DEPENDS
  "${ZSPACE_TOOLSETS_ROOT}/include/zspace/*.h")

add_library(zSpace_Toolsets ${ZSPACE_TOOLSETS_LIBRARY_TYPE}
  ${ZSPACE_TOOLSETS_SOURCES}
  ${ZSPACE_TOOLSETS_HEADERS}
)
add_library(zspace::toolsets ALIAS zSpace_Toolsets)
set_target_properties(zSpace_Toolsets PROPERTIES EXPORT_NAME toolsets)

target_include_directories(zSpace_Toolsets
  PUBLIC
    $<BUILD_INTERFACE:${ZSPACE_TOOLSETS_ROOT}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(zSpace_Toolsets PUBLIC zspace::core zspace::interface zspace::io)

target_compile_definitions(zSpace_Toolsets
  PUBLIC
    $<$<BOOL:${ZSPACE_TOOLSETS_BUILD_SHARED}>:ZSPACE_TOOLSETS_DYNAMIC_LIBRARY>
    $<$<NOT:$<BOOL:${ZSPACE_TOOLSETS_BUILD_SHARED}>>:ZSPACE_TOOLSETS_STATIC_LIBRARY>
  PRIVATE
    $<$<BOOL:${ZSPACE_TOOLSETS_BUILD_SHARED}>:ZSPACE_TOOLSETS_EXPORTS>
)

if(ZSPACE_TOOLSETS_WITH_TETGEN)
  find_path(ZSPACE_TETGEN_INCLUDE_DIR tetgen.h HINTS "${ZSPACE_TETGEN_ROOT}" PATH_SUFFIXES include REQUIRED)
  find_library(ZSPACE_TETGEN_LIBRARY NAMES tetgen tet HINTS "${ZSPACE_TETGEN_ROOT}" PATH_SUFFIXES lib REQUIRED)
  target_include_directories(zSpace_Toolsets SYSTEM PRIVATE "${ZSPACE_TETGEN_INCLUDE_DIR}")
  target_compile_definitions(zSpace_Toolsets PUBLIC ZSPACE_TOOLSETS_WITH_TETGEN PRIVATE TETLIBRARY)
  target_link_libraries(zSpace_Toolsets PRIVATE "${ZSPACE_TETGEN_LIBRARY}")
endif()

if(ZSPACE_TOOLSETS_SLICING_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  # Layer-level parallelism owns scheduling; avoid nested Eigen threads.
  target_compile_definitions(zSpace_Toolsets PRIVATE EIGEN_DONT_PARALLELIZE)
  target_link_libraries(zSpace_Toolsets PRIVATE OpenMP::OpenMP_CXX)
endif()

if(ZSPACE_TOOLSETS_SLICING_WITH_OPENUSD)
  if(NOT COMMAND zspace_link_openusd)
    if(NOT EXISTS "${ZSPACE_CORE_DIR}/cmake/ZSpaceOpenUSD.cmake")
      message(FATAL_ERROR "Slicer USD IO requires the zspace_core OpenUSD CMake helper and ZSPACE_OPENUSD_ROOT")
    endif()
    include("${ZSPACE_CORE_DIR}/cmake/ZSpaceOpenUSD.cmake")
  endif()
  zspace_link_openusd(zSpace_Toolsets)
  target_compile_definitions(zSpace_Toolsets PRIVATE ZSPACE_TS3DP_OPENUSD)
endif()

if(MSVC)
  target_compile_options(zSpace_Toolsets PRIVATE /MP /W3 /permissive- /wd4251)
  target_compile_definitions(zSpace_Toolsets PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

install(TARGETS zSpace_Toolsets
  EXPORT zspaceToolsetsTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION bin
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include
)
install(DIRECTORY "${ZSPACE_TOOLSETS_ROOT}/include/zspace" DESTINATION include)
install(EXPORT zspaceToolsetsTargets FILE zspaceToolsetsTargets.cmake NAMESPACE zspace:: DESTINATION lib/cmake/zspace_toolsets)

if(ZSPACE_TOOLSETS_BUILD_TESTS)
  include(CTest)
  add_executable(zspace_toolsets_face_tags tests/smoke/zspace_toolsets_face_tags.cpp)
  target_link_libraries(zspace_toolsets_face_tags PRIVATE zspace::toolsets)
  add_test(NAME zspace_toolsets.slicing.face_tags COMMAND zspace_toolsets_face_tags 0)
  if(ZSPACE_TOOLSETS_WITH_TETGEN)
    add_test(NAME zspace_toolsets.slicing.face_tags_harmonic COMMAND zspace_toolsets_face_tags 1)
  endif()
  add_executable(zspace_toolsets_printsyn tests/smoke/zspace_toolsets_printsyn.cpp)
  target_link_libraries(zspace_toolsets_printsyn PRIVATE zspace::toolsets)
  add_test(NAME zspace_toolsets.slicing.postprocessing COMMAND zspace_toolsets_printsyn)
  add_executable(zspace_toolsets_print_mesh tests/smoke/zspace_toolsets_print_mesh.cpp)
  target_link_libraries(zspace_toolsets_print_mesh PRIVATE zspace::toolsets)
  add_test(NAME zspace_toolsets.slicing.print_mesh COMMAND zspace_toolsets_print_mesh)
  add_executable(zspace_toolsets_slicing tests/smoke/zspace_toolsets_slicing.cpp)
  target_link_libraries(zspace_toolsets_slicing PRIVATE zspace::toolsets)
  add_test(NAME zspace_toolsets.slicing.topology COMMAND zspace_toolsets_slicing topology)
  add_test(NAME zspace_toolsets.slicing.invalid COMMAND zspace_toolsets_slicing invalid)
  add_test(NAME zspace_toolsets.slicing.features COMMAND zspace_toolsets_slicing features)
  add_test(NAME zspace_toolsets.slicing.parallel COMMAND zspace_toolsets_slicing parallel)
  add_executable(zspace_toolsets_sampling3d tests/smoke/zspace_toolsets_sampling3d.cpp)
  target_link_libraries(zspace_toolsets_sampling3d PRIVATE zspace::toolsets)
  add_test(NAME zspace_toolsets.slicing.sampling3d COMMAND zspace_toolsets_sampling3d)
  set_tests_properties(zspace_toolsets.slicing.parallel PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
  if(ZSPACE_TOOLSETS_WITH_TETGEN)
    add_executable(zspace_toolsets_iso_relaxation tests/smoke/zspace_toolsets_iso_relaxation.cpp)
    target_link_libraries(zspace_toolsets_iso_relaxation PRIVATE zspace::toolsets)
    add_test(NAME zspace_toolsets.slicing.iso_relaxation COMMAND zspace_toolsets_iso_relaxation)
    add_test(NAME zspace_toolsets.slicing.harmonic COMMAND zspace_toolsets_slicing harmonic)
    add_test(NAME zspace_toolsets.slicing.tet_refinement COMMAND zspace_toolsets_slicing tet_refinement)
  endif()
  if(NOT ZSPACE_TOOLSETS_SLICING_ONLY)
    add_executable(zspace_toolsets_force_polyhedra tests/smoke/zspace_toolsets_force_polyhedra.cpp)
    target_link_libraries(zspace_toolsets_force_polyhedra PRIVATE zspace::toolsets)
    add_test(NAME zspace_toolsets.force_polyhedra COMMAND zspace_toolsets_force_polyhedra)
    add_executable(zspace_toolsets_stairs tests/smoke/zspace_toolsets_stairs.cpp)
    target_link_libraries(zspace_toolsets_stairs PRIVATE zspace::toolsets)
    add_test(NAME zspace_toolsets.stairs COMMAND zspace_toolsets_stairs)
    add_executable(zspace_toolsets_panelisation tests/smoke/zspace_toolsets_panelisation.cpp)
    target_link_libraries(zspace_toolsets_panelisation PRIVATE zspace::toolsets)
    add_test(NAME zspace_toolsets.panelisation COMMAND zspace_toolsets_panelisation)
    add_executable(zspace_toolsets_ctubes tests/smoke/zspace_toolsets_ctubes.cpp)
    target_link_libraries(zspace_toolsets_ctubes PRIVATE zspace::toolsets)
    add_test(NAME zspace_toolsets.ctubes COMMAND zspace_toolsets_ctubes)
    add_executable(zspace_toolsets_vault tests/smoke/zspace_toolsets_vault.cpp)
    target_link_libraries(zspace_toolsets_vault PRIVATE zspace::toolsets)
    add_test(NAME zspace_toolsets.vault COMMAND zspace_toolsets_vault)
  endif()
endif()

message(STATUS "zspace_toolsets Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  zSpace source mode: ${ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE}")
message(STATUS "  zSpace core dir: ${ZSPACE_CORE_DIR}")

