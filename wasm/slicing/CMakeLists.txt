cmake_minimum_required(VERSION 3.21)
project(zspace_toolset_wasm LANGUAGES CXX)
if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "Configure with Emscripten/emcmake")
endif()
get_filename_component(TOOLSETS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(ZSPACE_CORE_DIR "${TOOLSETS_ROOT}/../zspace_core" CACHE PATH "Core source checkout")
set(ZSPACE_TETGEN_SOURCE_DIR "${TOOLSETS_ROOT}/wasm/cache/tetgen" CACHE PATH "TetGen source checkout")
# Core owns its implementation, compatibility shims, bindings and USD settings.
set(ZSPACE_WASM_OBJECT_LIBRARY ON)
set(ZSPACE_WASM_WITH_FIELDS ON)
set(ZSPACE_WASM_WITH_USD ON)
set(ZSPACE_WASM_THREADS ON)
set(ZSPACE_WASM_OPENMP ON)
add_subdirectory("${ZSPACE_CORE_DIR}/wasm" "${CMAKE_CURRENT_BINARY_DIR}/core")
set(sources)
set(harmonic OFF)
if(EXISTS "${ZSPACE_TETGEN_SOURCE_DIR}/tetgen.cxx" AND EXISTS "${ZSPACE_TETGEN_SOURCE_DIR}/predicates.cxx")
  set(harmonic ON)
  list(APPEND sources "${ZSPACE_TETGEN_SOURCE_DIR}/tetgen.cxx" "${ZSPACE_TETGEN_SOURCE_DIR}/predicates.cxx")
endif()
include("${TOOLSETS_ROOT}/cmake/SlicingSources.cmake")
zspace_slicing_sources(slicing_sources "${TOOLSETS_ROOT}" "${harmonic}")
list(APPEND sources ${slicing_sources})
list(APPEND sources "${TOOLSETS_ROOT}/src/zToolsets/fabrication/slicing/utils/IO.cpp")
add_executable(zspace_toolset ${sources} "${CMAKE_CURRENT_SOURCE_DIR}/zspace_slicing_bridge.cpp")
zspace_wasm_configure_compilation(zspace_toolset)
# The compilation helper already adds the executable's OpenUSD link closure.
# Consume Core's objects directly so its private whole-archive dependencies
# are not linked a second time through the object library interface.
target_sources(zspace_toolset PRIVATE $<TARGET_OBJECTS:zspace_core_wasm_objects>)
set_target_properties(zspace_toolset PROPERTIES
  SUFFIX ".js" RUNTIME_OUTPUT_DIRECTORY "${TOOLSETS_ROOT}/wasm/out")
target_compile_definitions(zspace_toolset PRIVATE ZSPACE_TOOLSETS_WITH_USD_IO)

target_include_directories(zspace_toolset PRIVATE
  "${TOOLSETS_ROOT}/include" "${TOOLSETS_ROOT}/wasm/bridge")
target_compile_definitions(zspace_toolset PRIVATE ZSPACE_TOOLSETS_STATIC_LIBRARY)
if(harmonic)
  target_include_directories(zspace_toolset PRIVATE "${ZSPACE_TETGEN_SOURCE_DIR}")
  target_compile_definitions(zspace_toolset PRIVATE ZSPACE_TOOLSETS_WITH_TETGEN TETLIBRARY)
endif()

option(ZSPACE_WASM_USE_PCH "Use precompiled Eigen headers" ON)
if(ZSPACE_WASM_USE_PCH)
  target_precompile_headers(zspace_toolset PRIVATE <depends/Eigen/Dense>)
endif()
option(ZSPACE_WASM_FAST_LINK "Skip binaryen optimization for development" OFF)
set(zspace_wasm_link_opt -O2)
if(ZSPACE_WASM_FAST_LINK)
  set(zspace_wasm_link_opt -O0)
endif()
target_link_options(zspace_toolset PRIVATE ${zspace_wasm_link_opt} -fexceptions -pthread -fopenmp "-sMODULARIZE=1" "-sEXPORT_ES6=1"
  "-sENVIRONMENT=web,worker,node" "-sFORCE_FILESYSTEM=1" "--use-port=zlib" "-sALLOW_MEMORY_GROWTH=1" "-sSTACK_SIZE=2097152"
  "-sPTHREAD_POOL_SIZE=(Module['slicingThreadCount']=Math.max(1,Math.floor(Module['slicingThreadCount']||globalThis.navigator?.hardwareConcurrency||1)))-1" "-sPTHREAD_POOL_SIZE_STRICT=2"
  "-sDEFAULT_PTHREAD_STACK_SIZE=2097152"
  "-sEXPORTED_RUNTIME_METHODS=['HEAPF32','HEAPU32','UTF8ToString','ccall','PThread','FS']")
message(STATUS "Slicing WASM harmonic support: ${harmonic}")

add_custom_command(TARGET zspace_toolset POST_BUILD
  COMMAND "${CMAKE_COMMAND}" -E true
  BYPRODUCTS "${TOOLSETS_ROOT}/wasm/out/zspace_toolset.wasm")
