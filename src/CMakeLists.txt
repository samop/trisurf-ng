# --- Compute TS_VERSION like Automake: git describe ... ---
set(GITVERSION "unknown")
find_package(Git QUIET)
if (Git_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GITVERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()

# Common compile flags from AM_CFLAGS / AM_CPPFLAGS
set(TRISURF_COMMON_CFLAGS
  -Wall
  -fgnu89-inline
  -march=native
  -fopenmp
  -fopt-info-vec-optimized
  -fopt-info-vec-missed
)

# Define TS_VERSION="..."
# (CMake will handle escaping correctly here)
add_compile_definitions(TS_VERSION="${GITVERSION}")

# --- Library: libtrisurf.la -> libtrisurf ---
add_library(trisurf SHARED
  general.c
  vertex.c
  bond.c
  triangle.c
  cell.c
  vesicle.c
  initial_distribution.c
  io.c
  dumpstate.c
  frame.c
  energy.c
  timestep.c
  vertexmove.c
  bondflip.c
  poly.c
  stats.c
  sh.c
  shcomplex.c
  constvol.c
  snapshot.c
  restore.c
  cluster.c
)

# Keep output name as "trisurf" (libtrisurf.so)
set_target_properties(trisurf PROPERTIES OUTPUT_NAME "trisurf")

target_compile_options(trisurf PRIVATE ${TRISURF_COMMON_CFLAGS})
target_include_directories(trisurf PRIVATE ${TRISURF_COMMON_INCLUDES})

target_link_libraries(trisurf PRIVATE ${TRISURF_COMMON_LIBS})



# -- static library ---

if(ENABLE_STATIC)
# Static library
add_library(trisurf_static STATIC
  general.c
  vertex.c
  bond.c
  triangle.c
  cell.c
  vesicle.c
  initial_distribution.c
  io.c
  dumpstate.c
  frame.c
  energy.c
  timestep.c
  vertexmove.c
  bondflip.c
  poly.c
  stats.c
  sh.c
  shcomplex.c
  constvol.c
  snapshot.c
  restore.c
  cluster.c
)

set_target_properties(trisurf_static PROPERTIES
  OUTPUT_NAME trisurf
)

target_compile_options(trisurf_static PRIVATE ${TRISURF_COMMON_CFLAGS})
target_include_directories(trisurf_static PRIVATE ${LIBXML2_INCLUDES})
target_link_libraries(trisurf_static PRIVATE ${TRISURF_COMMON_LIBS})
endif()




# --- Executables ---
add_executable(trisurf_bin main.c)
set_target_properties(trisurf_bin PROPERTIES OUTPUT_NAME "trisurf")
target_compile_options(trisurf_bin PRIVATE ${TRISURF_COMMON_CFLAGS})
target_include_directories(trisurf_bin PRIVATE ${TRISURF_COMMON_INCLUDES})
target_link_libraries(trisurf_bin PRIVATE trisurf ${TRISURF_COMMON_LIBS})

add_executable(tsmeasure tsmeasure.c)
target_compile_options(tsmeasure PRIVATE ${TRISURF_COMMON_CFLAGS})
target_include_directories(tsmeasure PRIVATE ${TRISURF_COMMON_INCLUDES})
target_link_libraries(tsmeasure PRIVATE trisurf ${TRISURF_COMMON_LIBS})

add_executable(tspoststat tspoststat.c)
target_compile_options(tspoststat PRIVATE ${TRISURF_COMMON_CFLAGS})
target_include_directories(tspoststat PRIVATE ${TRISURF_COMMON_INCLUDES})
target_link_libraries(tspoststat PRIVATE trisurf ${TRISURF_COMMON_LIBS})



# --- Install (Automake bin_PROGRAMS + lib_LTLIBRARIES) ---
install(TARGETS trisurf_bin tsmeasure tspoststat
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS trisurf
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# --- Install headers (pkginclude_HEADERS -> includedir/${PACKAGE}) ---
install(FILES
  general.h
  vertex.h
  bond.h
  triangle.h
  cell.h
  vesicle.h
  initial_distribution.h
  io.h
  dumpstate.h
  frame.h
  energy.h
  timestep.h
  vertexmove.h
  bondflip.h
  poly.h
  stats.h
  sh.h
  shcomplex.h
  constvol.h
  snapshot.h
  restore.h
  cluster.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

