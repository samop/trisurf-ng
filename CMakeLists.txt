cmake_minimum_required(VERSION 3.16)

project(trisurf_ng
  VERSION 1.0
  LANGUAGES C
)

include(GNUInstallDirs)

# --- Option equivalent to ./configure --enable-debug ---
option(ENABLE_DEBUG "Enable debugging support" OFF)

if (ENABLE_DEBUG)
  message(STATUS "Debugging enabled")
  add_compile_options(-g -O2 -DDEBUG)
else()
  message(STATUS "Debugging disabled")
  add_compile_options(-O3 -ffast-math)
endif()

# --- config.h (AC_CONFIG_HEADERS) ---
#configure_file(
#  ${CMAKE_SOURCE_DIR}/config.h.in
#  ${CMAKE_BINARY_DIR}/config.h
#)
#add_compile_definitions(HAVE_CONFIG_H=1)
#include_directories(${CMAKE_BINARY_DIR})

# --- Dependencies ---
find_library(M_LIB m REQUIRED)
find_library(CONFUSE_LIB confuse REQUIRED)

# GSL: try find_package first, fallback to find_library
find_package(GSL REQUIRED)
if (GSL_FOUND)
  set(GSL_LIBS GSL::gsl GSL::gslcblas)
else()
  find_library(GSL_LIB gsl REQUIRED)
  find_library(GSLCBLAS_LIB gslcblas REQUIRED)
  set(GSL_LIBS ${GSL_LIB} ${GSLCBLAS_LIB})
endif()

find_package(ZLIB REQUIRED)
#find_package(CURL REQUIRED)

# libxml2: prefer CMake package, fallback to pkg-config
find_package(LibXml2 REQUIRED)
if (LibXml2_FOUND)
  set(LIBXML2_TARGET LibXml2::LibXml2)
  set(LIBXML2_INCLUDE_DIRS "")
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(LIBXML2 REQUIRED libxml-2.0)
  set(LIBXML2_TARGET ${LIBXML2_LIBRARIES})
  set(LIBXML2_INCLUDE_DIRS ${LIBXML2_INCLUDE_DIRS})
endif()

# Make common link list available to subdirs
set(TRISURF_COMMON_LIBS
  ${CONFUSE_LIB}
  ${M_LIB}
  ${GSL_LIBS}
  ZLIB::ZLIB
  #  CURL::libcurl
  ${LIBXML2_TARGET}
)

set(TRISURF_COMMON_INCLUDES
  ${LIBXML2_INCLUDE_DIRS}
)

add_subdirectory(src)

