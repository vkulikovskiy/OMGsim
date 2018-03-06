# (1)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(OMGsim)


# check for x64 or x86
get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)
if (${LIB64} STREQUAL "TRUE")
    set(LIBSUFFIX 64)
else()
    set(LIBSUFFIX "")
endif()

# (2) ROOT
find_package(ROOT)   #this does not work for old ROOT installations on Lyon
if (ROOT_FOUND)
    message("ROOT found using find_package method. Using in cmake ${ROOT_USE_FILE}")
    include(${ROOT_USE_FILE})
else (ROOT_FOUND)
    message("ROOT was not found using find_package method... falling back to FindROOT.cmake")
    include (FindROOT.cmake)       #We stored locally this file slighly modified to search ROOT
endif(ROOT_FOUND)
message("We use ROOT_INCLUDE_DIRS which is ${ROOT_INCLUDE_DIRS} instead ROOT_INCLUDE_DIR is ${ROOT_INCLUDE_DIR}") 

# (3) Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# (4) organize project libraries
file(GLOB sources */src/[a-zA-Z0-9]*.cc)
file(GLOB headers */include/)

# get rid of test dirs
foreach(adir IN LISTS sources)
  if (adir MATCHES .*/tests/+src.*)
    list (REMOVE_ITEM sources ${adir})
  endif ()
endforeach ()

foreach(adir IN LISTS headers)
  if (adir MATCHES .*/tests/+include/)
    list (REMOVE_ITEM headers ${adir})
  endif ()
endforeach ()

file(GLOB test_sources tests/src/[a-zA-Z0-9]*.cc)
file(GLOB test_headers tests/include/)

#message (${sources})
#message (${test_sources})

file(GLOB PROGRAMS programs/[a-zA-Z1-9]*.cxx)
file(GLOB TESTS tests/[a-zA-Z1-9]*.cxx)

#bstr: set include directories, include ROOT
include_directories (${headers} ${ROOT_INCLUDE_DIRS})
if (DEFINED ENV{X11HEADERS_PATH})
    include_directories (/opt/X11/include)
endif()
    

#bstr: set the libonfig include dir and create linker command. /usr/local is already included,
#bstr: so there is no need to include this separately. Otherwise specify below
#include_directories(/usr/local/include)

#On Lyon libconfig++ libraries are installed in private/project folder, defined in LIBCONFIG_ROOT
if(DEFINED ENV{LIBCONFIG_ROOT})
  message("LIBCONFIG_ROOT is defined to $ENV{LIBCONFIG_ROOT}")
  SET(LIBCONFLIB "-L/usr/local/lib -L$ENV{LIBCONFIG_ROOT}/lib -lconfig++ -lconfig")
else(DEFINED ENV{LIBCONFIG_ROOT})
  message("LIBCONFIG_ROOT is not defined $ENV{LIBCONFIG_ROOT}")
  SET(LIBCONFLIB "-L/usr/local/lib -lconfig++ -lconfig")
endif(DEFINED ENV{LIBCONFIG_ROOT})


# (5) set the needed cflags
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -Wno-variadic-macros -Wshadow -pipe -std=c++11 -DPROJECT_SOURCE_DIR=\\\"${PROJECT_SOURCE_DIR}\\\"")


if(USE_FISSION_NEW)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/${CMAKE_SOURCE_DIR}/FissionLib/fission_v1.9.1/lib -lFission")
  include_directories(${CMAKE_SOURCE_DIR}/FissionLib/fission_v1.9.1/include)
  add_definitions(-DFISSION_NEW)
endif()


if (DEBUG)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
else ()
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif ()

if (FULL_OPTIMIZATION)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
endif (FULL_OPTIMIZATION)

if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.8.9")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=auto")
endif(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.8.9")


SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__CXX_MIN_VERSION=${CMAKE_CXX_COMPILER_VERSION}")

# (6) finish with all the binaries
if(ENABLE_TESTING)
  add_library(libOMG SHARED ${sources})
  target_link_libraries(libOMG ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${LIBCONFLIB})
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_TESTING=1")
  foreach(aprog IN LISTS TESTS)
    get_filename_component(aprogname ${aprog} NAME_WE)
    add_executable(${aprogname} ${aprog} ${headers} ${test_headers} ${test_sources})
    target_link_libraries(${aprogname} libOMG ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${LIBCONFLIB})
  endforeach ()
else(ENABLE_TESTING)
  add_library(libOMG SHARED ${sources})
  target_link_libraries(libOMG ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${LIBCONFLIB})
  foreach(aprog IN LISTS PROGRAMS)
    get_filename_component(aprogname ${aprog} NAME_WE)
    add_executable(${aprogname} ${aprog} ${headers})
    target_link_libraries(${aprogname} libOMG ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${LIBCONFLIB})
  endforeach ()
endif(ENABLE_TESTING)

file(COPY ${PROJECT_SOURCE_DIR}/scripts DESTINATION ${PROJECT_BINARY_DIR})
