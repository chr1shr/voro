cmake_minimum_required(VERSION 3.10)

project(voro++ VERSION 0.4.6 LANGUAGES CXX)
set(SOVERSION "0")

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

########################################################################
# User input options                                                   #
########################################################################
option(VORO_BUILD_SHARED_LIBS "Build shared libs" ON)
include(GNUInstallDirs)

########################################################################
#Find external packages
########################################################################
find_package(Doxygen)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

######################################
# Include the following subdirectory # 
######################################

file(GLOB VORO_SOURCES src/*.cc)
file(GLOB NOT_VORO_SOURCES src/v_base_wl.cc src/cmd_line.cc src/voro++.cc src/v_base_wl_2d.cc src/v_base_wl_3d.cc)
list(REMOVE_ITEM VORO_SOURCES ${NOT_VORO_SOURCES})
add_library(voro++ ${VORO_SOURCES})
# target_compile_options(voro++ PRIVATE -w)
target_link_libraries(voro++ PRIVATE "${OpenMP_CXX_FLAGS}")
set_target_properties(voro++ PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src"
  SOVERSION ${SOVERSION})
install(TARGETS voro++ EXPORT VORO_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
#for voto++.hh
target_include_directories(voro++ PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_executable(cmd_line src/cmd_line.cc)
target_link_libraries(cmd_line PRIVATE voro++)
#cannot have two targets with the same name, so renaming cmd_line
set_target_properties(cmd_line PROPERTIES OUTPUT_NAME voro++ 
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src") 
install(TARGETS cmd_line RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

file(GLOB EXAMPLE_SOURCES examples/*/*.cc)
file(GLOB NOT_EXAMPLE examples/no_release/*.cc)
list(REMOVE_ITEM EXAMPLE_SOURCES ${NOT_EXAMPLE})

foreach(SOURCE ${EXAMPLE_SOURCES})
  string(REGEX REPLACE "^.*/([^/]*)\\.cc$" "\\1" PROGNAME "${SOURCE}")
  if (NOT PROGNAME STREQUAL ellipsoid) #ellipsoid is broken
    string(REGEX REPLACE "^.*/(examples/.*)/${PROGNAME}\\.cc$" "\\1" DIRNAME "${SOURCE}")
    add_executable(${PROGNAME} ${SOURCE})
    # target_compile_options(${PROGNAME} PRIVATE -w)
    target_link_libraries(${PROGNAME} PRIVATE voro++)
    set_target_properties(${PROGNAME} PROPERTIES 
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DIRNAME}" )
  endif()
endforeach(SOURCE)

file(GLOB_RECURSE VORO_HEADERS src/voro++.hh)
install(FILES ${VORO_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_SOURCE_DIR}/man/voro++.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
# no external deps for we can use target file as config file
install(EXPORT VORO_Targets FILE VOROConfig.cmake NAMESPACE VORO:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VORO)
include(CMakePackageConfigHelpers)
write_basic_package_version_file("VOROConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/VOROConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VORO)

if (DOXYGEN_FOUND)
  add_custom_target(doxygen COMMAND ${DOXYGEN_EXECUTABLE} src/Doxyfile 
    COMMENT "Build doxygen documentation")
endif (DOXYGEN_FOUND)
