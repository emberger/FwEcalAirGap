# $Id: CMakeLists.txt 100922a 2016-11-03 10:50:34Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(B4c)

set(CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#set(CMAKE_PREFIX_PATH $ENV{ROOTSYS}/cmake)
#find_package(ROOT REQUIRED COMPONENTS RIO Net)


set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)

#find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Matrix Minuit Minuit2)
#find_package(ROOT REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Minuit Minuit2)
set(ROOT_LIBRARIES ${ROOT_LIBRARIES})

#find_package(ROOT REQUIRED COMPONENTS libMinuit libMinuit2)
#set(ROOT_LIBRARIES ${ROOT_LIBRARIES})
#include(${ROOT_USE_FILE})



include_directories(${ROOT_INCLUDE_DIR})
#include_directories(${ROOT_INCLUDE_DIR}/Minuit2)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})



#find_package(Boost COMPONENTS system filesystem REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#------------------------------------------------------
#Find and include Eigen

#list(APPEND CMAKE_MODULE_PATH $ENV{EIGEN})
#find_package(Eigen3 REQUIRED)

#SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
#IF( NOT EIGEN3_INCLUDE_DIR )
#    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
#ENDIF()
#INCLUDE_DIRECTORIES (${EIGEN3_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleB4c exampleB4c.cc ${sources} ${headers})
add_executable(Analysis AnalysisMain.cc ${sources} ${headers})
add_executable(AngularResolution AngularResolution.cc ${sources} ${headers})
add_executable(PlotQuantiles PlotQuantiles.cc ${sources} ${headers})
add_executable(EResAndShowerprofile EResAndShowerprofile.cc ${sources} ${headers})

target_link_libraries(exampleB4c ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )#RootLib)
target_link_libraries(Analysis  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )#RootLib)
target_link_libraries(AngularResolution ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )# ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(PlotQuantiles ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
target_link_libraries(EResAndShowerprofile ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4c. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB4C_SCRIPTS
  exampleB4c.out
  exampleB4.in
  gui.mac
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${EXAMPLEB4C_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
 #add_executable(Analysis AnalysisMain.cc ${sources} ${headers})
   )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleB4c DESTINATION bin)
install(TARGETS AngularResolution DESTINATION bin)
install(TARGETS Analysis DESTINATION bin)
install(TARGETS PlotQuantiles DESTINATION bin)
install(TARGETS EResAndShowerprofile DESTINATION bin)

