#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(DoseWater)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/src
                    ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
add_executable(DoseWater main.cc ${sources} ${headers})
target_link_libraries(DoseWater ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
set(SCRIPTS
  init_vis.mac
  vis.mac
  run.mac
  )

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS DoseWater DESTINATION bin)
