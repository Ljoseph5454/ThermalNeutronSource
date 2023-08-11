cmake_minimum_required(VERSION 3.10...3.26)
project(B1) #LANGUAGUES CXX)
#set(CMAKE_CXX_STANDARD 11)
find_package(NCrystal REQUIRED COMPONENTS GEANT4BINDINGS)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(DMX DMX.cc ${sources} ${headers})
target_link_libraries(DMX ${Geant4_LIBRARIES})
target_link_libraries(DMX NCrystal::G4NCrystal)

add_custom_target(B1 DEPENDS DMX)

install(TARGETS DMX DESTINATION bin)


