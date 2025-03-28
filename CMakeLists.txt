########################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( MarlinACTS )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 3 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### BUILD OPTIONS ###########################################################

set(PACKAGE_INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(PACKAGE_INSTALL_INC_DIR include CACHE PATH "Installation directory for headers")

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDCore )
INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}  )

FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
LINK_LIBRARIES( ${LCIO_LIBRARIES} )

FIND_PACKAGE( ROOT REQUIRED)
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )

FIND_PACKAGE( Acts REQUIRED COMPONENTS Core PluginJson PluginTGeo PluginDD4hep )

INCLUDE(GNUInstallDirs)

### DATA ###########################################################
ADD_SUBDIRECTORY(data)

### LIBRARY #################################################################
ADD_DEFINITIONS( "-Wno-long-long -Wno-effc++ -Wno-unused-parameter -Wno-narrowing" )

# sources
ADD_LIBRARY(MarlinACTS SHARED src/ACTSBaseTracker.cpp
                              src/ACTSCKFBaseTracker.cpp
                              src/ACTSCKFSeededTracker.cpp
                              src/ACTSCKFTruthTracker.cpp
                              src/CollectionHandler.cpp
                              src/GeometryIdMappingTool.cpp
                              src/GeometryIdSelector.cpp
                              src/Helpers.cpp)

TARGET_INCLUDE_DIRECTORIES(MarlinACTS PRIVATE
                                      ./include
                                      ${CMAKE_CURRENT_BINARY_DIR}
                                      ${Marlin_INCLUDE_DIRS}
                                      ${ROOT_INCLUDE_DIRS}
                                      ${DD4hep_INCLUDE_DIRS}
                                      ${LCIO_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(MarlinACTS ${ROOT_LIBRARIES}
                                 DD4hep::DDCore
                                 ActsCore
                                 ActsPluginTGeo
                                 ActsPluginJson
                                 ActsPluginDD4hep)

INSTALL( TARGETS MarlinACTS DESTINATION ${PACKAGE_INSTALL_LIB_DIR} )
INSTALL( DIRECTORY MarlinACTS DESTINATION ${PACKAGE_INSTALL_INC_DIR}
         FILES_MATCHING PATTERN "*.hxx" )


