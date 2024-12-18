# --------------------------------------------------------------------
# Installation directories
# --------------------------------------------------------------------
include(GNUInstallDirs)

# Define installation directories
set(INSTALL_BIN_DIR
  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
  CACHE PATH "Installation directory for executables"
)
set(INSTALL_LIB_DIR
  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  CACHE PATH "Installation directory for libraries"
)

# --------------------------------------------------------------------
# Installation of dependencies
# --------------------------------------------------------------------
# Function to install shared library dependencies
function(install_shared_libs TARGET)
  if(NOT CMAKE_CROSSCOMPILING)
    get_target_property(TARGET_TYPE ${TARGET} TYPE)
    if(TARGET_TYPE STREQUAL "EXECUTABLE")
      # Get runtime dependencies
      include(GetPrerequisites)
      get_prerequisites(${TARGET} DEPENDENCIES 1 1 "" "")

      # Install each dependency
      foreach(DEPENDENCY ${DEPENDENCIES})
        if(EXISTS "${DEPENDENCY}")
          file(COPY ${DEPENDENCY} DESTINATION ${INSTALL_LIB_DIR} FOLLOW_SYMLINK_CHAIN)
        endif()
      endforeach()
    endif()
  endif()
endfunction()

# --------------------------------------------------------------------
# Boost libraries installation
# --------------------------------------------------------------------
# Get Boost library locations
get_target_property(BOOST_DATEIME_LIB Boost::date_time IMPORTED_LOCATION_RELEASE)
get_target_property(BOOST_FILESYSTEM_LIB Boost::filesystem IMPORTED_LOCATION_RELEASE)
get_target_property(BOOST_PROGRAM_OPTIONS_LIB Boost::program_options IMPORTED_LOCATION_RELEASE)
get_target_property(BOOST_SYSTEM_LIB Boost::system IMPORTED_LOCATION_RELEASE)
get_target_property(BOOST_REGEX_LIB Boost::regex IMPORTED_LOCATION_RELEASE)
get_target_property(BOOST_TIMER_LIB Boost::timer IMPORTED_LOCATION_RELEASE)
get_target_property(BOOST_CHRONO_LIB Boost::chrono IMPORTED_LOCATION_RELEASE)

# Install Boost libraries
install(FILES 
  ${BOOST_DATEIME_LIB}
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_PROGRAM_OPTIONS_LIB}
  ${BOOST_SYSTEM_LIB}
  ${BOOST_REGEX_LIB}
  ${BOOST_TIMER_LIB}
  ${BOOST_CHRONO_LIB}
  DESTINATION ${INSTALL_LIB_DIR}
  OPTIONAL
)
