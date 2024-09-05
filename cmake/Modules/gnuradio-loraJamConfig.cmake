find_package(PkgConfig)

PKG_CHECK_MODULES(PC_GR_LORAJAM gnuradio-loraJam)

FIND_PATH(
    GR_LORAJAM_INCLUDE_DIRS
    NAMES gnuradio/loraJam/api.h
    HINTS $ENV{LORAJAM_DIR}/include
        ${PC_LORAJAM_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GR_LORAJAM_LIBRARIES
    NAMES gnuradio-loraJam
    HINTS $ENV{LORAJAM_DIR}/lib
        ${PC_LORAJAM_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-loraJamTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GR_LORAJAM DEFAULT_MSG GR_LORAJAM_LIBRARIES GR_LORAJAM_INCLUDE_DIRS)
MARK_AS_ADVANCED(GR_LORAJAM_LIBRARIES GR_LORAJAM_INCLUDE_DIRS)
