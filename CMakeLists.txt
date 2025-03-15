cmake_minimum_required(VERSION 3.5.1)
project(frbsearchlib VERSION 1.0.1)
include(GNUInstallDirs)

file(GLOB frbsearchlib_sources "src/*.cpp")
file(GLOB frbsearchlib_headers "src/*.h")
add_library(frbsearchlib SHARED ${frbsearchlib_sources})

add_definitions(-D_UNIX -D_ENABLE_PROFILER_ -D_ENABLE_TRACING_ -D_ENABLE_VISUALIZATION_ -D_MONTE_CARLO_ -D_NO_ROOT_ -I${CMAKE_SOURCE_DIR}/src)

set_target_properties(frbsearchlib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${frbsearchlib_headers}")
 
configure_file(frbsearchlib.pc.in frbsearchlib.pc @ONLY)


install(TARGETS frbsearchlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES ${CMAKE_BINARY_DIR}/frbsearchlib.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


find_library(CFITSIO_LIB cfitsio HINTS ENV LD_LIBRARY_PATH)
find_library(FFTW3_LIB fftw3 HINTS ENV FFTW_LIB REQUIRED)
find_library(LIBNOVA_LIB nova PATHS ENV LD_LIBRARY_PATH)


# required to properly link the dynamic library :
target_link_libraries(frbsearchlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

add_executable(create_dynaspec apps/main_dynaspec.cpp)
target_link_libraries(create_dynaspec frbsearchlib msfitslib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

add_executable(dynaspec_search apps/main_test_dynaspec.cpp)
target_link_libraries(dynaspec_search frbsearchlib msfitslib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

add_executable(calcstat apps/calcstat.cpp)
target_link_libraries(calcstat frbsearchlib msfitslib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

add_executable(dedisperse_fits apps/main_dedisperse_fits.cpp)
target_link_libraries(dedisperse_fits frbsearchlib msfitslib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)



# INSTALLATION:
# install(TARGETS noise_mapper frbsearchlib calcfits_bg dump_lc avg_images homeopatic_avg image_differencer ux2sid_file ux2sid sid2ux radec2azh date2date print_sun flip_vertical flip_horizontal create_dynaspec dynaspec_search RUNTIME DESTINATION bin)
