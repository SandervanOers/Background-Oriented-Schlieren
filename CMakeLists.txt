SET(CMAKE_LEGACY_CYGWIN_WIN32 0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
cmake_minimum_required (VERSION 2.8)
project (DIC)
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable (DIC DisplayImage.cpp)
set (LIBSRC DisplayImage.hpp PointsWithValue.hpp 
PositionDirection.hpp pixeltranslation.hpp pixeltranslation.cpp 
nonlineariteration.hpp 
nonlineariteration.cpp calculateN.hpp calculateN.cpp InputOut.hpp InputOut.cpp coeff.h coeff.c interpol.h interpol.c)
add_library(coeff STATIC ${LIBSRC})
target_link_libraries( DIC coeff ${OpenCV_LIBS} -lpthread)
set(CMAKE_CXX_FLAGS "-g3 -pg -O0 -Wall -Wextra -Wconversion -pedantic")
set(CMAKE_C_FLAGS "-g3 -pg -O0 -Wall -Wextra -Wconversion -pedantic")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -pg")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g3 -pg -O0")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g3 -pg -O0")
