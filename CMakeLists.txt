cmake_minimum_required(VERSION 3.16)
project(android-gpos)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    message(STATUS "optional:-std=c++11")
endif(CMAKE_COMPILER_IS_GNUCXX)

include_directories(${PROJECT_BINARY_DIR})

aux_source_directory(src SRC_ANDROID)

if(WIN32)
    set(LIBS Winmm)
else()
    set(LIBS m pthread)
endif()

add_library(agpos ${SRC_ANDROID})
target_link_libraries(agpos ${LIBS} gfortran quadmath)

add_executable(test_agpos test_agpos.c)
target_link_libraries(test_agpos agpos)







