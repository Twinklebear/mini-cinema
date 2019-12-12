cmake_minimum_required(VERSION 3.5)
project(minicinema)

if (NOT WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

find_package(Threads REQUIRED)
find_package(ospray 2.0 REQUIRED)
find_package(MPI REQUIRED)
find_package(VTK REQUIRED)

add_executable(mini_cinema main.cpp loader.cpp util.cpp)

set_target_properties(mini_cinema PROPERTIES
	CXX_STANDARD 14
	CXX_STANDARD_REQUIRED ON)

target_include_directories(mini_cinema PUBLIC
    ${VTK_INCLUDE_DIRS})

target_link_libraries(mini_cinema PUBLIC
    ospray::ospray
    ospcommon::ospcommon
    MPI::MPI_CXX
    ${VTK_LIBRARIES})

