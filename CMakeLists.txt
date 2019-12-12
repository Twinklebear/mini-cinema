cmake_minimum_required(VERSION 3.5)
project(minicinema)

if (NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

find_package(Threads REQUIRED)
find_package(ospray 2.0 REQUIRED)
find_package(MPI REQUIRED)
find_package(VTK)

add_executable(mini_cinema main.cpp loader.cpp util.cpp)

set_target_properties(mini_cinema PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON)


target_link_libraries(mini_cinema PUBLIC
    ospray::ospray
    ospcommon::ospcommon
    MPI::MPI_CXX)

if (${VTK_FOUND})
    target_compile_definitions(mini_cinema PUBLIC -DVTK_FOUND=1)
    target_include_directories(mini_cinema PUBLIC
        ${VTK_INCLUDE_DIRS})

    target_link_libraries(mini_cinema PUBLIC
        ${VTK_LIBRARIES})
else()
    message(WARNING "VTK not found, but is required for explicit isosurfaces. "
        "Isosurface support will be disabled")
endif()

