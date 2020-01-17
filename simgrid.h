#pragma once

#include <iostream>
#include <vector>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>

#include <libIS/is_client.h>

using namespace ospcommon::math;

struct SimBrick {
    box3f bounds;
    vec3i id = vec3i(-1);
    int owner = -1;

    SimBrick() = default;
    SimBrick(const box3f &b, int owner);
};

struct SimGrid {
    std::vector<SimBrick> bricks;
    vec3i grid = vec3i(0);

    int brick_owner(const vec3i &id) const;

    const SimBrick &brick_at(const box3f &bounds) const;
};

std::ostream &operator<<(std::ostream &os, const SimGrid &sg);

SimGrid reconstruct_grid(const std::vector<is::SimState> &regions,
                         const int mpi_rank,
                         const int mpi_size);

