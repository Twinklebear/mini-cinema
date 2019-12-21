#pragma once

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>
#include "json.hpp"

#include <libIS/is_client.h>

using namespace ospray;
using namespace ospcommon;
using namespace ospcommon::math;
using json = nlohmann::json;

struct VolumeBrick {
    // the volume data itself
    cpp::Volume brick;
    // the bounds of the owned portion of data
    box3f bounds;
    // the full bounds of the owned portion + ghost voxels
    box3f ghost_bounds;

    vec3i dims;
    // full dims includes ghost voxels
    vec3i full_dims;

    vec2f value_range;

    std::shared_ptr<is::Array> voxel_data;
};

struct Camera {
    vec3f pos;
    vec3f dir;
    vec3f up;

    Camera(const vec3f &pos, const vec3f &dir, const vec3f &up);
};

struct Isosurface {
    size_t n_triangles = 0;
    cpp::Geometry geometry;
};

bool compute_divisor(int x, int &divisor);

/* Compute an X x Y x Z grid to have 'num' grid cells,
 * only gives a nice grid for numbers with even factors since
 * we don't search for factors of the number, we just try dividing by two
 */
vec3i compute_grid(int num);

enum GhostFace { NEITHER_FACE = 0, POS_FACE = 1, NEG_FACE = 2 };

/* Compute which faces of this brick we need to specify ghost voxels
 * for to have correct interpolation at brick boundaries.
 */
std::array<int, 3> compute_ghost_faces(const vec3i &brick_id, const vec3i &grid);

VolumeBrick load_volume_brick(json &config,
                              const is::SimState &region,
                              const int mpi_rank,
                              const int mpi_size);

std::vector<Camera> load_cameras(const std::vector<json> &camera_set,
                                 const box3f &world_bounds);

std::vector<cpp::TransferFunction> load_colormaps(const std::vector<std::string> &files,
                                                  const vec2f &value_range);

std::vector<Isosurface> extract_isosurfaces(const json &config,
                                            const VolumeBrick &brick,
                                            const int mpi_rank,
                                            const vec2f &value_range);
