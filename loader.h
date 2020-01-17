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

enum GhostFace { NEITHER_FACE = 0, POS_FACE = 1, NEG_FACE = 2 };

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

using DataArray = std::shared_ptr<std::vector<uint8_t>>;

struct GhostData {
    vec3i ghost_faces = vec3i(NEITHER_FACE);
    // Faces are -/+x, -/+y, -/+z
    std::array<DataArray, 6> faces = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    // Note: just doing faces to start, edges/corners likely not visible
    // std::array<std::shared_ptr<std::vector<uint8_t>>, 12> edges;
    // std::array<std::shared_ptr<std::vector<uint8_t>>, 8> corners;

    GhostData() = default;
    GhostData(const vec3i &brick_id, const vec3i &grid, const is::Field &field);
};

struct ParticleBrick {
    // Particle geometry
    cpp::Geometry geom;

    size_t num_particles;
    // bounds of the owned portion of the particles
    box3f bounds;
    // full bounds of the owned portion + ghost zones
    box3f ghost_bounds;
    // The particle data array
    std::shared_ptr<is::Array> data;

    ParticleBrick() = default;
    ParticleBrick(const is::SimState &region);
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

/* Compute which faces of this brick we need to specify ghost voxels
 * for to have correct interpolation at brick boundaries.
 */
vec3i compute_ghost_faces(const vec3i &brick_id, const vec3i &grid);

std::vector<VolumeBrick> load_volume_bricks(json &config,
                                            const std::vector<is::SimState> &regions,
                                            const int mpi_rank,
                                            const int mpi_size);

std::vector<Camera> load_cameras(const std::vector<json> &camera_set,
                                 const box3f &world_bounds);

std::vector<cpp::TransferFunction> load_colormaps(const std::vector<std::string> &files,
                                                  const vec2f &value_range);

// Load an RGB8 texture
cpp::Texture load_texture(const std::string &file);

std::vector<Isosurface> extract_isosurfaces(const json &config,
                                            const VolumeBrick &brick,
                                            const int mpi_rank,
                                            const vec2f &value_range);

// Load the set of particle bricks for the regions. Returns an empty vector if
// no particles in any region
std::vector<ParticleBrick> load_particle_bricks(const std::vector<is::SimState> &regions,
                                                const float radius);
