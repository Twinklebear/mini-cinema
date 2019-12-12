#include "loader.h"
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>
#include "json.hpp"
#include "util.h"

Camera::Camera(const vec3f &pos, const vec3f &dir, const vec3f &up)
    : pos(pos), dir(dir), up(up)
{
}

bool compute_divisor(int x, int &divisor)
{
    const int upper = std::sqrt(x);
    for (int i = 2; i <= upper; ++i) {
        if (x % i == 0) {
            divisor = i;
            return true;
        }
    }
    return false;
}

vec3i compute_grid(int num)
{
    vec3i grid(1);
    int axis = 0;
    int divisor = 0;
    while (compute_divisor(num, divisor)) {
        grid[axis] *= divisor;
        num /= divisor;
        axis = (axis + 1) % 3;
    }
    if (num != 1) {
        grid[axis] *= num;
    }
    return grid;
}
std::array<int, 3> compute_ghost_faces(const vec3i &brick_id, const vec3i &grid)
{
    std::array<int, 3> faces = {NEITHER_FACE, NEITHER_FACE, NEITHER_FACE};
    for (size_t i = 0; i < 3; ++i) {
        if (brick_id[i] < grid[i] - 1) {
            faces[i] |= POS_FACE;
        }
        if (brick_id[i] > 0) {
            faces[i] |= NEG_FACE;
        }
    }
    return faces;
}

VolumeBrick load_volume_brick(const json &json, const int mpi_rank, const int mpi_size)
{
    VolumeBrick brick;

    const std::string volume_file = json["volume"].get<std::string>();
    const vec3i volume_dims = get_vec<int, 3>(json["dimensions"]);
    const vec3i grid = compute_grid(mpi_size);
    const vec3i brick_id(
        mpi_rank % grid.x, (mpi_rank / grid.x) % grid.y, mpi_rank / (grid.x * grid.y));

    const vec3i brick_volume_dims = volume_dims / grid;

    const vec3f brick_lower = brick_id * brick_volume_dims;
    const vec3f brick_upper = brick_id * brick_volume_dims + brick_volume_dims;

    brick.bounds = box3f(brick_lower, brick_upper);

    vec3i brick_full_dims = brick_volume_dims;
    vec3i brick_read_offset = brick_lower;
    brick.ghost_bounds = brick.bounds;
    {
        const auto ghost_faces = compute_ghost_faces(brick_id, grid);
        for (size_t i = 0; i < 3; ++i) {
            if (ghost_faces[i] & NEG_FACE) {
                brick_full_dims[i] += 1;
                brick.ghost_bounds.lower[i] -= 1.f;  // todo: base on grid spacing
                brick_read_offset[i] -= 1;
            }
            if (ghost_faces[i] & POS_FACE) {
                brick_full_dims[i] += 1;
                brick.ghost_bounds.upper[i] += 1.f;  // todo: base on grid spacing
            }
        }
    }

    brick.brick = cpp::Volume("structured_regular");
    brick.brick.setParam("dimensions", brick_full_dims);
    brick.brick.setParam("gridOrigin", brick.ghost_bounds.lower);
    brick.brick.setParam("gridSpacing", get_vec<int, 3>(json["grid_spacing"]));

    // Load the sub-bricks using MPI I/O
    {
        size_t voxel_size = 0;
        MPI_Datatype voxel_type;
        const std::string voxel_type_string = json["data_type"].get<std::string>();
        if (voxel_type_string == "uint8") {
            voxel_type = MPI_UNSIGNED_CHAR;
            voxel_size = 1;
        } else if (voxel_type_string == "uint16") {
            voxel_type = MPI_UNSIGNED_SHORT;
            voxel_size = 2;
        } else if (voxel_type_string == "float32") {
            voxel_type = MPI_FLOAT;
            voxel_size = 4;
        } else if (voxel_type_string == "float64") {
            voxel_type = MPI_DOUBLE;
            voxel_size = 8;
        } else {
            throw std::runtime_error("Unrecognized voxel type " + voxel_type_string);
        }

        const size_t n_voxels =
            size_t(brick_full_dims.x) * size_t(brick_full_dims.y) * size_t(brick_full_dims.z);
        brick.voxel_data = std::make_shared<std::vector<uint8_t>>(n_voxels * voxel_size, 0);
        MPI_Datatype brick_type;
        MPI_Type_create_subarray(3,
                                 &volume_dims.x,
                                 &brick_full_dims.x,
                                 &brick_read_offset.x,
                                 MPI_ORDER_FORTRAN,
                                 voxel_type,
                                 &brick_type);
        MPI_Type_commit(&brick_type);

        MPI_File file_handle;
        MPI_File_open(
            MPI_COMM_WORLD, volume_file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handle);
        MPI_File_set_view(file_handle, 0, voxel_type, brick_type, "native", MPI_INFO_NULL);
        MPI_File_read_all(
            file_handle, brick.voxel_data->data(), n_voxels, voxel_type, MPI_STATUS_IGNORE);
        MPI_File_close(&file_handle);
        MPI_Type_free(&brick_type);

        brick.brick.setParam("voxelData", cpp::Data(*brick.voxel_data, true));
    }

    // Set the clipping box of the volume to clip off the ghost voxels
    brick.brick.setParam("volumeClippingBoxLower", brick.bounds.lower);
    brick.brick.setParam("volumeClippingBoxUpper", brick.bounds.upper);
    brick.brick.commit();

    brick.model = cpp::VolumetricModel(brick.brick);
    cpp::TransferFunction tfn("piecewise_linear");
    std::vector<vec3f> colors = {vec3f(0.f, 0.f, 1.f), vec3f(1.f, 0.f, 0.f)};
    std::vector<float> opacities = {0.0001f, 1.f};

    tfn.setParam("color", cpp::Data(colors));
    tfn.setParam("opacity", cpp::Data(opacities));
    // color the bricks by their rank, we pad the range out a bit to keep
    // any brick from being completely transparent
    vec2f value_range = vec2f(0, 256);
    tfn.setParam("valueRange", value_range);
    tfn.commit();
    brick.model.setParam("transferFunction", tfn);
    brick.model.setParam("samplingRate", 0.5f);
    brick.model.commit();

    brick.group = cpp::Group();
    brick.group.setParam("volume", cpp::Data(brick.model));
    brick.group.commit();

    brick.instance = cpp::Instance(brick.group);
    brick.instance.commit();

    return brick;
}

std::vector<Camera> load_cameras(const std::vector<json> &camera_set,
                                 const box3f &world_bounds)
{
    std::vector<Camera> cameras;
    for (size_t i = 0; i < camera_set.size(); ++i) {
        const auto &c = camera_set[i];
        if (c.find("orbit") != c.end()) {
            const float orbit_radius = length(world_bounds.size()) * 0.75f;
            std::cout << "orbit radius: " << orbit_radius
                      << ", points: " << c["orbit"].get<int>() << "\n";
            auto orbit_points = generate_fibonacci_sphere(c["orbit"].get<int>(), orbit_radius);
            for (const auto &p : orbit_points) {
                cameras.emplace_back(p + world_bounds.center(), normalize(-p), vec3f(0, 1, 0));
            }
        } else {
            cameras.emplace_back(get_vec<float, 3>(c["pos"]),
                                 get_vec<float, 3>(c["dir"]),
                                 get_vec<float, 3>(c["up"]));
        }
    }
    return cameras;
}

