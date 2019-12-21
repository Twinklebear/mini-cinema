#include <chrono>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>
#include "json.hpp"
#include "stb_image.h"
#include "util.h"

#ifdef VTK_FOUND
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkFlyingEdges3D.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#endif

#include "loader.h"

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

VolumeBrick load_volume_brick(json &config,
                              const is::SimState &region,
                              const int mpi_rank,
                              const int mpi_size)
{
    using namespace std::chrono;
    VolumeBrick brick;

    const std::string field_name = config["field"].get<std::string>();
    const auto it = region.fields.find(field_name);
    if (it == region.fields.end()) {
        std::cerr << "[error]: Requested field " << field_name
                  << " was not found in the simulation\n";
        throw std::runtime_error("[error]: Field not found");
    }
    const is::Field &field = it->second;
    brick.voxel_data = field.array;
    brick.dims = vec3i(field.dims[0], field.dims[1], field.dims[2]);
    brick.bounds = box3f(vec3f(region.local.min.x, region.local.min.y, region.local.min.z),
                         vec3f(region.local.max.x, region.local.max.y, region.local.max.z));
    const vec3f spacing = get_vec<int, 3>(config["spacing"]);
    brick.full_dims = brick.dims;
    // TODO: This is assuming the sim isn't giving us ghost zone data
    brick.ghost_bounds = brick.bounds;

    brick.brick = cpp::Volume("structured_regular");
    brick.brick.setParam("dimensions", brick.full_dims);
    brick.brick.setParam("gridSpacing", spacing);

    const size_t n_voxels = brick.dims.long_product();
    cpp::Data osp_data;
    if (field.dataType == UINT8) {
        osp_data =
            cpp::Data(n_voxels, static_cast<const uint8_t *>(brick.voxel_data->data()), true);
        config["type"] = "uint8";
    } else if (field.dataType == FLOAT) {
        osp_data =
            cpp::Data(n_voxels, static_cast<const float *>(brick.voxel_data->data()), true);
        config["type"] = "float32";
    } else if (field.dataType == DOUBLE) {
        osp_data =
            cpp::Data(n_voxels, static_cast<const double *>(brick.voxel_data->data()), true);
        config["type"] = "float64";
    } else {
        std::cerr << "[error]: Unsupported voxel type\n";
        throw std::runtime_error("[error]: Unsupported voxel type");
    }
    brick.brick.setParam("data", osp_data);

    // If the value range wasn't provided, compute it
    if (config.find("value_range") == config.end()) {
        vec2f value_range;
        if (field.dataType == UINT8) {
            value_range = compute_value_range(
                reinterpret_cast<uint8_t *>(brick.voxel_data->data()), n_voxels);
        } else if (field.dataType == FLOAT) {
            value_range = compute_value_range(
                reinterpret_cast<float *>(brick.voxel_data->data()), n_voxels);
        } else if (field.dataType == DOUBLE) {
            value_range = compute_value_range(
                reinterpret_cast<double *>(brick.voxel_data->data()), n_voxels);
        } else {
            std::cerr << "[error]: Unsupported voxel type\n";
            throw std::runtime_error("[error]: Unsupported voxel type");
        }

        vec2f global_value_range;
        MPI_Allreduce(
            &value_range.x, &global_value_range.x, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(
            &value_range.y, &global_value_range.y, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

        if (mpi_rank == 0) {
            std::cout << "Computed value range: " << global_value_range << "\n";
        }
        config["value_range"] = {global_value_range.x, global_value_range.y};
    }

    // Set the clipping box of the volume to clip off the ghost voxels
    brick.brick.setParam("volumeClippingBoxLower", brick.bounds.lower);
    brick.brick.setParam("volumeClippingBoxUpper", brick.bounds.upper);
    brick.brick.commit();
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

std::vector<cpp::TransferFunction> load_colormaps(const std::vector<std::string> &files,
                                                  const vec2f &value_range)
{
    std::vector<cpp::TransferFunction> colormaps;
    for (const auto &f : files) {
        cpp::TransferFunction tfn("piecewise_linear");
        int x, y, n;
        uint8_t *data = stbi_load(f.c_str(), &x, &y, &n, 4);
        if (!data) {
            std::cerr << "[error]: failed to load image: " << f << "\n" << std::flush;
            throw std::runtime_error("Failed to load " + f);
        }

        std::vector<vec3f> colors;
        std::vector<float> opacities;
        for (int i = 0; i < x; ++i) {
            colors.emplace_back(
                data[i * 4] / 255.f, data[i * 4 + 1] / 255.f, data[i * 4 + 2] / 255.f);
            // If no alpha in the image, generate a linear ramp
            if (n == 3) {
                opacities.emplace_back(static_cast<float>(i) / x);
            } else {
                opacities.emplace_back(data[i * 4 + 3] / 255.f);
            }
        }
        stbi_image_free(data);

        std::cout << "value range: " << value_range << "\n";
        tfn.setParam("color", cpp::Data(colors));
        tfn.setParam("opacity", cpp::Data(opacities));
        tfn.setParam("valueRange", value_range);
        tfn.commit();
        colormaps.push_back(tfn);
    }
    return colormaps;
}

std::vector<cpp::Geometry> extract_isosurfaces(const json &config,
                                               const VolumeBrick &brick,
                                               const int mpi_rank,
                                               const vec2f &value_range)
{
    using namespace std::chrono;

    std::vector<cpp::Geometry> isosurfaces;
#ifdef VTK_FOUND
    const std::string voxel_type_string = config["type"].get<std::string>();
    vtkSmartPointer<vtkDataArray> data_array = nullptr;
    if (voxel_type_string == "uint8") {
        auto arr = vtkSmartPointer<vtkUnsignedCharArray>::New();
        arr->SetArray(
            const_cast<uint8_t *>(static_cast<const uint8_t *>(brick.voxel_data->data())),
            brick.voxel_data->size(),
            1);
        data_array = arr;
    } else if (voxel_type_string == "float32") {
        auto arr = vtkSmartPointer<vtkFloatArray>::New();
        arr->SetArray(
            const_cast<float *>(static_cast<const float *>(brick.voxel_data->data())),
            brick.voxel_data->size(),
            1);
        data_array = arr;
    } else if (voxel_type_string == "float64") {
        auto arr = vtkSmartPointer<vtkDoubleArray>::New();
        arr->SetArray(
            const_cast<double *>(static_cast<const double *>(brick.voxel_data->data())),
            brick.voxel_data->size(),
            1);
        data_array = arr;
    } else {
        throw std::runtime_error("Unrecognized voxel type " + voxel_type_string);
    }

    vtkSmartPointer<vtkImageData> img_data = vtkSmartPointer<vtkImageData>::New();
    img_data->SetDimensions(brick.full_dims.x, brick.full_dims.y, brick.full_dims.z);
    img_data->GetPointData()->SetScalars(data_array);

    // For simulations isovalues may be best treated as [0, 1] picking a t value along
    // the range of the simulation's data values. Otherwise we have to know the value range
    // ahead of time, which we may not know.
    auto isovals = config["isovalue"].get<std::vector<float>>();
    for (const auto &v : isovals) {
        const float val = lerp(v, value_range.x, value_range.y);
        auto start = high_resolution_clock::now();
        vtkSmartPointer<vtkFlyingEdges3D> fedges = vtkSmartPointer<vtkFlyingEdges3D>::New();
        fedges->SetInputData(img_data);
        fedges->SetNumberOfContours(1);
        fedges->SetValue(0, val);
        fedges->SetComputeNormals(false);
        fedges->Update();
        vtkPolyData *isosurf = fedges->GetOutput();

        std::vector<vec3f> vertices;
        std::vector<vec3ui> indices;
        vertices.reserve(isosurf->GetNumberOfCells());
        indices.reserve(isosurf->GetNumberOfCells());
        for (size_t i = 0; i < isosurf->GetNumberOfCells(); ++i) {
            vtkTriangle *tri = dynamic_cast<vtkTriangle *>(isosurf->GetCell(i));
            if (tri->ComputeArea() == 0.0) {
                continue;
            }
            vec3ui tids;
            for (size_t v = 0; v < 3; ++v) {
                const double *pt = isosurf->GetPoint(tri->GetPointId(v));
                const vec3f vert(pt[0], pt[1], pt[2]);
                tids[v] = vertices.size();
                vertices.push_back(vert);
            }
            indices.push_back(tids);
        }
        auto end = high_resolution_clock::now();

        std::cout << "(rank " << mpi_rank << "): Isosurface at " << val << " has "
                  << indices.size() << " triangles, extracted in "
                  << duration_cast<milliseconds>(end - start).count() << "ms\n";
        if (indices.empty()) {
            continue;
        }

        cpp::Geometry mesh("mesh");
        mesh.setParam("vertex.position", cpp::Data(vertices));
        mesh.setParam("index", cpp::Data(indices));
        mesh.commit();
        isosurfaces.push_back(mesh);
    }
#else
    std::cerr << "[warning]: Scene requested isosurface geometry, but app was not compiled "
                 "with VTK to support explicit isosurfaces.\n";
#endif
    return isosurfaces;
}

