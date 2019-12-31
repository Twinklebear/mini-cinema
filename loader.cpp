#include "loader.h"
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

VolumeBrick load_volume_brick(json &config, const int mpi_rank, const int mpi_size)
{
    using namespace std::chrono;
    VolumeBrick brick;

    const std::string volume_file = config["volume"].get<std::string>();
    const vec3i volume_dims = get_vec<int, 3>(config["size"]);
    const vec3f spacing = get_vec<int, 3>(config["spacing"]);
    const vec3i grid = compute_grid(mpi_size);
    const vec3i brick_id(
        mpi_rank % grid.x, (mpi_rank / grid.x) % grid.y, mpi_rank / (grid.x * grid.y));

    brick.dims = volume_dims / grid;

    const vec3f brick_lower = brick_id * brick.dims;
    const vec3f brick_upper = brick_id * brick.dims + brick.dims;

    brick.bounds = box3f(brick_lower, brick_upper);

    brick.full_dims = brick.dims;
    vec3i brick_read_offset = brick_lower;
    brick.ghost_bounds = brick.bounds;
    {
        const auto ghost_faces = compute_ghost_faces(brick_id, grid);
        for (size_t i = 0; i < 3; ++i) {
            if (ghost_faces[i] & NEG_FACE) {
                brick.full_dims[i] += 1;
                brick.ghost_bounds.lower[i] -= spacing[i];
                brick_read_offset[i] -= 1;
            }
            if (ghost_faces[i] & POS_FACE) {
                brick.full_dims[i] += 1;
                brick.ghost_bounds.upper[i] += spacing[i];
            }
        }
    }

    brick.brick = cpp::Volume("structured_regular");
    brick.brick.setParam("dimensions", brick.full_dims);
    brick.brick.setParam("gridSpacing", spacing);

    const std::string voxel_type = config["type"].get<std::string>();
    const size_t n_voxels = brick.full_dims.long_product();

    auto start = high_resolution_clock::now();

    brick.voxel_data = load_raw_volume(
        volume_file, voxel_type, volume_dims, brick.full_dims, brick_read_offset);

    auto end = high_resolution_clock::now();
    if (mpi_rank == 0) {
        std::cout << "Loading volume brick took "
                  << duration_cast<milliseconds>(end - start).count() << "ms\n";
    }

    cpp::Data osp_data;
    if (voxel_type == "uint8") {
        osp_data = cpp::Data(vec3ul(brick.full_dims), brick.voxel_data->data(), true);
    } else if (voxel_type == "uint16") {
        osp_data = cpp::Data(vec3ul(brick.full_dims),
                             reinterpret_cast<uint16_t *>(brick.voxel_data->data()),
                             true);
    } else if (voxel_type == "float32") {
        osp_data = cpp::Data(vec3ul(brick.full_dims),
                             reinterpret_cast<float *>(brick.voxel_data->data()),
                             true);
    } else if (voxel_type == "float64") {
        osp_data = cpp::Data(vec3ul(brick.full_dims),
                             reinterpret_cast<double *>(brick.voxel_data->data()),
                             true);
    } else {
        std::cerr << "[error]: Unsupported voxel type\n";
        throw std::runtime_error("[error]: Unsupported voxel type");
    }
    brick.brick.setParam("data", osp_data);

    // If the value range wasn't provided, compute it
    if (config.find("value_range") == config.end()) {
        start = high_resolution_clock::now();

        vec2f value_range;
        if (voxel_type == "uint8") {
            value_range = compute_value_range(brick.voxel_data->data(), n_voxels);
        } else if (voxel_type == "uint16") {
            value_range = compute_value_range(
                reinterpret_cast<uint16_t *>(brick.voxel_data->data()), n_voxels);
        } else if (voxel_type == "float32") {
            value_range = compute_value_range(
                reinterpret_cast<float *>(brick.voxel_data->data()), n_voxels);
        } else if (voxel_type == "float64") {
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

        end = high_resolution_clock::now();

        if (mpi_rank == 0) {
            std::cout << "Computed value range: " << global_value_range << "\n"
                      << "Value range computation took "
                      << duration_cast<milliseconds>(end - start).count() << "ms\n";
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
                                               const bool isosurface_full_volume)
{
    using namespace std::chrono;

    std::vector<cpp::Geometry> isosurfaces;
#ifdef VTK_FOUND
    std::shared_ptr<std::vector<uint8_t>> voxel_data = brick.voxel_data;
    vec3i brick_dims = brick.full_dims;

    const std::string voxel_type_string = config["type"].get<std::string>();
    if (isosurface_full_volume) {
        const std::string volume_file = config["volume"].get<std::string>();
        const vec3i volume_dims = get_vec<int, 3>(config["size"]);
        brick_dims = volume_dims;
        voxel_data = load_raw_volume(
            volume_file, voxel_type_string, volume_dims, volume_dims, vec3i(0));
    }

    vtkSmartPointer<vtkDataArray> data_array = nullptr;
    if (voxel_type_string == "uint8") {
        auto arr = vtkSmartPointer<vtkUnsignedCharArray>::New();
        arr->SetArray(voxel_data->data(), voxel_data->size(), 1);
        data_array = arr;
    } else if (voxel_type_string == "uint16") {
        auto arr = vtkSmartPointer<vtkUnsignedShortArray>::New();
        arr->SetArray(reinterpret_cast<uint16_t *>(voxel_data->data()),
                      voxel_data->size() / sizeof(uint16_t),
                      1);
        data_array = arr;
    } else if (voxel_type_string == "float32") {
        auto arr = vtkSmartPointer<vtkFloatArray>::New();
        arr->SetArray(reinterpret_cast<float *>(voxel_data->data()),
                      voxel_data->size() / sizeof(float),
                      1);
        data_array = arr;
    } else if (voxel_type_string == "float64") {
        auto arr = vtkSmartPointer<vtkDoubleArray>::New();
        arr->SetArray(reinterpret_cast<double *>(voxel_data->data()),
                      voxel_data->size() / sizeof(double),
                      1);
        data_array = arr;
    } else {
        throw std::runtime_error("Unrecognized voxel type " + voxel_type_string);
    }

    const vec3f grid_spacing = get_vec<float, 3>(config["spacing"]);
    vtkSmartPointer<vtkImageData> img_data = vtkSmartPointer<vtkImageData>::New();
    img_data->SetDimensions(brick_dims.x, brick_dims.y, brick_dims.z);
    img_data->SetSpacing(grid_spacing.x, grid_spacing.y, grid_spacing.z);
    // For whole volume iso we need to "untranslate" it so it's placed properly by the
    // instance transform we'll apply later, since I didn't want to add a separate instance
    if (isosurface_full_volume) {
        img_data->SetOrigin(-brick.ghost_bounds.lower.x,
                            -brick.ghost_bounds.lower.y,
                            -brick.ghost_bounds.lower.z);
    }
    img_data->GetPointData()->SetScalars(data_array);

    auto isovals = config["isovalue"].get<std::vector<float>>();
    for (const auto &val : isovals) {
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

std::shared_ptr<std::vector<uint8_t>> load_raw_volume(const std::string &file,
                                                      const std::string &dtype,
                                                      const vec3i &vol_dims,
                                                      const vec3i &brick_dims,
                                                      const vec3i &brick_offset)
{
    size_t voxel_size = 0;
    MPI_Datatype voxel_type;
    if (dtype == "uint8") {
        voxel_type = MPI_UNSIGNED_CHAR;
        voxel_size = 1;
    } else if (dtype == "uint16") {
        voxel_type = MPI_UNSIGNED_SHORT;
        voxel_size = 2;
    } else if (dtype == "float32") {
        voxel_type = MPI_FLOAT;
        voxel_size = 4;
    } else if (dtype == "float64") {
        voxel_type = MPI_DOUBLE;
        voxel_size = 8;
    } else {
        throw std::runtime_error("Unrecognized voxel type " + dtype);
    }

    const size_t n_voxels = brick_dims.long_product();
    auto voxel_data = std::make_shared<std::vector<uint8_t>>(n_voxels * voxel_size, 0);

    // MPI still uses 32-bit signed ints for counts of objects, so we have to split reads
    // of large data up so the count doesn't overflow. This assumes each X-Y slice is within
    // that size limit and reads chunks
    const size_t n_chunks = n_voxels / std::numeric_limits<int32_t>::max() +
                            (n_voxels % std::numeric_limits<int32_t>::max() > 0 ? 1 : 0);

    MPI_File file_handle;
    auto rc = MPI_File_open(
        MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_handle);
    if (rc != MPI_SUCCESS) {
        std::cerr << "[error]: Failed to open file " << file
                  << ". MPI Error: " << get_mpi_error(rc) << "\n";
        throw std::runtime_error("Failed to open " + file);
    }
    for (size_t i = 0; i < n_chunks; ++i) {
        const size_t chunk_thickness = brick_dims.z / n_chunks;
        const vec3i chunk_offset(
            brick_offset.x, brick_offset.y, brick_offset.z + i * chunk_thickness);
        vec3i chunk_dims = vec3i(brick_dims.x, brick_dims.y, chunk_thickness);
        if (i * chunk_thickness + chunk_thickness >= brick_dims.z) {
            chunk_dims.z = brick_dims.z - i * chunk_thickness;
        }
        const size_t byte_offset = i * chunk_thickness * brick_dims.y * brick_dims.x;
        const int chunk_voxels = chunk_dims.long_product();

        MPI_Datatype brick_type;
        MPI_Type_create_subarray(3,
                                 &vol_dims.x,
                                 &chunk_dims.x,
                                 &chunk_offset.x,
                                 MPI_ORDER_FORTRAN,
                                 voxel_type,
                                 &brick_type);
        MPI_Type_commit(&brick_type);

        MPI_File_set_view(file_handle, 0, voxel_type, brick_type, "native", MPI_INFO_NULL);
        rc = MPI_File_read_all(file_handle,
                               voxel_data->data() + byte_offset,
                               chunk_voxels,
                               voxel_type,
                               MPI_STATUS_IGNORE);
        if (rc != MPI_SUCCESS) {
            std::cerr << "[error]: Failed to read all voxels from file. MPI Error: "
                      << get_mpi_error(rc) << "\n";
            throw std::runtime_error("Failed to read all voxels from file");
        }
        MPI_Type_free(&brick_type);
    }
    MPI_File_close(&file_handle);
    return voxel_data;
}
