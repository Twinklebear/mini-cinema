#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace ospray;
using namespace ospcommon;
using namespace ospcommon::math;

int mpi_rank = 0;
int mpi_size = 0;
box3f world_bounds;

struct VolumeBrick {
    // the volume data itself
    cpp::Volume brick;
    cpp::VolumetricModel model;
    cpp::Group group;
    cpp::Instance instance;
    // the bounds of the owned portion of data
    box3f bounds;
    // the full bounds of the owned portion + ghost voxels
    box3f ghost_bounds;
};

VolumeBrick load_volume_brick();
void render_images(const std::vector<std::string> &args);

int main(int argc, char **argv)
{
    {
        int thread_capability = 0;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_capability);
        if (thread_capability != MPI_THREAD_MULTIPLE) {
            std::cerr << "[error]: Thread multiple is needed for asynchronous "
                      << "rendering but is not available.\n";
            return 1;
        }
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::cout << "rank " << mpi_rank << "/" << mpi_size << "\n";

    {
        ospLoadModule("mpi");
        cpp::Device device("mpi_distributed");
        device.commit();
        device.setCurrent();

        // set an error callback to catch any OSPRay errors and exit the application
        ospDeviceSetErrorFunc(ospGetCurrentDevice(), [](OSPError error, const char *msg) {
            std::cerr << "[OSPRay error]: " << msg << std::endl << std::flush;
            std::exit(error);
        });

        render_images(std::vector<std::string>(argv, argv + argc));
    }

    ospShutdown();
    MPI_Finalize();
    return 0;
}

void render_images(const std::vector<std::string> &args)
{
    VolumeBrick brick = load_volume_brick();

    vec3f cam_pos(world_bounds.center().x, world_bounds.center().y, -32);
    vec3f cam_dir(0, 0, 1);
    vec3f cam_up(0, 1, 0);
    cpp::Camera camera("perspective");
    camera.setParam("aspect", 1.f);
    camera.setParam("position", cam_pos);
    camera.setParam("direction", cam_dir);
    camera.setParam("up", cam_up);
    camera.commit();

    cpp::World world;
    world.setParam("instance", cpp::Data(brick.instance));
    world.setParam("regions", cpp::Data(brick.bounds));
    world.commit();

    // create and setup an ambient light
    cpp::Light ambient_light("ambient");
    ambient_light.commit();

    cpp::Renderer renderer("mpi_raycast");
    renderer.setParam("light", cpp::Data(ambient_light));
    renderer.commit();

    cpp::FrameBuffer framebuffer(vec2i(1024, 1024), OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
    framebuffer.clear();
    framebuffer.renderFrame(renderer, camera, world);

    if (mpi_rank == 0) {
        uint32_t *fb = (uint32_t *)framebuffer.map(OSP_FB_COLOR);
        stbi_write_jpg("test.jpg", 1024, 1024, 4, fb, 90);
    }
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

// Compute an X x Y x Z grid to have 'num' grid cells,
// only gives a nice grid for numbers with even factors since
// we don't search for factors of the number, we just try dividing by two
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

VolumeBrick load_volume_brick()
{
    const vec3i grid = compute_grid(mpi_size);
    const vec3i brick_id(
        mpi_rank % grid.x, (mpi_rank / grid.x) % grid.y, mpi_rank / (grid.x * grid.y));
    // The bricks are 64^3 + 1 layer of ghost voxels on each axis
    const vec3i brick_volume_dims = vec3i(32);
    const vec3i brick_ghost_dims = vec3i(brick_volume_dims + 2);

    // The grid is over the [0, grid * brick_volume_dims] box
    world_bounds = box3f(vec3f(0.f), vec3f(grid * brick_volume_dims));
    const vec3f brick_lower = brick_id * brick_volume_dims;
    const vec3f brick_upper = brick_id * brick_volume_dims + brick_volume_dims;

    VolumeBrick brick;
    brick.bounds = box3f(brick_lower, brick_upper);
    // we just put ghost voxels on all sides here, but a real application
    // would change which faces of each brick have ghost voxels dependent
    // on the actual data
    brick.ghost_bounds = box3f(brick_lower - vec3f(1.f), brick_upper + vec3f(1.f));

    brick.brick = cpp::Volume("structured_regular");

    brick.brick.setParam("dimensions", brick_ghost_dims);

    // we use the grid origin to place this brick in the right position inside
    // the global volume
    brick.brick.setParam("gridOrigin", brick.ghost_bounds.lower);

    // generate the volume data to just be filled with this rank's id
    const size_t n_voxels = brick_ghost_dims.x * brick_ghost_dims.y * brick_ghost_dims.z;
    std::vector<uint8_t> volumeData(n_voxels, static_cast<uint8_t>(mpi_rank));
    brick.brick.setParam("voxelData", cpp::Data(volumeData));

    // Set the clipping box of the volume to clip off the ghost voxels
    brick.brick.setParam("volumeClippingBoxLower", brick.bounds.lower);
    brick.brick.setParam("volumeClippingBoxUpper", brick.bounds.upper);
    brick.brick.commit();

    brick.model = cpp::VolumetricModel(brick.brick);
    cpp::TransferFunction tfn("piecewise_linear");
    std::vector<vec3f> colors = {vec3f(0.f, 0.f, 1.f), vec3f(1.f, 0.f, 0.f)};
    std::vector<float> opacities = {0.05f, 1.f};

    tfn.setParam("color", cpp::Data(colors));
    tfn.setParam("opacity", cpp::Data(opacities));
    // color the bricks by their rank, we pad the range out a bit to keep
    // any brick from being completely transparent
    vec2f value_range = vec2f(0, mpi_size);
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

