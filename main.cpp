#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>
#include "json.hpp"
#include "loader.h"
#include "util.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace ospray;
using namespace ospcommon;
using namespace ospcommon::math;
using json = nlohmann::json;

int mpi_rank = 0;
int mpi_size = 0;
box3f world_bounds;

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
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

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
    if (args.size() < 2) {
        std::cout << "[error]: A config file to render is required\n";
        return;
    }

    json config;
    {
        std::ifstream cfg_file(args[1].c_str());
        cfg_file >> config;
    }
    if (mpi_rank == 0) {
        std::cout << "Rendering Config:\n" << config.dump(4) << "\n";
    }

    // Prefix the volume file name with the path to the config file
    {
        const std::string base_path = get_file_basepath(args[1]);
        if (base_path != args[1]) {
            config["volume"] = base_path + "/" + config["volume"].get<std::string>();
        }
    }

    VolumeBrick brick = load_volume_brick(config, mpi_rank, mpi_size);

    world_bounds = box3f(vec3f(0), get_vec<float, 3>(config["dimensions"]));
    const vec2i img_size = get_vec<int, 2>(config["image_size"]);

    cpp::Camera camera("perspective");
    camera.setParam("aspect", static_cast<float>(img_size.x) / img_size.y);

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

    cpp::FrameBuffer framebuffer(img_size, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
    auto camera_set = config["camera"].get<std::vector<json>>();
    for (size_t i = 0; i < camera_set.size(); ++i) {
        const vec3f cam_pos = get_vec<float, 3>(camera_set[i]["pos"]);
        const vec3f cam_dir = get_vec<float, 3>(camera_set[i]["dir"]);
        const vec3f cam_up = get_vec<float, 3>(camera_set[i]["up"]);
        camera.setParam("position", cam_pos);
        camera.setParam("direction", cam_dir);
        camera.setParam("up", cam_up);
        camera.commit();

        framebuffer.clear();
        // TODO: render async
        framebuffer.renderFrame(renderer, camera, world);

        if (mpi_rank == 0) {
            uint32_t *fb = (uint32_t *)framebuffer.map(OSP_FB_COLOR);
            const std::string fname = "mini-cinema_" + std::to_string(i) + ".jpg";
            stbi_write_jpg(fname.c_str(), img_size.x, img_size.y, 4, fb, 90);
        }
    }
}

