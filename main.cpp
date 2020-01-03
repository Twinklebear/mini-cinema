#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <thread>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>
#include <tbb/task_group.h>
#include <tbb/tbb.h>
#include <unistd.h>
#include "json.hpp"
#include "loader.h"
#include "profiling.h"
#include "util.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

using namespace ospray;
using namespace ospcommon;
using namespace ospcommon::math;
using namespace std::chrono;
using json = nlohmann::json;

int mpi_rank = 0;
int mpi_size = 0;
box3f world_bounds;
bool save_images = true;
bool detailed_cpu_stats = false;

const std::string USAGE =
    "./mini_cinema <config.json> [options]\n"
    "Options:\n"
    "  -prefix <name>       Provide a prefix to prepend to the image file names.\n"
    "  -fif <N>             Restrict the number of frames being rendered in parallel.\n"
    "  -no-output           Don't save output images (useful for benchmarking).\n"
    "  -detailed-stats      Record and print statistics about CPU use, thread pinning, etc.\n"
    "  -h                   Print this help.";

struct AsyncRender {
    vec2i img_size;
    cpp::FrameBuffer fb = nullptr;
    cpp::Future future = nullptr;
    std::string output_file;

    AsyncRender() = default;
    AsyncRender(cpp::Renderer renderer,
                cpp::Camera camera,
                cpp::World world,
                vec2i img_size,
                std::string output_file);

    AsyncRender(const AsyncRender &) = delete;
    AsyncRender &operator=(const AsyncRender &) = delete;

    bool finished() const;

    void save_image() const;
};

AsyncRender::AsyncRender(cpp::Renderer renderer,
                         cpp::Camera camera,
                         cpp::World world,
                         vec2i img,
                         std::string output_file)
    : img_size(img), output_file(output_file)
{
    fb = cpp::FrameBuffer(img_size, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
    fb.setParam("timeCompositingOverhead", 0);
    fb.commit();
    future = fb.renderFrame(renderer, camera, world);
}

bool AsyncRender::finished() const
{
    return ospIsReady(future.handle(), OSP_TASK_FINISHED);
}

void AsyncRender::save_image() const
{
    uint32_t *img = (uint32_t *)fb.map(OSP_FB_COLOR);
    stbi_write_jpg(output_file.c_str(), img_size.x, img_size.y, 4, img, 90);
    fb.unmap(img);
}

void render_images(const std::vector<std::string> &args);
void process_finished_renders(std::vector<std::shared_ptr<AsyncRender>> &renders,
                              tbb::task_group &tasks);

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
        std::cerr << "[error]: A config file to render is required\n";
        std::cout << USAGE << "\n";
        return;
    }

    json config;
    int frames_in_flight = -1;
    std::string prefix;
    for (size_t i = 1; i < args.size(); ++i) {
        if (args[i] == "-prefix") {
            prefix = args[++i] + "-";
        } else if (args[i] == "-fif") {
            frames_in_flight = std::stof(args[++i]);
            if (frames_in_flight <= 0) {
                std::cerr << "[error]: Frames in flight must be >= 1\n";
                throw std::runtime_error("Frames in flight must be >= 1");
            }
        } else if (args[i] == "-no-output") {
            save_images = false;
        } else if (args[i] == "-detailed-stats") {
            detailed_cpu_stats = true;
        } else if (args[i] == "-h") {
            std::cout << USAGE << "\n";
            return;
        } else {
            std::ifstream cfg_file(args[i].c_str());
            if (!cfg_file) {
                std::cerr << "[error]: Failed to open config file " << args[i] << "\n";
                throw std::runtime_error("Failed to open input config file");
            }
            cfg_file >> config;
        }
    }

    if (mpi_rank == 0) {
        std::cout << "Rendering Config: " << config.dump() << "\n" << std::flush;
    }

    char hostname[HOST_NAME_MAX + 1] = {0};
    gethostname(hostname, HOST_NAME_MAX);
    if (detailed_cpu_stats) {
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < mpi_size; ++i) {
            if (i == mpi_rank) {
                std::cout << "rank " << mpi_rank << "/" << mpi_size << " on " << hostname
                          << "\n"
                          << get_file_content("/proc/self/status") << "\n=========\n"
                          << std::flush;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    VolumeBrick brick = load_volume_brick(config, mpi_rank, mpi_size);

    world_bounds = box3f(vec3f(0), get_vec<float, 3>(config["size"]));
    const vec2f value_range = get_vec<float, 2>(config["value_range"]);
    const vec2i img_size = get_vec<int, 2>(config["image_size"]);
    const auto colormaps =
        load_colormaps(config["colormap"].get<std::vector<std::string>>(), value_range);
    const auto camera_set =
        load_cameras(config["camera"].get<std::vector<json>>(), world_bounds);

    // create and setup an ambient light
    cpp::Light ambient_light("ambient");
    ambient_light.commit();

    cpp::Renderer renderer("mpi_raycast");
    renderer.setParam("light", cpp::Data(ambient_light));
    if (config.find("background_color") != config.end()) {
        renderer.setParam("bgColor", get_vec<float, 3>(config["background_color"]));
    }
    if (config.find("spp") != config.end()) {
        renderer.setParam("spp", config["spp"].get<int>());
    }
    bool isosurface_full_volume = false;
    if (config.find("ao") != config.end()) {
        // If we want AO we compute the isosurface of the entire volume on each rank
        // so that it's available for AO rays
        isosurface_full_volume = true;
        renderer.setParam("aoSamples", config["ao"].get<int>());
    }
    renderer.setParam("volumeSamplingRate", 1.f);
    renderer.commit();

    // We don't use the implicit isosurfaces geometry because I want to test
    // on non-volume objects, e.g. explicit triangle surfaces
    std::vector<cpp::Geometry> isosurfaces;
    if (config.find("isovalue") != config.end()) {
        isosurfaces = extract_isosurfaces(config, brick, mpi_rank, isosurface_full_volume);
    }

    cpp::Material material("scivis", "default");
    if (config.find("isosurface_color") != config.end()) {
        auto color = config["isosurface_color"].get<std::vector<float>>();
        material.setParam("Kd", vec3f(color[0], color[1], color[2]));
        if (color.size() == 4) {
            material.setParam("d", color[3]);
        }
    } else {
        material.setParam("Kd", vec3f(1.f));
    }
    material.commit();

    if (mpi_rank == 0) {
        std::cout << "Beginning rendering\n";
    }

    const std::string fmt_string =
        "%0" + std::to_string(static_cast<int>(std::log10(camera_set.size())) + 1) + "d";
    std::string fmt_out_buf(static_cast<int>(std::log10(camera_set.size())) + 1, '0');

    ProfilingPoint start;
    tbb::task_group tasks;
    std::vector<std::shared_ptr<AsyncRender>> active_renders;
    for (size_t k = 0; k < std::max(isosurfaces.size(), size_t(1)); ++k) {
        cpp::GeometricModel geom_model;
        if (!isosurfaces.empty()) {
            geom_model = cpp::GeometricModel(isosurfaces[k]);
            geom_model.setParam("material", material);
            geom_model.commit();
        }

        for (size_t j = 0; j < colormaps.size(); ++j) {
            cpp::VolumetricModel model(brick.brick);
            model.setParam("transferFunction", colormaps[j]);
            model.commit();

            cpp::Group group;
            group.setParam("volume", cpp::Data(model));
            if (geom_model) {
                group.setParam("geometry", cpp::Data(geom_model));
            }
            group.commit();

            cpp::Instance instance(group);
            // We apply a translation to the instance to place it correctly in
            // the distributed world
            auto transform = affine3f::translate(brick.ghost_bounds.lower);
            instance.setParam("xfm", transform);
            instance.commit();

            cpp::World world;
            world.setParam("instance", cpp::Data(instance));
            world.setParam("regions", cpp::Data(brick.bounds));
            world.commit();

            for (size_t i = 0; i < camera_set.size(); ++i) {
                cpp::Camera camera("perspective");
                camera.setParam("aspect", static_cast<float>(img_size.x) / img_size.y);
                camera.setParam("position", camera_set[i].pos);
                camera.setParam("direction", camera_set[i].dir);
                camera.setParam("up", camera_set[i].up);
                camera.commit();

                std::string fname = prefix + "mini-cinema";
                if (!isosurfaces.empty()) {
                    fname += "-iso" + std::to_string(k);
                }
                std::sprintf(&fmt_out_buf[0], fmt_string.c_str(), static_cast<int>(i));
                fname += "-tfn" + std::to_string(j) + "-cam" + fmt_out_buf + ".jpg";

                active_renders.push_back(
                    std::make_shared<AsyncRender>(renderer, camera, world, img_size, fname));
                process_finished_renders(active_renders, tasks);

                while (active_renders.size() >= frames_in_flight) {
                    process_finished_renders(active_renders, tasks);
                    std::this_thread::sleep_for(std::chrono::milliseconds(5));
                }
            }
        }
    }
    while (!active_renders.empty()) {
        process_finished_renders(active_renders, tasks);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    ProfilingPoint end;
    if (mpi_rank == 0) {
        std::cout << "All renders completed in " << elapsed_time_ms(start, end) << "ms\n";
    }
    if (detailed_cpu_stats) {
        for (int i = 0; i < mpi_size; ++i) {
            if (i == mpi_rank) {
                std::cout << "rank " << mpi_rank << "/" << mpi_size
                          << ", CPU: " << cpu_utilization(start, end) << "%\n";
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    // Wait for image writes to complete
    tasks.wait();
    end = ProfilingPoint();
    if (mpi_rank == 0) {
        std::cout << "Full run (renders + saving images) completed in "
                  << elapsed_time_ms(start, end) << "ms\n";
    }
}

void process_finished_renders(std::vector<std::shared_ptr<AsyncRender>> &renders,
                              tbb::task_group &tasks)
{
    auto done = std::stable_partition(
        renders.begin(), renders.end(), [](const std::shared_ptr<AsyncRender> &a) {
            return !a->finished();
        });

    if (save_images) {
        for (auto it = done; it != renders.end(); ++it) {
            auto r = *it;
            if (mpi_rank == 0) {
                tasks.run([r]() { r->save_image(); });
            }
        }
    }

    renders.erase(done, renders.end());
}
