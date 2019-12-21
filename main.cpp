#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <thread>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>
#include <tbb/task_group.h>
#include <tbb/tbb.h>
#include "json.hpp"
#include "loader.h"
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

const std::string USAGE =
    "./mini_cinema <config.json> [options]\n"
    "Options:\n"
    "  -prefix <name>       Provide a prefix to prepend to the image file names.\n"
    "  -fif <N>             Restrict the number of frames being rendered in parallel.\n"
    "  -server <host>       Specify rank 0 of the simulation to query data from.\n"
    "  -port <port>         Specify the port rank 0 of the simulation is listening for "
    "clients on\n"
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

    bool finished() const;

    void save_image() const;
};

AsyncRender::AsyncRender(cpp::Renderer renderer,
                         cpp::Camera camera,
                         cpp::World world,
                         vec2i img,
                         std::string output_file)
    : img_size(img),
      fb(img_size, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM),
      future(fb.renderFrame(renderer, camera, world)),
      output_file(output_file)
{
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
void process_finished_renders(std::vector<AsyncRender> &renders, tbb::task_group &tasks);

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
        std::cerr << "[error]: A config file to render is required\n";
        std::cout << USAGE << "\n";
        return;
    }

    json config;
    int frames_in_flight = -1;
    std::string prefix;
    std::string server;
    int port = -1;
    int num_queries = 10;
    for (size_t i = 1; i < args.size(); ++i) {
        if (args[i] == "-prefix") {
            prefix = args[++i] + "-";
        } else if (args[i] == "-fif") {
            frames_in_flight = std::stof(args[++i]);
            if (frames_in_flight <= 0) {
                std::cerr << "[error]: Frames in flight must be >= 1\n";
                return;
            }
        } else if (args[i] == "-server") {
            server = args[++i];
        } else if (args[i] == "-port") {
            port = std::stoi(args[++i]);
        } else if (args[i] == "-n") {
            num_queries = std::stoi(args[++i]);
        } else if (args[i] == "-h") {
            std::cout << USAGE << "\n";
            return;
        } else {
            std::ifstream cfg_file(args[i].c_str());
            cfg_file >> config;

            // Prefix the colormap file names with the path to the config file
            const std::string base_path = get_file_basepath(args[i]);
            if (base_path != args[i]) {
                for (size_t j = 0; j < config["colormap"].size(); ++j) {
                    config["colormap"][j] =
                        base_path + "/" + config["colormap"][j].get<std::string>();
                }
            }
        }
    }
    if (server.empty() || port < 0) {
        std::cerr
            << "[error]: A simulation server (-server) and port (-port) must be specified\n";
        throw std::runtime_error("Missing server or port arguments");
    }

    if (mpi_rank == 0) {
        std::cout << "Rendering Config:\n" << config.dump(4) << "\n";
    }

    // Connect to the simulation
    is::client::connect(server, port, MPI_COMM_WORLD);

    const vec2i img_size = get_vec<int, 2>(config["image_size"]);

    cpp::Material material("scivis", "default");
    if (config.find("isosurface_color") != config.end()) {
        material.setParam("Kd", get_vec<float, 3>(config["isosurface_color"]));
    } else {
        material.setParam("Kd", vec3f(1.f));
    }
    material.commit();

    // create and setup an ambient light
    cpp::Light ambient_light("ambient");
    ambient_light.commit();

    // WILL TODO: mpi raycast progressive refinement/jittering is disabled?
    // though after accumulation i do see some artifacts in scivis too
    cpp::Renderer renderer("mpi_raycast");
    renderer.setParam("light", cpp::Data(ambient_light));
    if (config.find("background_color") != config.end()) {
        renderer.setParam("bgColor", get_vec<float, 3>(config["background_color"]));
    }
    if (config.find("spp") != config.end()) {
        renderer.setParam("spp", config["spp"].get<int>());
    }
    renderer.setParam("volumeSamplingRate", 1.f);
    renderer.commit();

    const std::string fmt_string =
        "%0" + std::to_string(static_cast<int>(std::log10(num_queries)) + 1) + "d";
    std::string fmt_out_buf(static_cast<int>(std::log10(num_queries)) + 1, '0');
    for (int t = 0; t < num_queries; ++t) {
        if (!is::client::sim_connected()) {
            break;
        }
        // Query the data from the simulation
        auto regions = is::client::query();
        std::cout << "Rank " << mpi_rank << " received " << regions.size() << "regions\n";

        auto start = high_resolution_clock::now();
        if (mpi_rank == 0) {
            std::cout << "Beginning rendering\n";
        }

        // TODO: support for multiple regions, we'll need to pass the regions param to
        // the world.
        std::vector<VolumeBrick> bricks;
        for (auto &region : regions) {
            bricks.push_back(load_volume_brick(config, region, mpi_rank, mpi_size));
        }
        world_bounds = box3f(
            vec3f(regions[0].world.min.x, regions[0].world.min.y, regions[0].world.min.z),
            vec3f(regions[0].world.max.x, regions[0].world.max.y, regions[0].world.max.z));

        if (config.find("value_range") == config.end()) {
            vec2f local_value_range = bricks[0].value_range;
            for (size_t i = 1; i < bricks.size(); ++i) {
                local_value_range.x = std::min(bricks[i].value_range.x, local_value_range.x);
                local_value_range.y = std::max(bricks[i].value_range.y, local_value_range.y);
            }

            vec2f global_value_range;
            MPI_Allreduce(&local_value_range.x,
                          &global_value_range.x,
                          1,
                          MPI_FLOAT,
                          MPI_MIN,
                          MPI_COMM_WORLD);
            MPI_Allreduce(&local_value_range.y,
                          &global_value_range.y,
                          1,
                          MPI_FLOAT,
                          MPI_MAX,
                          MPI_COMM_WORLD);

            if (mpi_rank == 0) {
                std::cout << "Computed value range: " << global_value_range << "\n";
            }
            config["value_range"] = {global_value_range.x, global_value_range.y};
        }

        const vec2f value_range = get_vec<float, 2>(config["value_range"]);
        const auto colormaps =
            load_colormaps(config["colormap"].get<std::vector<std::string>>(), value_range);
        const auto camera_set =
            load_cameras(config["camera"].get<std::vector<json>>(), world_bounds);
        // We don't use the implicit isosurfaces geometry because I want to test
        // on non-volume objects, e.g. explicit triangle surfaces
        std::vector<std::vector<Isosurface>> isosurfaces;
        size_t num_isovalues = 0;
        if (config.find("isovalue") != config.end()) {
            num_isovalues = config["isovalue"].size();
            for (const auto &b : bricks) {
                isosurfaces.push_back(extract_isosurfaces(config, b, mpi_rank, value_range));
            }
        }

        tbb::task_group tasks;
        std::vector<AsyncRender> active_renders;
        for (size_t k = 0; k < std::max(num_isovalues, size_t(1)); ++k) {
            for (size_t j = 0; j < colormaps.size(); ++j) {
                std::vector<cpp::Instance> instances;
                std::vector<box3f> region_bounds;
                for (size_t rid = 0; rid < bricks.size(); ++rid) {
                    const auto &b = bricks[rid];
                    cpp::VolumetricModel model(b.brick);
                    model.setParam("transferFunction", colormaps[j]);
                    model.commit();

                    cpp::Group group;
                    group.setParam("volume", cpp::Data(model));
                    if (!isosurfaces.empty() && !isosurfaces[rid].empty() &&
                        isosurfaces[rid][k].n_triangles) {
                        cpp::GeometricModel gm(isosurfaces[rid][k].geometry);
                        gm.setParam("material", material);
                        gm.commit();
                        group.setParam("geometry", cpp::Data(gm));
                    }
                    group.commit();

                    cpp::Instance instance(group);
                    // We apply a translation to the instance to place it correctly in
                    // the distributed world
                    auto transform = affine3f::translate(b.ghost_bounds.lower);
                    instance.setParam("xfm", transform);
                    instance.commit();
                    instances.push_back(instance);
                    region_bounds.push_back(b.bounds);
                }

                cpp::World world;
                world.setParam("instance", cpp::Data(instances));
                world.setParam("regions", cpp::Data(region_bounds));
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
                    std::sprintf(&fmt_out_buf[0], fmt_string.c_str(), static_cast<int>(t));
                    fname += "-tfn" + std::to_string(j) + "-cam" + std::to_string(i) + "-df" +
                             fmt_out_buf + ".jpg";

                    active_renders.emplace_back(renderer, camera, world, img_size, fname);
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
        auto end = high_resolution_clock::now();
        if (mpi_rank == 0) {
            std::cout << "All renders completed in "
                      << duration_cast<milliseconds>(end - start).count() << "ms\n";
        }
        // Wait for image writes to complete
        tasks.wait();
        end = high_resolution_clock::now();
        if (mpi_rank == 0) {
            std::cout << "Full run (renders + saving images) completed in "
                      << duration_cast<milliseconds>(end - start).count() << "ms\n";
        }
    }
    is::client::disconnect();
}

void process_finished_renders(std::vector<AsyncRender> &renders, tbb::task_group &tasks)
{
    auto done = std::stable_partition(
        renders.begin(), renders.end(), [](const AsyncRender &a) { return !a.finished(); });

    for (auto it = done; it != renders.end(); ++it) {
        AsyncRender r = *it;
        if (mpi_rank == 0) {
            tasks.run([r]() { r.save_image(); });
        }
    }

    renders.erase(done, renders.end());
}
