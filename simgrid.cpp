#include "simgrid.h"
#include <algorithm>
#include <mpi.h>

SimBrick::SimBrick(const box3f &b, int owner) : bounds(b), id(-1), owner(owner) {}

std::ostream &operator<<(std::ostream &os, const SimGrid &sg)
{
    os << "SimGrid: " << sg.grid << ", bricks:\n";
    for (const auto &b : sg.bricks) {
        os << "\t(" << b.bounds << ", id: " << b.id << ", owner: " << b.owner << ")\n";
    }
    return os;
}

SimGrid reconstruct_grid(const std::vector<is::SimState> &regions,
                         const int mpi_rank,
                         const int mpi_size)
{
    // Build our list of local bricks
    std::vector<SimBrick> local_bricks;
    for (const auto &r : regions) {
        box3f b(vec3f(r.local.min.x, r.local.min.y, r.local.min.z),
                vec3f(r.local.max.x, r.local.max.y, r.local.max.z));
        local_bricks.emplace_back(b, mpi_rank);
    }

    SimGrid sim_grid;
    for (int i = 0; i < mpi_size; ++i) {
        if (i == mpi_rank) {
            int nbricks = local_bricks.size();
            MPI_Bcast(&nbricks, 1, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(
                local_bricks.data(), nbricks * sizeof(SimBrick), MPI_BYTE, i, MPI_COMM_WORLD);
            std::copy(
                local_bricks.begin(), local_bricks.end(), std::back_inserter(sim_grid.bricks));
        } else {
            int nbricks = 0;
            MPI_Bcast(&nbricks, 1, MPI_INT, i, MPI_COMM_WORLD);
            const size_t offset = sim_grid.bricks.size();
            sim_grid.bricks.resize(offset + nbricks);
            MPI_Bcast(&sim_grid.bricks[offset],
                      nbricks * sizeof(SimBrick),
                      MPI_BYTE,
                      i,
                      MPI_COMM_WORLD);
        }
    }

    // Sort into row-major ordering of the bricks. This assumes that the
    // bricks form some uniform grid
    std::sort(sim_grid.bricks.begin(),
              sim_grid.bricks.end(),
              [](const SimBrick &a, const SimBrick &b) {
                  return a.bounds.lower.z < b.bounds.lower.z ||
                         (a.bounds.lower.z == b.bounds.lower.z &&
                          a.bounds.lower.y < b.bounds.lower.y) ||
                         (a.bounds.lower.z == b.bounds.lower.z &&
                          a.bounds.lower.y == b.bounds.lower.y &&
                          a.bounds.lower.x < b.bounds.lower.x);
              });

    // Now determine the simulation processor grid dimensions
    vec3i brick_id(0);
    for (size_t i = 0; i < sim_grid.bricks.size(); ++i) {
        sim_grid.bricks[i].id = brick_id;

        if (i + 1 < sim_grid.bricks.size()) {
            if (sim_grid.bricks[i + 1].bounds.lower.z > sim_grid.bricks[i].bounds.lower.z) {
                brick_id.x = 0;
                brick_id.y = 0;
                ++brick_id.z;
            } else if (sim_grid.bricks[i + 1].bounds.lower.y > sim_grid.bricks[i].bounds.lower.y) {
                brick_id.x = 0;
                ++brick_id.y;
            } else {
                ++brick_id.x;
            }
        }
    }
    sim_grid.grid = brick_id + vec3i(1);

    return sim_grid;
}

