#include "dns_libis.h"

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <string>

#include "mpi.h"

#include "libIS/is_sim.h"

static libISSimState* libis_state;

void libis_send(double* udata, int* nx32, int* ny, int* nz32) {
#if 0  // dump raw data
    static int tstep = 0;
    int rank;
    MPI_Comm sim_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(sim_comm, &rank);

    int xsize = *nx32;
    int ysize = *ny;
    int zsize = *nz32;

    std::string name = "dns-rank" + std::to_string(rank) + "-t" +
                       std::to_string(tstep) + ".raw";
    FILE* fp = fopen(name.c_str(), "wb");
    if (fp) {
      fwrite(udata, sizeof(double), xsize * ysize * zsize, fp);
      fclose(fp);
    }
    ++tstep;
#endif
  libISProcess(libis_state);
}

void libis_initialize(double* udata, int* nx32, int* ny, int* nz32, int* Xyst,
                      int* Xysz, int* Xzst, int* Xzsz, double* LMFx,
                      double* LMFy, double* LMFz) {
  // NOTE: volume size should be (x,y,z) = (8pi, 2, 3pi)
  double world_size_x = *LMFx * 2.0 * M_PI * 1000;
  double world_size_y = *LMFy * 2.0 * 1000;
  double world_size_z = *LMFz * 2.0 * M_PI * 1000;

  double spacing_x = world_size_x / (*nx32);
  // TODO spacinngs vary along y-axis, use constant factor for now
  double spacing_y = world_size_y / (*ny);
  double spacing_z = world_size_z / (*nz32);

  double world_begin_x = 0.0;
  double world_begin_y = 0.0;
  double world_begin_z = 0.0;

  double world_end_x = world_size_x;
  double world_end_y = world_size_y;
  double world_end_z = world_size_z;

  double local_begin_x = world_begin_x;
  double local_begin_y = spacing_y * (*Xyst);
  double local_begin_z = spacing_z * (*Xzst);

  double local_end_x = world_end_x;
  double local_end_y = local_begin_y + (spacing_y * (*Xysz));
  double local_end_z = local_begin_z + (spacing_z * (*Xzsz));

  libISInit(MPI_COMM_WORLD, 29374);

  libis_state = libISMakeSimState();

  libISVec3f world_min{world_begin_x, world_begin_y, world_begin_z};
  libISVec3f world_max{world_end_x, world_end_y, world_end_z};

  libISBox3f world_bounds = libISMakeBox3f();
  libISBoxExtend(&world_bounds, &world_min);
  libISBoxExtend(&world_bounds, &world_max);
  libISSetWorldBounds(libis_state, world_bounds);

  libISVec3f local_min{local_begin_x, local_begin_y, local_begin_z};
  libISVec3f local_max{local_end_x, local_end_y, local_end_z};

  libISBox3f local_bounds = libISMakeBox3f();
  libISBoxExtend(&local_bounds, &local_min);
  libISBoxExtend(&local_bounds, &local_max);
  libISSetLocalBounds(libis_state, local_bounds);
  // libISSetGhostBounds(libis_state, ghost_bounds);
  //
  const uint64_t dimensions[3] = {*nx32, *Xysz, *Xzsz};
  libISSetField(libis_state, "field_u", dimensions, DOUBLE, udata);

#if 0  // dump size info
    int rank;
    MPI_Comm sim_comm = MPI_COMM_WORLD;
    MPI_Comm_rank(sim_comm, &rank);

    std::string name_size = "dns-rank" + std::to_string(rank) + ".size";
    FILE* fp_size = fopen(name_size.c_str(), "wt");
    if (fp_size) {
      fprintf(fp_size, "%d %d %d %f %f %f %f %f %f", dimensions[0],
              dimensions[1], dimensions[2], local_min.x, local_min.y,
              local_min.z, local_max.x, local_max.y, local_max.z);
      fclose(fp_size);
    }
#endif
}

void libis_finalize() {
  libISFreeSimState(libis_state);
  libISFinalize();
}

