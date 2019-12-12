#include <iostream>
#include <cstdio>
#include <mpi.h>
#include <ospray/ospray.h>
#include <ospray/ospray_cpp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace ospray;

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

    int mpi_rank = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::cout << "rank " << mpi_rank << "/" << mpi_size << "\n";

    MPI_Finalize();

    return 0;
}

