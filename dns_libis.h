#ifndef DNS_LIBIS_H_
#define DNS_LIBIS_H_

extern "C" {

/**
 * Configure libIS.
 * \param udata An 1D array containing 3D volume data organized in xyz order.
 * \param nx32 Numer of voxels in x direction.
 * \param ny Number of voxels in y direction.
 * \param nz32 Number of voxels in z direction.
 * \param Xyst An index value of the first voxel along the y-axis for an
 * x-pencil sub-volume. \param Xysz Number of volxels along the y-axis for an
 * x-pencil sub-volume. \param Xzst An index value of the first voxel along the
 * z-axis for an x-pencil sub-volume. \param Xzsz Number of volxels along the
 * z-axis for an x-pencil sub-volume. \param LMFx A user-defined volume
 * stretching factor in x direction. \param LMFy A user-defined volume
 * stretching factor in y direction. \param LMFz A user-defined volume
 * stretching factor in z direction.
 */
void libis_initialize(double* udata, int* nx32, int* ny, int* nz32, int* Xyst,
                      int* Xysz, int* Xzst, int* Xzsz, double* LMFx,
                      double* LMFy, double* LMFz);
/**
 * Send assigned volume data.
 * \param udata An 1D array containing 3D volume data organized in xyz order.
 * \param nx32 Numer of voxels in x direction.
 * \param ny Number of voxels in y direction.
 * \param nz32 Number of voxels in z direction.
 */
void libis_send(double* udata, int* nx32, int* ny, int* nz32);

/**
 * Free libIS simulation state disconnect form libIS.
 */
void libis_finalize();
}

#endif
