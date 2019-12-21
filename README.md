# A Cinema Mini-app Example Using OSPRay

This app provides an example of performing
[Cinema](https://github.com/Kitware/cinema)-style
in situ distributed rendering using OSPRay's [mpi_module](https://github.com/Kitware/cinema)
and the new asynchronous rendering API.

The libIS version of the app takes a JSON config file
describing the desired rendering configurations, but pulls data
from a simulation using [libIS]().
See the [libis example config](configs/libis_example_sim.json), which is commented below.
To run the app you specify the config file, the simulation host running libIS
to connect to (`-server` and `-port`), and the number of timesteps to query and render (`-n`).

```
mpirun -n <N> ./mini-cinema ../configs/libis_example_sim.json \
    -server localhost -port 29374 -n 10
```

```json
{
    "The field name to render from the list of fields sent back by the simulation"
    "field": "field_one",

    "The physical size of each voxel"
    "spacing": [1, 1, 1],

    "The dimensions of the image to render"
    "image_size": [512, 512],

    "Samples per-pixel to take when rendering each image"
    "spp": 2,

    "A list of camera positions or orbit counts to render"
    "camera": [
        {
            "A specific camera position can be specified by the
            camera position, direction and up vector"
            "pos": [128, 128, -256],
            "dir": [0, 0, 1],
            "up": [0, 1, 0]
        },
        {
            "A set of camera positions orbiting a sphere around
            the data set can be specified as an orbit, and the
            number of points on this orbit to render"
            "orbit": 2
        }
    ],

    "A list of colormaps to render with, can be RGB or full RGBA.
    RGBA colormaps will use the alpha channel to set opacity in
    the transfer function. Paths should be relative to the directory
    containing the JSON config file"
    "colormap": [
        "ice_fire.png",
        "jet.png",
        "paraview_cool_warm.png"
    ],

    "A set of isovalues to render with (requires VTK, since it tests
    using explicit geometry). As the value range may be unknown these
    are specified as locations between value_min and value_max to compute
    the isosurface on: isoval = lerp(value_min, value_max, t)"
    "isovalue": [0.2, 0.8],

    "Color to use for the isosufaces"
    "isosurface_color": [0.2, 0.5, 1.0],

    "Image background color, values in [0, 1]"
    "background_color": [0, 0, 0]
}
```

## Dependencies

- Intel TBB
- MPI
- [OSPRay 2.0](https://github.com/ospray/ospray/tree/release-2.0.x), with the [MPI module](https://github.com/ospray/module_mpi/). This will likely need to be built from source
- [libIS](https://github.com/ospray/libIS/)
- VTK (optional, for isosurface extraction)

