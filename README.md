## Build Instructions

This project uses CMake to build. You only need to run the following commands to build the project:

```
mkdir build
cd build
cmake ..
````

Then you can build the project using your preferred build system. For example, if you are using Visual Studio, you can open the generated solution file and build the project from there.

This project uses CMake's FetchContent to download dependencies. Thus internet connection is required to download the dependencies.

This project has been tested on Windows 10 with Visual Studio 2022. It should work on other platforms and compilers, but it has not been tested.

## Camera Controls

- `WASD` to move the camera.
- `Up` and `Down` arrow keys to zoom in and out.
- `R` to reset the camera.

## About Homeworks

The declaration of classes of homeworks is in the `include/GeometryProcessing.h`. The implementation is in the folder `src/GeometryProcessing/`. The visualization is in the folder `examples/`.

Modify `examples/CMakelists.txt` to choose which visualization to build.

Currently the build target is set to visualize the **cubic stylization** algorithm.