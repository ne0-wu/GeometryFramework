## Build Instructions

This project uses CMake to build. You only need to run the following commands to build the project:

```
mkdir build
cd build
cmake ..
````

Then you can build the project using your preferred build system. For example, if you are using Visual Studio, you can open the generated solution file and build the project from there.

This project uses CMake's FetchContent to download OpenMesh and GLFW. You don't need to download or install any dependencies manually. However, you will need an internet connection to download the dependencies.

This project has been tested on Windows 10 with Visual Studio 2022. It should work on other platforms and compilers, but it has not been tested. Note that the `CMakeLists.txt` automatically copies required obj files and shader files to the build directory, so you can build and run directly. This part has not been tested for other platforms or compilers.

## Controls

- `WASD` to move the camera.
- `Up` and `Down` arrow keys to zoom in and out.
- `R` to reset the camera.
-  `Control+S` to save the current frame as a PNG image.
 
## About homeworks

- **Homework 1** Press `Q` to simplify the mesh with QEM. The target is to reduce 10% of the vertices each time.
- **Homework 2** No user interface has been implemented. The program will automatically open 'meshes/spot_simplified.obj', sample a point cloud on it and reconstruct the mesh using Poisson reconstruction. The result will be saved as 'output.obj' in 'build/HelloOpenGL/{Release or Debug or ...}'.