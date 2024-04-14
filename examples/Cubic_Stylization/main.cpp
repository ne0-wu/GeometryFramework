#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"

#include "GeometryProcessing.h"

int main()
{
    Scene scene;

    // --------------------------------------------------

    // Load mesh and fit into unit ball
    Mesh mesh("meshes/camelhead.obj");

    scene.addMesh(mesh);

    // --------------------------------------------------

    // auto initialStateKeyQ = glfwGetKey(scene.window, GLFW_KEY_Q);
    auto initialStateKeyP = glfwGetKey(scene.window, GLFW_KEY_P);
    int numIter = -2;

    // render loop
    while (!scene.shouldClose())
    {
        // input
        // -----
        scene.processInput();

        // render
        // ------
        scene.draw();

        scene.swapBuffers();
        scene.pollEvents();
    }

    return 0;
}