#pragma once

#include <vector>

#include "Mesh.h"
#include "Shader.h"

class GLFWwindow;

class Scene
{
private:
    GLFWwindow *window;
    std::vector<Mesh> meshes;

public:
    Scene();
    ~Scene();

    bool shouldClose();
    void processInput();
    void pollEvents();

    void swapBuffers();

    void draw();
    void update();
};
