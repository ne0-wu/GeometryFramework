#include <iostream>
#include <filesystem>
#include <vector>

#include "Mesh.h"
#include "GeometryProcessing.h"

int main()
{
    Mesh mesh("meshes/cathead.obj");

    CubicStylization cubic_style(mesh, 0.5, 100);
    cubic_style.stylize();
    cubic_style.get_stylized_mesh().save("result.obj");

    return 0;
}