#include "Mesh.h"
#include "GeometryProcessing.h"

int main()
{
    Mesh mesh("meshes/spot_quadrangulated.obj");

    QEMSimplification qem(mesh);

    std::cout << "Number of vertices: " << qem.getMesh().numVertices() << std::endl;

    qem.simplify(500);

    qem.getMesh().save("output.obj");

    return 0;
}