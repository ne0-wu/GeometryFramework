#include <iostream>
#include <filesystem>
#include <vector>

#include "Renderer/Scene.h"
#include "Mesh.h"
#include "Renderer/Window.h"

#include "GeometryProcessing.h"

int main()
{
	Mesh originalMesh("meshes/spot_simplified_2.obj");

	SpectralSimplification sp_simp(originalMesh, 100);

	sp_simp.simplify(originalMesh.n_vertices() * 4 / 5);
	sp_simp.getMesh().save("output.obj");

	return 0;
}