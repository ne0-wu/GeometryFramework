#include <random>

#include <Eigen/Core>

#include "Mesh.h"

struct PointCloudData
{
    Eigen::Vector3d point;
    Eigen::Vector3d normal;

    PointCloudData(Eigen::Vector3d point, Eigen::Vector3d normal)
        : point(point), normal(normal) {}
};

class PointCloud
{
    enum class SAMPLE_TYPE
    {
        SAMPLE_ON_VERTICES,
        SAMPLE_ON_FACES
    };

public:
    std::vector<PointCloudData> data;

    PointCloud(Mesh &mesh, SAMPLE_TYPE sampleType = SAMPLE_TYPE::SAMPLE_ON_VERTICES)
    {
        mesh.update_normals();

        switch (sampleType)
        {
        case SAMPLE_TYPE::SAMPLE_ON_VERTICES:
            data.reserve(mesh.numVertices());
            for (auto v : mesh.vertices())
                data.push_back(PointCloudData(mesh.point(v), mesh.normal(v)));
            break;
        case SAMPLE_TYPE::SAMPLE_ON_FACES:
            data.reserve(mesh.numFaces());
            for (auto f : mesh.faces())
            {
                for (auto v : mesh.fv_range(f))
                    data.push_back(PointCloudData(mesh.point(v), mesh.normal(f)));
            }
            break;
        }
    }

    PointCloud(Mesh &mesh, int numSamples)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0, 1);

        mesh.update_normals();
        data.reserve(numSamples);

        int eachFace = numSamples / mesh.numFaces();
        int remaining = numSamples % mesh.numFaces();

        std::vector<Mesh::FaceHandle> faceHandles;
        faceHandles.reserve(mesh.numFaces());

        auto sampleOnFace = [&](Mesh::FaceHandle f)
        {
            double u = dis(gen);
            double v = dis(gen);
            if (u + v > 1)
            {
                u = 1 - u;
                v = 1 - v;
            }
            double w = 1 - u - v;

            auto v_it = mesh.fv_iter(f);
            Eigen::Vector3d p0 = mesh.point(v_it);
            Eigen::Vector3d n0 = mesh.normal(v_it++);
            Eigen::Vector3d p1 = mesh.point(v_it);
            Eigen::Vector3d n1 = mesh.normal(v_it++);
            Eigen::Vector3d p2 = mesh.point(v_it);
            Eigen::Vector3d n2 = mesh.normal(v_it);

            Eigen::Vector3d point = u * p0 + v * p1 + w * p2;
            Eigen::Vector3d normal = u * n0 + v * n1 + w * n2;

            data.push_back(PointCloudData(point, normal));
        };

        for (auto f : mesh.faces())
        {
            faceHandles.push_back(f);
            for (int i = 0; i < eachFace; i++)
                sampleOnFace(f);
        }

        std::shuffle(faceHandles.begin(), faceHandles.end(), gen);
        for (auto f = faceHandles.begin(); f != faceHandles.end() && remaining > 0; f++, remaining--)
            sampleOnFace(*f);
    }

    void reserve(int numPoints)
    {
        data.reserve(numPoints);
    }
};