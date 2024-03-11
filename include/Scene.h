#pragma once

#include <vector>

#include <Eigen/Core>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "Mesh.h"
#include "Shader.h"

class GLFWwindow;

enum controlType
{
	FPS,
	ORBIT
};

struct Camera
{
	Eigen::Vector3f position = Eigen::Vector3f(1.5, -1.5, 0);
	Eigen::Vector3f direction = Eigen::Vector3f(-1, 1, 0).normalized();
	Eigen::Vector3f up = Eigen::Vector3f(0, 0, 1); // opposite to gravity direction
	float fovY = 45.0f;
	float near = 0.1f;
	float far = 100.0f;
	float aspect = 1.0f;
	controlType type = ORBIT;
} defaultCamera;

// struct Camera
// {
// 	Eigen::Vector3f position = Eigen::Vector3f(0, 0, -1);
// 	Eigen::Vector3f direction = Eigen::Vector3f(0, 0, 1).normalized();
// 	Eigen::Vector3f up = Eigen::Vector3f(1, 0, 0); // opposite to gravity direction
// 	float fovY = 45.0f;
// 	float near = 0.1f;
// 	float far = 100.0f;
// 	float aspect = 2.0f;
// } defaultCamera;

class Scene
{
private:
	float lastFrameTime = 0.0f;
	Camera camera = defaultCamera;

	struct
	{
		int KEYBOARD_S_STATE;
	} initialState;

	static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
	{
		glViewport(0, 0, width, height);
	}

public:
	GLFWwindow *window;
	std::vector<Mesh> meshes;

	Scene()
	{
		// glfw: initialize and configure
		// ------------------------------
		glfwInit();

		// Open 3.3 Core
		// -------------
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

		// glfw window creation
		// --------------------
		window = glfwCreateWindow(1000, 1000, "HelloOpenGL", NULL, NULL);
		if (window == NULL)
		{
			std::cout << "Failed to create GLFW window" << std::endl;
			glfwTerminate();
		}
		glfwMakeContextCurrent(window);
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback); // Resize Window

		// glad: load all OpenGL function pointers
		// ---------------------------------------
		if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		{
			std::cout << "Failed to initialize GLAD" << std::endl;
			// return -1;
		}

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glEnable(GL_DEPTH_TEST);
	}

	~Scene()
	{
		glfwTerminate();
	}

	bool shouldClose()
	{
		return glfwWindowShouldClose(window);
	}

	void pollEvents()
	{
		glfwPollEvents();
	}

	void swapBuffers()
	{
		glfwSwapBuffers(window);
	}

	Eigen::Matrix4f modelMatrix()
	{
		Eigen::Matrix4f model;
		model.setIdentity();
		model.block<3, 1>(0, 3) = -camera.position;
		return model;
	}

	Eigen::Matrix4f viewMatrix()
	{
		Eigen::Matrix4f view;

		view.setIdentity();
		auto right = camera.up.cross(camera.direction).normalized();
		auto top = camera.direction.cross(right);
		view.block<1, 3>(0, 0) = right.transpose();
		view.block<1, 3>(1, 0) = top.transpose();
		view.block<1, 3>(2, 0) = -camera.direction.transpose();

		return view;
	}

	Eigen::Matrix4f projectionMatrix()
	{
		Eigen::Matrix4f projection = Eigen::Matrix4f::Zero();

		float tanHalfFovY = tan(camera.fovY / 2);
		float depth = camera.far - camera.near;

		projection << 1 / (camera.aspect * tanHalfFovY), 0, 0, 0,
			0, 1 / tanHalfFovY, 0, 0,
			0, 0, -(camera.far + camera.near) / depth, -2 * camera.far * camera.near / depth,
			0, 0, -1, 0;

		return projection;
	}

	void processInput()
	{
		// Control + S to save the current mesh
		if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
		{
			if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS && initialState.KEYBOARD_S_STATE == GLFW_RELEASE)
			{
				saveFrame("output_numVertices=" + std::to_string(meshes[0].numVertices()) + ".png");
				initialState.KEYBOARD_S_STATE = GLFW_PRESS;
			}
			// skip camera control
			return;
		}
		initialState.KEYBOARD_S_STATE = glfwGetKey(window, GLFW_KEY_S);

		// Camera control
		float currentFrameTime = glfwGetTime();
		float deltaTime = currentFrameTime - lastFrameTime;
		lastFrameTime = currentFrameTime;

		float cameraSpeed;
		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
			cameraSpeed = 10.0f * deltaTime;
		else
			cameraSpeed = 3.0f * deltaTime;

		float rotationSpeed = deltaTime;

		switch (camera.type)
		{
		case controlType::FPS:
			// WASD to move the camera forward, left, backward, right
			if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
				camera.position += cameraSpeed * camera.direction.normalized();
			if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
				camera.position -= cameraSpeed * camera.direction.normalized();
			if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
				camera.position += cameraSpeed * camera.direction.cross(camera.up).normalized();
			if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
				camera.position -= cameraSpeed * camera.direction.cross(camera.up).normalized();
			if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
				camera.position += cameraSpeed * camera.up.normalized();
			if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
				camera.position -= cameraSpeed * camera.up.normalized();

			// Arrow keys to rotate the camera
			if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
				camera.direction = (camera.direction * cos(rotationSpeed) + camera.up * sin(rotationSpeed)).normalized();
			if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
				camera.direction = (camera.direction * cos(rotationSpeed) - camera.up * sin(rotationSpeed)).normalized();
			if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
				camera.direction = (camera.direction * cos(rotationSpeed) + camera.direction.cross(camera.up) * sin(rotationSpeed)).normalized();
			if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
				camera.direction = (camera.direction * cos(rotationSpeed) - camera.direction.cross(camera.up) * sin(rotationSpeed)).normalized();

			break;
		case controlType::ORBIT:
			float radius = camera.position.norm();
			float theta = atan2(camera.position.y(), camera.position.x());
			float phi = acos(camera.position.z() / radius);

			// WASD to move the camera north (up), west (left), south (down), east (right)
			if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
				phi += cameraSpeed;
			if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
				phi -= cameraSpeed;
			if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
				theta -= cameraSpeed;
			if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
				theta += cameraSpeed;

			// Arrow keys to zoom in and out
			if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
				radius -= cameraSpeed;
			if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
				radius += cameraSpeed;

			phi = std::max(0.2f, std::min(phi, 3.14f - 0.2f));
			theta = fmod(theta, 2 * 3.14159f);
			camera.position.x() = radius * sin(phi) * cos(theta);
			camera.position.y() = radius * sin(phi) * sin(theta);
			camera.position.z() = radius * cos(phi);
			camera.direction = -camera.position.normalized();

			break;
		}

		// Reset camera to default position
		if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
			camera = defaultCamera;

		// Press Esc to close the window
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);
	}

	void draw()
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		for (auto &mesh : meshes)
		{
			mesh.bindShader();

			int modelLoc = glGetUniformLocation(mesh.getShaderID(), "model"),
				viewLoc = glGetUniformLocation(mesh.getShaderID(), "view"),
				projectionLoc = glGetUniformLocation(mesh.getShaderID(), "projection");

			glUniformMatrix4fv(modelLoc, 1, GL_FALSE, modelMatrix().data());
			glUniformMatrix4fv(viewLoc, 1, GL_FALSE, viewMatrix().data());
			glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, projectionMatrix().data());

			mesh.draw();
		}
	}

	void update();

	void addMesh(Mesh &mesh)
	{
		meshes.push_back(mesh);
		meshes.back().initializeGLMesh();
		meshes.back().initializeShader("shaders/basic.vert", "shaders/basic.frag");
	}

	void saveFrame(std::string filename)
	{
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);

		std::vector<unsigned char> pixels(3 * width * height);
		glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

		stbi_flip_vertically_on_write(true);
		stbi_write_png(filename.c_str(), width, height, 3, pixels.data(), 0);
	}
};