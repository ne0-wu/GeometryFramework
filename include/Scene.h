#pragma once

#include <vector>

#include <Eigen/Core>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "Mesh.h"
#include "Shader.h"

class GLFWwindow;

class Scene
{
private:
	GLFWwindow *window;
	std::vector<Mesh> meshes;
	struct Camera
	{
		Eigen::Vector3d position = Eigen::Vector3d(2, 2, 1);
		Eigen::Vector3d direction = Eigen::Vector3d(-1, -1, 0).normalized();
		Eigen::Vector3d up = Eigen::Vector3d(0, 0, 1); // opposite to gravity direction
	} camera;

	static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
	{
		glViewport(0, 0, width, height);
	}

public:
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

	Eigen::Matrix4d modelMatrix()
	{
		Eigen::Matrix4d model;
		model.setIdentity();
		model.block<3, 1>(0, 3) = -camera.position;
		return model;
	}

	Eigen::Matrix4d viewMatrix()
	{
		Eigen::Matrix4d view;
		view.setIdentity();
		auto right = camera.up.cross(camera.direction);
		auto top = camera.direction.cross(right);
		view.block<1, 3>(0, 0) = right.transpose();
		view.block<1, 3>(1, 0) = top.transpose();
		view.block<1, 3>(2, 0) = -camera.direction.transpose();
		return view;
	}

	Eigen::Matrix4d projectionMatrix() {}

	void processInput()
	{
		// press Esc to close the window
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);
	}

	void draw();
	void update();
};