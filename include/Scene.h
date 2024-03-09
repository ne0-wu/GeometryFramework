#pragma once

#include <vector>

#include <Eigen/Core>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "Mesh.h"
#include "Shader.h"

class GLFWwindow;

// struct Camera
// {
// 	Eigen::Vector3f position = Eigen::Vector3f(2, 2, 0);
// 	Eigen::Vector3f direction = Eigen::Vector3f(-1, -1, 0).normalized();
// 	Eigen::Vector3f up = Eigen::Vector3f(0, 0, 1); // opposite to gravity direction
// 	float fov = 45.0f;
// 	float near = 0.1f;
// 	float far = 100.0f;
// 	float aspect = 2.0f;
// } defaultCamera;

struct Camera
{
	Eigen::Vector3f position = Eigen::Vector3f(0, 0, -1);
	Eigen::Vector3f direction = Eigen::Vector3f(0, 0, 1).normalized();
	Eigen::Vector3f up = Eigen::Vector3f(1, 0, 0); // opposite to gravity direction
	float fov = 45.0f;
	float near = 0.1f;
	float far = 100.0f;
	float aspect = 2.0f;
} defaultCamera;

class Scene
{
private:
	GLFWwindow *window;
	std::vector<Mesh> meshes;
	float lastFrameTime = 0.0f;
	Camera camera;

	static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
	{
		glViewport(0, 0, width, height);
	}

public:
	Shader shaderProgram;

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

		shaderProgram = Shader("shaders/basic.vert", "shaders/basic.frag");
		shaderProgram.use();
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
		auto right = camera.up.cross(camera.direction);
		auto top = camera.direction.cross(right);
		view.block<1, 3>(0, 0) = right.transpose();
		view.block<1, 3>(1, 0) = top.transpose();
		view.block<1, 3>(2, 0) = -camera.direction.transpose();
		return view;
	}

	Eigen::Matrix4f projectionMatrix()
	{
		Eigen::Matrix4f projection;
		projection.setIdentity();
		projection(0, 0) = 1 / tan(camera.fov / 2);
		projection(1, 1) = 1 / tan(camera.fov / 2);
		projection(2, 2) = -(camera.far + camera.near) / (camera.far - camera.near);
		projection(2, 3) = -2 * camera.far * camera.near / (camera.far - camera.near);
		projection(3, 2) = -1;
		projection(3, 3) = 0;
		return projection;
	}

	void processInput()
	{
		float currentFrameTime = glfwGetTime();
		float deltaTime = currentFrameTime - lastFrameTime;
		lastFrameTime = currentFrameTime;

		float cameraSpeed;
		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
			cameraSpeed = 10.0f * deltaTime;
		else
			cameraSpeed = 3.0f * deltaTime;

		float rotationSpeed = deltaTime;

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

		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
			camera.direction = (camera.direction * cos(rotationSpeed) + camera.up * sin(rotationSpeed)).normalized();
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
			camera.direction = (camera.direction * cos(rotationSpeed) - camera.up * sin(rotationSpeed)).normalized();
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
			camera.direction = (camera.direction * cos(rotationSpeed) - camera.direction.cross(camera.up) * sin(rotationSpeed)).normalized();
		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
			camera.direction = (camera.direction * cos(rotationSpeed) + camera.direction.cross(camera.up) * sin(rotationSpeed)).normalized();

		if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
			camera = defaultCamera;

		int modelLoc = glGetUniformLocation(shaderProgram.getID(), "model"),
			viewLoc = glGetUniformLocation(shaderProgram.getID(), "view"),
			projectionLoc = glGetUniformLocation(shaderProgram.getID(), "projection");

		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, modelMatrix().data());
		glUniformMatrix4fv(viewLoc, 1, GL_FALSE, viewMatrix().data());
		glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, projectionMatrix().data());

		// press Esc to close the window
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);
	}

	void draw();
	void update();
};