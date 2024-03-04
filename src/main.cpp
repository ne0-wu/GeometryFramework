#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <filesystem>
#include <vector>

#include <GLMesh.h>
#include <Shader.h>

std::vector<GLfloat> vertices = {
	0.5f, 0.5f, 0.0f,	// Top right
	0.5f, -0.5f, 0.0f,	// Bottom right
	-0.5f, -0.5f, 0.0f, // Bottom left
	-0.5f, 0.5f, 0.0f	// Top left
};

std::vector<GLuint> indices = {
	0, 1, 3, // First triangle
	1, 2, 3	 // Second triangle
};

class GLWindow
{
public:
	GLWindow()
	{
		// glfw: initialize and configure
		// ------------------------------
		glfwInit();
		// Open 3.3 Core
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

		// glfw window creation
		// --------------------
		window = glfwCreateWindow(800, 600, "HelloOpenGL", NULL, NULL);
		if (window == NULL)
		{
			std::cout << "Failed to create GLFW window" << std::endl;
			glfwTerminate();
		}
		glfwMakeContextCurrent(window);
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback); // Resize Window
	}

	bool shouldClose()
	{
		return glfwWindowShouldClose(window);
	}

	void processInput()
	{
		// press Esc to close the window
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);
	}

	void swapBuffers()
	{
		glfwSwapBuffers(window);
	}

private:
	GLFWwindow *window;

	static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
	{
		glViewport(0, 0, width, height);
	}
};

int main()
{
	std::filesystem::path currentPath = std::filesystem::current_path();
	std::cout << "Current working directory: " << currentPath << std::endl;

	GLWindow glWindow;

	// glad: load all OpenGL function pointers
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	Shader shaderProgram("shaders/basic.vert", "shaders/basic.frag");

	// GLMesh glMesh(vertices, indices);

	GLMesh glMesh("meshes/bunny.obj");

	// render loop
	while (!glWindow.shouldClose())
	{
		// input
		// -----
		glWindow.processInput();

		// render
		// ------
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		shaderProgram.use();

		glMesh.draw();

		glWindow.swapBuffers();
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}