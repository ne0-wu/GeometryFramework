#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <filesystem>

#include <Mesh.h>
#include <Shader.h>

float vertices[] = {
	0.5f, 0.5f, 0.0f,	// 右上角
	0.5f, -0.5f, 0.0f,	// 右下角
	-0.5f, -0.5f, 0.0f, // 左下角
	-0.5f, 0.5f, 0.0f	// 左上角
};

unsigned int indices[] = {
	0, 1, 3, // 第一个三角形
	1, 2, 3	 // 第二个三角形
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
	GLWindow glWindow;

	// glad: load all OpenGL function pointers
	// ---------------------------------------
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	Shader shaderProgram("shaders/basic.vert", "shaders/basic.frag");

	// set up vertex data (and buffer(s)) and configure vertex attributes
	// ------------------------------------------------------------------
	GLuint VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);
	// bind the Vertex Array Object
	glBindVertexArray(VAO);
	// bind the Vertex Buffer Object, copy the vertices into the buffer
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	// bind the Element Buffer Object, copy the indices into the buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	// set the vertex attributes pointers
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
	// enable vertex attributes pointers
	glEnableVertexAttribArray(0);

	// ------------------------------

	glBindVertexArray(0); // unbind VAO

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

		// draw our first triangle
		shaderProgram.use();
		glBindVertexArray(VAO);								 // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0); // set the count to 6 since we're drawing 6 vertices now (2 triangles); not 3!
		glBindVertexArray(0);								 // no need to unbind it every time

		glWindow.swapBuffers(); // swap buffers
		glfwPollEvents();		// poll IO events
	}

	glfwTerminate();
	return 0;
}