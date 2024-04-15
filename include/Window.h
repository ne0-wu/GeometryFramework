#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

class Window
{
protected:
	GLFWwindow *window;

public:
	Window(int width = 1280, int height = 720, const char *title = "Window")
	{
		setup_GLFW(width, height, title);
		setup_GLAD();
		setup_ImGui();
	}

	~Window()
	{
		glfwDestroyWindow(window);
		glfwTerminate();
	}

	bool shouldClose()
	{
		return glfwWindowShouldClose(window);
	}

	void swapBuffers()
	{
		glfwSwapBuffers(window);
	}

	void pollEvents()
	{
		glfwPollEvents();
	}

	void processInput()
	{
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, true);
		}
	}

private:
	static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
	{
		glViewport(0, 0, width, height);
	}

	void setup_GLFW(int width, int height, const char *title)
	{
		// Initialize GLFW
		if (!glfwInit())
		{
			throw std::runtime_error("Failed to initialize GLFW");
		}

		// Set OpenGL version to 3.3
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

		// Create a windowed mode window and its OpenGL context
		window = glfwCreateWindow(width, height, title, NULL, NULL);
		if (!window)
		{
			glfwTerminate();
			throw std::runtime_error("Failed to create window");
		}

		// Make the window's context current
		glfwMakeContextCurrent(window);

		// Set the callback function for window resize
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	}

	void setup_GLAD()
	{
		if (!gladLoadGL(glfwGetProcAddress))
		{
			throw std::runtime_error("Failed to initialize GLAD");
		}
	}

	void setup_ImGui()
	{
		// Setup Dear ImGui context
		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		ImGuiIO &io = ImGui::GetIO();
		(void)io;
		io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
		io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;  // Enable Gamepad Controls

		// Setup Dear ImGui style
		ImGui::StyleColorsDark();
		// ImGui::StyleColorsLight();

		// Setup Platform/Renderer bindings
		ImGui_ImplGlfw_InitForOpenGL(window, true);
		ImGui_ImplOpenGL3_Init("#version 330");
	}
};