#include <chrono>
#include <iostream>
#include <string>

class TickTock
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> stop;
    std::string message;

public:
    TickTock()
        : start(std::chrono::high_resolution_clock::now()),
          message("Operation") {}

    TickTock(std::string message)
        : start(std::chrono::high_resolution_clock::now()),
          message(message) {}

    void tick()
    {
        start = std::chrono::high_resolution_clock::now();
    }

    void tock()
    {
        stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        std::cout << message << " took " << duration.count() << " seconds." << std::endl;
    }
};