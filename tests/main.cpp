#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <filesystem>
#include <iostream>

std::filesystem::path g_repo_root;

std::filesystem::path find_repo_root(std::filesystem::path start)
{
    while (!start.empty())
    {
        if (std::filesystem::exists(start / ".git") || std::filesystem::exists(start / "CMakeLists.txt"))
        {
            return start;
        }
        start = start.parent_path();
    }
    throw std::runtime_error("Could not find project root from path: " + start.string());
}

int main(int argc, char *argv[])
{
    try
    {
        std::filesystem::path exe_path = argv[0];

        if (!exe_path.is_absolute())
        {
            exe_path = std::filesystem::current_path() / exe_path;
        }

        try
        {
            exe_path = std::filesystem::canonical(exe_path);
        }
        catch (...)
        {
            exe_path = std::filesystem::absolute(argv[0]); // fallback
        }

        g_repo_root = find_repo_root(exe_path.parent_path());

        if (argc == 1 || std::string(argv[1]) != "--list-tests")
        {
            std::cout << "[DEBUG] Detected g_repo_root = " << g_repo_root << std::endl;
        }

        return Catch::Session().run(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
}
