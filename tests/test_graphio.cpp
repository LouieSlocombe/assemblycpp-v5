#include <catch2/catch_all.hpp>
#include "graphio.h"
#include "molGraph.h"
#include "test_utils.h" // for CoutSilencer
#include <fstream>
#include <filesystem>

extern std::filesystem::path g_repo_root;

TEST_CASE("graphio correctly parses graphio_test", "[graphio]")
{

    std::filesystem::path path = g_repo_root / "tests/data/graphio_test.txt";
    std::ifstream molFile(path);
    REQUIRE(molFile.is_open());

    molGraph mg;
    {
        CoutSilencer silence; // suppress molecule information on stdout from graphio
        graphio(molFile, mg);
    }
    REQUIRE(mg.mg.size() == 12);
}
