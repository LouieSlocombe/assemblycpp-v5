#include <catch2/catch_all.hpp>
#include "molfileParser.h"
#include "molGraph.h"
#include "test_utils.h" // for CoutSilencer
#include <fstream>
#include <algorithm>
#include <cctype>
#include <string>
#include <filesystem>

extern std::filesystem::path g_repo_root;

inline void trim(std::string &s)
{
    // Left trim
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
                                    { return !std::isspace(ch); }));
    // Right trim
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
                         { return !std::isspace(ch); })
                .base(),
            s.end());
}

TEST_CASE("molfileParser correctly parses tryptophan.mol", "[molfileParser]")
{
    std::filesystem::path path = g_repo_root / "tests/data/tryptophan.mol";
    std::ifstream molFile(path);
    REQUIRE(molFile.is_open());

    molGraph mg;
    std::string molfile_name;

    {
        CoutSilencer silence; // suppress molecule information on stdout from molfileParser
        molfile_name = molfileParser(molFile, mg);
    }

    trim(molfile_name);
    REQUIRE(molfile_name == "tryptophan");
    REQUIRE(mg.mg.size() == 15);
    REQUIRE(mg.totalBonds == 16);
    REQUIRE(mg.degree(0) == 2);
    REQUIRE(mg.atype(0) == "C");
    REQUIRE(mg.atype(14) == "O");
}
