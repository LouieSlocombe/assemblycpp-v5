#include <catch2/catch_all.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>
#include <optional>
#include <cstdlib>

namespace fs = std::filesystem;

struct IntegrationCase
{
    std::string base_name;
    int expected_index;
};

std::vector<IntegrationCase> read_csv(const std::filesystem::path &path)
{
    std::ifstream f(path);
    std::vector<IntegrationCase> cases;
    std::string line;

    while (std::getline(f, line))
    {
        if (line.empty())
            continue;

        std::istringstream iss(line);
        std::string base;
        int expected;

        if (iss >> base >> expected)
        {
            cases.push_back({base, expected});
        }
        else
        {
            std::cerr << "[SKIP] Malformed line: " << line << "\n";
        }
    }

    return cases;
}

struct AssemblyResult
{
    int index;
    double time_ms;
};

std::optional<AssemblyResult> extract_result(const std::string &out_file)
{
    std::ifstream f(out_file);
    std::string line;

    std::optional<int> index;
    std::optional<double> time;

    std::regex index_pattern(R"(has assembly index:\s*(\d+))");
    std::regex time_pattern(R"(time to completion:\s*(\d+))");
    std::smatch match;

    while (std::getline(f, line))
    {
        if (!index && std::regex_search(line, match, index_pattern))
        {
            index = std::stoi(match[1]);
        }
        else if (!time && std::regex_search(line, match, time_pattern))
        {
            time = std::stod(match[1]); // milliseconds
        }

        if (index && time)
            break;
    }

    if (index && time)
    {
        return AssemblyResult{*index, *time};
    }
    return std::nullopt;
}
void cleanup_outputs(const fs::path &directory)
{
    for (const auto &entry : fs::directory_iterator(directory))
    {
        std::string name = entry.path().filename().string();
        if ((name.size() >= 3 && name.compare(name.size() - 3, 3, "Out") == 0) ||
            (name.size() >= 7 && name.compare(name.size() - 7, 7, "Pathway") == 0))
        {
            fs::remove(entry.path());
        }
    }
}

TEST_CASE("Integration tests for ./build/bin/assembly")
{
    extern std::filesystem::path g_repo_root;

    fs::path binary;

#ifdef _WIN32
    binary = g_repo_root / "build/bin/assembly.exe";
    if (!fs::exists(binary))
    {
        binary = g_repo_root / "build/bin/Release/assembly.exe";
        if (!fs::exists(binary))
        {
            std::ostringstream oss;
            oss << "Could not find 'assembly.exe'. Checked:\n"
                << " - build/bin/assembly\n"
                << " - " << binary << "\n"
                << "Did you forget to build in Release mode for MSVC?";
            FAIL(oss.str());
        }
    }
#else
    binary = g_repo_root / "build/bin/assembly";
    if (!fs::exists(binary))
    {
        FAIL("Could not find 'assembly' in build/bin.\n"
             "Did you forget to build it?");
    }
#endif

    const auto csv_file = g_repo_root / "tests/integration/integration_test_data.csv";
    const auto mol_dir = g_repo_root / "tests/integration/molfiles";

    const auto failure_log_path = g_repo_root / "integration_failures.txt";
    const auto timing_csv_path = g_repo_root / "integration_timings.csv";

    std::ofstream log_file(failure_log_path);
    std::ofstream timing_csv(timing_csv_path);

    if (!log_file)
    {
        std::cerr << "[ERROR] Could not open failure log at: " << failure_log_path << "\n";
    }
    if (!timing_csv)
    {
        std::cerr << "[ERROR] Could not open timing log at: " << timing_csv_path << "\n";
    }

    timing_csv << "molecule,assembly_index,expected_index,time_seconds,status\n";

    auto test_cases = read_csv(csv_file);
    std::cout << "Loaded test cases: " << test_cases.size() << "\n";

    int passed = 0;
    int failed = 0;
    int total = test_cases.size();
    int index = 0;

    for (const auto &[base, expected_index] : test_cases)
    {
        index++;
        cleanup_outputs(mol_dir);

        std::string mol_path = (mol_dir / base).string();
        std::string out_file = (mol_dir / (base + "Out")).string();

        // clang-format off
        #ifdef _WIN32
                const std::string null_redirect = " > NUL";
        #else
                const std::string null_redirect = " > /dev/null";
        #endif
        // clang-format on

        std::string cmd = binary.string() + " " + mol_path + null_redirect;

        std::cout << "[ " << index << " / " << total << " ] " << base << "\r" << std::flush;

        auto t0 = std::chrono::steady_clock::now();
        int ret_code = std::system(cmd.c_str());
        auto t1 = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = t1 - t0;

        std::string status;
        std::string index_str;

        if (ret_code != 0)
        {
            status = "ERROR";
            std::string msg = "[SEGFAULT/ERROR] " + base + " (exit code " + std::to_string(ret_code) + ")";
            log_file << msg << "\n";
            INFO(msg);
            CHECK(false);
            failed++;
        }
        else if (!fs::exists(out_file))
        {
            status = "MISSING";
            std::string msg = "[MISSING OUTPUT] " + base + " (expected " + out_file + ")";
            log_file << msg << "\n";
            INFO(msg);
            CHECK(false);
            failed++;
        }
        else
        {
            auto maybe_result = extract_result(out_file);
            if (!maybe_result.has_value())
            {
                status = "PARSE_ERROR";
                std::string msg = "[PARSE ERROR] " + base + " (could not extract index or time)";
                log_file << msg << "\n";
                INFO(msg);
                CHECK(false);
                failed++;
            }
            else
            {
                index_str = std::to_string(maybe_result->index);
                if (maybe_result->index != expected_index)
                {
                    status = "FAIL";
                    std::ostringstream msg;
                    msg << "[FAIL] " << base << " → expected: " << expected_index
                        << ", got: " << maybe_result->index;
                    log_file << msg.str() << "\n";
                    INFO(msg.str());
                    CHECK(false);
                    failed++;
                }
                else
                {
                    status = "PASS";
                    passed++;
                }
                fs::remove(out_file);

                timing_csv << base << "," << index_str << "," << expected_index << ","
                           << maybe_result->time_ms / 1000.0 << "," << status << "\n"; // convert to seconds
            }
        }
    }
    cleanup_outputs(mol_dir);

    log_file.close();
    timing_csv.close();
    if (failed == 0)
    {
        std::cout << "\nAll integration tests passed on " << passed << " molecules.\n";
        fs::remove(failure_log_path); // keep timings, but remove failure log if clean
    }
    else
    {
        std::cout << "\n"
                  << failed << " molecules failed. See " << failure_log_path << " for details.\n";
    }

    REQUIRE(failed == 0);
}