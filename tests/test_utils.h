#pragma once

#include <iostream>
#include <fstream>
#include <streambuf>

#ifdef _WIN32
#define DEV_NULL "NUL"
#else
#define DEV_NULL "/dev/null"
#endif

// RAII-based std::cout silencer
struct CoutSilencer
{
    std::streambuf *old_buf;
    std::ofstream null_stream;

    CoutSilencer() : null_stream(DEV_NULL)
    {
        old_buf = std::cout.rdbuf();
        std::cout.rdbuf(null_stream.rdbuf());
    }

    ~CoutSilencer()
    {
        std::cout.rdbuf(old_buf);
    }
};
