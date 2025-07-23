/**
 * @file graphio.h
 * @brief Loads molGraph object from input file
 */
#pragma once
#include <fstream> // for ifstream
struct molGraph;

/**
 * @brief Graph I/O
 *
 * @param ifs input file
 * @param mg output molGraph
 */
void graphio(std::ifstream &ifs, molGraph &mg);