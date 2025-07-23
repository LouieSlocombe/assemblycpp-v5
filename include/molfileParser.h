/**
 * @file molfileParser.h
 * @brief convert input molfile to molGraph object
 */
#pragma once
#include <fstream> // for ifstream
#include <string>  // for string
struct molGraph;

/**
 * @brief Takes a molfile and turns it into a molGraph and returns the name of the molecule
 *
 * @param molfile the .mol input
 * @param mg the molGraph output of the function
 * @return string name of molecule
 */
std::string molfileParser(std::ifstream &molfile, molGraph &mg);