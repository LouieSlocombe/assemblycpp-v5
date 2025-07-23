/**
 * @file ioflag.h
 * @brief TODO: add brief description
 */
#pragma once
#include <string>        // for string, basic_string
#include <unordered_map> // for unordered_map

/**
 * @brief runTime flag
 *
 */
void owTime(std::string &_runtimeMax);

/**
 * @brief subgraph enumeration limit flag
 *
 */
void owEnumMax(std::string &_enum_max);

/**
 * @brief pathway flag
 *
 */
void owPathway(std::string &_isPathway);

/**
 * @brief hydrogen removal flag
 *
 */
void owRemoveHydrogens(std::string &_removeHydrogens);

/**
 * @brief disjoint graph JAI compensation flag
 *
 */
void owDisjointCompensate(std::string &_removeHydrogens);

/// For parsing flags
extern std::unordered_map<std::string, void (*)(std::string &)> fptrTable;

/**
 * @brief For parsing flags
 *
 */
void fillFptrTable();

/**
 * @brief For parsing flags
 *
 */
void flagParser(int argc, char **argv);