/**
 * @file assemblyCalculator.h
 * @brief Main functions that calculate the assembly index of a molecule
 *
 * The functions here are the entry point for the assembly index calculation.
 */
#pragma once
#include <string>

/**
 * @brief Function that takes a molfile fileName.mol and finds its assembly index
 *
 * This function will open the file, parse it, and then calls improvedBnb
 *
 * @param fileName The molfile in question
 */
void assemblyCalculator(std::string &fileName);