/**
 * @file pathwayGenerator.h
 * @brief code relating to reconstrucing assembly pathways produced by the program
 */
#pragma once
#include <fstream>            // for ofstream
#include <vector>             // for vector
#include "globalPrimitives.h" // for standardBitset, edgeL
#include "molGraph.h"         // for molGraph

/**
 * @brief Struct representing a recovered molecular state in a reconstruction pathway
 */
struct reconstructedState
{
    std::vector<molGraph> list;
    int assemblyIndex;
};

/**
 * @brief Struct representing a pair of edge masks used during pathway reconstruction
 */
struct reconstructedEdgelist
{
    std::vector<standardBitset> list;
    int assemblyIndex;
};

/**
 * @brief TODO: Add description
 */
extern std::vector<reconstructedState> minAssemblyPathway;

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing edges to output
 * @param ofs Output file stream
 */
void printMaskAsEdgeList(standardBitset mask, std::ofstream &ofs);

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param ofs Output file stream
 */
void printMaskAsEdgeList(std::ofstream &ofs);

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param maskList List of two bitsets (matching pair)
 * @param ofs Output file stream
 */
void printMatching(std::vector<standardBitset> &maskList, std::ofstream &ofs);

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param ofs Output file stream
 */
void printOriginalGraph(std::ofstream &ofs);

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing edges to remove from original molecule
 * @param ofs Output file stream
 */
void printRemnantGraph(standardBitset mask, std::ofstream &ofs);

/**
 * @brief Pathway reconstruction function, which outputs the original graph, remnants and duplicates to a file
 *
 * Outputs the pathway to a file whose name is stored in globals::moleculeName, which is the molecule name appended with "Pathway", e.g. "aspirinPathway"
 *
 * @param removedEdges Edges that were removed during assembly
 */
void recoverPathway2(std::vector<edgeL> &removedEdges);