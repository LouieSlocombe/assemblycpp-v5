/**
 * @file globalPrimitives.h
 * @brief Various global variables and tyledefs
 */

#pragma once

#include <ctime>         // for clock_t
#include <bitset>        // for bitset
#include <string>        // for string
#include <unordered_map> // for unordered_map
#include <utility>       // for pair
#include <vector>        // for vector

typedef std::vector<int> vi;
typedef std::vector<bool> vb;
typedef std::pair<int, int> pii;

constexpr int BITSET_LENGTH = 512;
constexpr int MAX_INT = 2147483647;
constexpr int HASH_DEPTH_MAX = 7;
using standardBitset = std::bitset<BITSET_LENGTH>;

template <typename T1, typename T2, typename T3>
struct triple
{
    T1 a;
    T2 b;
    T3 c;
    triple() {}
    triple(T1 &_a, T2 &_b, T3 &_c) : a(_a), b(_b), c(_c) {}
};

extern int DISCHARGE_FREQUENCY, ENUM_MAX;
extern int minAIfound;
extern int recursiveCount;
extern int assemblyIx;

extern std::unordered_map<std::string, int> atypeHash;
extern std::vector<double> coords;
extern std::string moleculeName;
extern standardBitset allEdges;
extern volatile bool interruptFlag;
extern clock_t startTime;
extern unsigned long long runTimeMax;

typedef triple<int, int, int> iii;
typedef triple<short, short, short> edgeL;
extern unsigned int totalBonds;
extern std::vector<edgeL> removedEdges;
extern std::vector<edgeL> originalEdgeList, univEdgeList;

/// Hash table for edgelists for pathway algorithm
extern std::unordered_map<standardBitset, pii> bitsetHashTable;

extern bool isPathway, removeHydrogens, disjointCompensation;