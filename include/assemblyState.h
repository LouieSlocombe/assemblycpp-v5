/**
 * @file assemblyState.h
 * @brief Structs relating to the assembly state and paths etc
 */
#pragma once

#include <cstddef>            // for size_t
#include <string_view>        // for hash
#include <unordered_set>      // for unordered_set
#include <vector>             // for vector, operator==, allocator
#include "globalPrimitives.h" // for vi, standardBitset

/**
 * @brief vector<int> hash, used to hash assembly states
 */
template <>
struct std::hash<vi>
{
    size_t operator()(const vi &v) const
    {
        std::size_t seed = v.size();
        for (auto &i : v)
        {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

/**
 * @brief Struct which is inserted into the hash table to store assembly states and recover the pathway
 */
struct assemblyPath
{
    /// vi represneting canonical indices of the fragments of the assembly state
    vi key;
    /// @brief number of duplicated bonds found so far
    int sumDupBonds;
    /// @brief needed to reconstruct the pathway
    unsigned short match, duplicate;
    /// @brief Assembly state from which this state is generated. needed to reconstruct the pathway
    assemblyPath *parent;
};

/**
 * @brief Wrapper for pathway hash table because C++ unordered_map does not guarantee pointer will remain unchanged
 *
 */
struct apWrapper
{
    assemblyPath *ap;

    bool operator==(const apWrapper &ap2) const
    {
        return ap->key == ap2.ap->key;
    }
};

/// Pointer for the minimum assembly path
extern assemblyPath *minAssemblyPath;

/// Hash table for assembly states for pathway algorithm
extern std::unordered_set<apWrapper> pathAssemblyMap;

/**
 * @brief vi hash called by apWrapper
 *
 */
template <>
struct std::hash<apWrapper>
{
    size_t operator()(const apWrapper &ap) const
    {
        std::size_t seed = ap.ap->key.size();
        for (auto &i : ap.ap->key)
        {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

/**
 * @brief Assembly state data structure. Records the current state of this assembly pathway
 */
struct assemblyState
{
    /// @brief each mask represents a separate fragment as a boolean edge list
    std::vector<standardBitset> masks;
    /// @brief number of duplicated bonds
    int sumDupBonds = 0;
    /// @brief index of the state
    int ix = 0;
    /// @brief path that was used to generate this state
    assemblyPath *apPtr = nullptr;

    /**
     * @brief Return the maximum fragment size, by counting the number of set bits in the first mask
     *
     * @return int (maximum size of a fragment)
     */
    int maxFragSizeF();

    /**
     * @brief Old branch and bound heuristic, used only during initial enumeration
     *
     * @return int (maximum duplicatable bonds value). Lower bound MA is total
     * bonds - 1 - this value.
     */
    int maxDupBonds();

    /**
     * @brief Calculates the maximum number of duplicatable bonds if given a maximum fragment size and a maximal fragment list
     *
     * @param sizeListMain The bitset counts of each fragment
     * @param maxFragSize The maximum allowed fragment size
     * @param sizeList The bitset counts of all duplicates in each fragment
     * @return int The maximum number of duplicatable bonds
     */
    int maxDupBonds(vi &sizeListMain, int maxFragSize, vi &sizeList);

    /**
     * @brief Basically same function as above but takes a vector of bitsets and uses the count to find sizeList
     *
     * @param targetMasks The vector of bitsets used in place of sizeList
     */
    int maxDupBonds(vi &sizeListMain, int maxFragSize, std::vector<standardBitset> &targetMasks);

    /**
     * @brief Like the function above but finds the maximum duplicate bonds for a vector of vector of bitsets
     *
     * @param fragSizeList The result vector
     * @param targetMasks The vector of vector of bitsets
     */
    void maxDupBonds(vi &fragSizeList, int maxFragSize, std::vector<std::vector<standardBitset>> &targetMasks);

    /**
     * @brief The simple branch and bound from v4
     *
     * @param maxFragSize The maximum allowed fragment size
     * @return int The maximum number of duplicatable bonds
     */
    int maxDupBonds(int maxFragSize);

    /**
     * @brief Old branch and bound. Only used during initial enumeration
     *
     * @return int The Lower bound
     */
    int lowBoundAI();

    /**
     * @brief Function for calculating lower bound on MA given an estimate and a maximum allowed fragment size
     *
     * @param maxFragSize The maximum allowed fragment size
     * @param estimate The estimate
     * @return int The lower bound
     */
    int lowBoundAI(int maxFragSize, int estimate);

    /**
     * @brief Upper bound MA given the sum of duplicatable bonds
     *
     * @return int The upper bound
     */
    int AI();

    /**
     * @brief Calculates the hash for the assembly state by returning a vi which is subjected to the vi hash
     *
     * @return vi The vector<int> to be hashed
     */
    vi assemblyHashCalculator();

    /**
     * @brief Prints the assembly state
     *
     */
    void print();
};
