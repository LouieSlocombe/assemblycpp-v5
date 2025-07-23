/**
 * @file duplicateMatching.h
 * @brief code relating to duplicates and potential duplciates
 */
#pragma once
#include <stddef.h>           // for size_t
#include <bitset>             // for bitset, operator&
#include <map>                // for map
#include <unordered_map>      // for unordered_map
#include <unordered_set>      // for unordered_set
#include <utility>            // for pair
#include <vector>             // for vector
#include "globalPrimitives.h" // for standardBitset, vb

/**
 * @brief Struct for storing a potential duplicate
 */
struct potentialDuplicate
{
    /// @brief mask representing edge list of potential duplicate
    standardBitset mask;
    /// @brief the index from the canonise function and index of the fragment
    int idx, fragment;
    potentialDuplicate() {}

    potentialDuplicate(standardBitset &_mask, int _fragment, int _idx) : mask(_mask), fragment(_fragment), idx(_idx) {}
};

/**
 * @brief Struct for storing a potential duplicate during the initial enumeration before the construction of the DAG
 *
 */
struct initialPotentialDuplicate : potentialDuplicate
{
    /// mask representing presence of specific atoms in the potential duplicate
    standardBitset atomMask = 0;
    /// mask representing the edge list of the parent fragment
    standardBitset fragMask = 0;
    /// Is the potential duplicate cyclic?
    bool isCyclic = 0;

    /**
     * @brief Construct a new potential Duplicate object
     *
     * @param x edge to be set
     * @param _fragMask Boolean edgelist of the fragment the duplicate is part of
     * @param _fragment Index of the fragment in its assembly state
     */
    initialPotentialDuplicate(int x, standardBitset &_fragMask, size_t _fragment);

    /**
     * @brief TODO: document
     */
    void generate(std::vector<initialPotentialDuplicate> &q, size_t fragment, std::unordered_set<standardBitset> &maskMap);

    /**
     * @brief Generate potential matches originating from this fragment and update the DAG
     *
     * @param q Potential duplicates which are isomorphic to this.mask
     * @param fragment Index of the fragment in its assembly state
     * @param maskMap Hash table of boolean edgelists of all fragments taken before to avoid repetition
     */
    void generateDAG(std::vector<initialPotentialDuplicate> &q, size_t fragment, std::unordered_set<standardBitset> &maskMap,
                     std::vector<std::unordered_map<standardBitset, std::pair<int, std::vector<standardBitset>>>> &tempDag);
};

/**
 * @brief Struct containing pair of valid duplicates for subsequent fragmentation
 *
 */
struct validMatchings
{
    /// @brief first and second duplicate masks
    standardBitset first, second;
    /// @brief frag 1 and frag 2 are indices of first, second respectively. maskFragSize is the maximum size of these fragments
    int frag1, frag2, maxFragSize;
    validMatchings() {}
    validMatchings(standardBitset &_first, standardBitset &_second, int _frag1, int _frag2, int _maxFragSize) : first(_first), second(_second), frag1(_frag1), frag2(_frag2), maxFragSize(_maxFragSize) {}
};

template <typename potentialDuplicate>
struct duplicateSet
{
    /// @brief bitset count of the duplicates in the duplicate set
    size_t size;
    /// @brief fragment masks from the assembly states
    std::vector<standardBitset> maskList;
    /// @brief list of potential duplicates
    std::vector<potentialDuplicate> list;
    duplicateSet() {}
    duplicateSet(size_t _size, size_t fragments)
    {
        size = _size;
        maskList.resize(fragments, 0);
    }

    /**
     * @brief Insert a potential duplicate into the list
     */
    void insert(const potentialDuplicate &m)
    {
        list.push_back(m);
        maskList[m.fragment] |= m.mask;
    }

    /**
     * @brief Check if there is a valid matching in the current maskList
     *
     * @return true if there is such a matching
     * @return false otherwise
     */
    bool isValid()
    {
        int count = 0, last = 0;
        for (size_t i = 0; i < maskList.size(); i++)
        {
            if (maskList[i] != 0)
            {
                count++;
                last = i;
            }
        }
        if (count > 1)
            return true;
        if (maskList[last].count() < (size << 1))
            return false;
        return true;
    }

    /**
     * @brief Generate matchings from pairable duplicates from the maskList
     *
     * @param v the output list of valid matchings
     * @return true if any valid matchings exist
     * @return false otherwise
     */
    bool generateMatchings(std::vector<validMatchings> &v)
    {
        bool output = 0;
        for (size_t i = 0; i < list.size(); i++)
        {
            size_t frag = list[i].fragment;
            if (size > 0)
            {
                for (size_t j = i + 1; j < list.size(); j++)
                {
                    if (frag == list[j].fragment)
                    {
                        if ((list[i].mask & list[j].mask) == 0)
                        {
                            validMatchings p(list[i].mask, list[j].mask, frag, frag, size);
                            v.push_back(p);
                        }
                    }
                    else
                    {
                        validMatchings p(list[i].mask, list[j].mask,
                                         list[i].fragment, list[j].fragment, size);
                        v.push_back(p);
                    }
                }
            }
        }
        return output;
    }
};

/**
 * @brief Set of boolean edgelists which are isomorphic
 *
 */
struct initialDuplicateSet : duplicateSet<initialPotentialDuplicate>
{
    using duplicateSet::duplicateSet;
    /**
     * @brief Generate size + 1 matchings from the current set and populate the DAG during the initial enumeration
     *
     * @param q list of potential duplicates
     * @param maskMap hash map to prevent bitsets from being added to the DAG more than once
     * @param tempDag the temporary DAG
     * @return true if any valid matchings exist
     * @return false otherwise
     */
    bool dagPopulator(std::vector<initialPotentialDuplicate> &q,
                      std::unordered_set<standardBitset> &maskMap,
                      std::vector<std::unordered_map<standardBitset, std::pair<int, std::vector<standardBitset>>>> &tempDag)
    {
        bool output = 0;
        vb alive(list.size(), 0);
        for (size_t i = 0; i < list.size(); i++)
        {
            size_t frag = list[i].fragment;
            if (size > 0)
            {
                for (size_t j = i + 1; j < list.size(); j++)
                {
                    if (frag == list[j].fragment)
                    {
                        if ((list[i].mask & list[j].mask) == 0)
                        {
                            alive[i] = 1;
                            alive[j] = 1;
                        }
                    }
                    else
                    {
                        alive[i] = 1;
                        alive[j] = 1;
                    }
                }
            }
            else
                alive[i] = 1;
            if (alive[i])
            {
                list[i].generateDAG(q, frag, maskMap, tempDag);
                output = 1;
            }
        }
        return output;
    }
};

/**
 * @brief Version of initialDuplicateSet which uses the DAG to search the duplicatable subgraph space more quickly
 */
struct dagDuplicateSet : duplicateSet<potentialDuplicate>
{
    bool dead = 1;
    using duplicateSet::duplicateSet;
};

/**
 * @brief Generate the next set of duplicates from the duplicate d
 *
 * @param d the duplicate from which the next set of duplicates is to be generated
 * @param stmap maps an integer corresponding to a unique index of each non-isomorphic graph to a set of potential duplicates
 * @param fragment the bitset corresponding to the fragment d is a part of
 * @param size the maximum allowed size of a duplicate
 * @param ordinal the maximum allowed index of a duplicate
 * @param frags the number of fragments in the assembly state
 * @return true if the canonical index of any duplicate is greater than the ordinal
 * @return false otherwise
 */
bool dagGenerate(potentialDuplicate &d, std::map<int, dagDuplicateSet> &stmap, standardBitset &fragment,
                 size_t size, int ordinal, size_t frags);

/**
 * @brief Generate the next set of duplicates from the duplicate set ds using the function dagGenerate
 *
 * @param ds the duplicate set from which the next set is to be generated
 * @param stmap maps an integer corresponding to a unique index of each non-isomorphic graph to a set of potential duplicates
 * @param takenMasks bitsets of all edges which could be part of a duplicate
 * @param stateMasks the bitsets of the original assembly state
 * @param ordinal the maximum allowed index of a duplicate
 * @param overweight true if the generation function has reached states which are larger than the ordinal
 * @param last true if this is to be the final iteration (previous iteration was overweight)
 * @return true if any valid duplicatable subgraphs found and not the final iteration
 * @return false
 */
bool dagDuplicateGenerator(dagDuplicateSet &ds, std::map<int, dagDuplicateSet> &stmap,
                           std::vector<standardBitset> &takenMasks, std::vector<standardBitset> &stateMasks, int ordinal, bool &overweight, bool last);