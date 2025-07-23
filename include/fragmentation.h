/**
 * @file fragmentation.h
 * @brief code relating to splitting the assembly state into fragments
 */
#pragma once
struct assemblyState;
struct validMatchings;

/**
 * @brief Splits the assembly state into fragments after a given duplicate is removed using the disjoint-set data structure
 *
 * @param _target The assembly state to be fragmented
 * @param validMatchings The matching with the duplicate pair. matching.first is retained and matching.second is deleted
 * @param _result The resulting assembly state
 */
void fragmentAssemblyState(assemblyState &_target, validMatchings &matching,
                           assemblyState &_result);

/**
 * @brief Empty the hash table
 *
 */
void clearPathMap();