/**
 * @file treeCanon.h
 * @brief tree canonisation code
 */
#pragma once
#include <string>             // for string
#include <vector>             // for vector
#include "globalPrimitives.h" // for pii
struct molGraph;

/**
 * @brief Subroutine of tree canonisation function
 *
 * @param mg Molecular graph (assumed to be a tree)
 * @param curr Current node index
 * @param parent Parent node index
 * @param type Bond type leading to this node
 * @return Canonical string for subtree rooted at curr
 */
std::string treeCanonRecursive(molGraph &mg, int curr, int parent, char type);

/**
 * @brief Subroutine of centroid finding function
 *
 * @param mg Molecular graph
 * @param weight Vector of subtree weights
 * @param curr Current node index
 * @param prev Previous node index
 * @return Subtree size rooted at curr
 */
int centroidDFS(molGraph &mg, std::vector<int> &weight, int curr, int prev);

/**
 * @brief Subroutine of tree canonisation function, finds centroid
 *
 * @param mg Molecular graph
 * @param currRoot Starting index to search for centroid
 * @return Pair of centroid(s) (second = -1 if only one)
 */
pii centroid(molGraph &mg, int currRoot);

/**
 * @brief Function to find a canonical string for a tree, using standard tree canonisation algorithm
 *
 * @param mg Acyclic molGraph
 * @param n Starting index
 * @return string (canonical form)
 */
std::string centroidTreeCanon(molGraph &mg, int n);