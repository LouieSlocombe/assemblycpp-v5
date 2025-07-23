
#include <climits>
#include "globalPrimitives.h"
int DISCHARGE_FREQUENCY = 30000000, ENUM_MAX = 50000000;
int minAIfound = -1;
int recursiveCount = 1;
int assemblyIx = 0;

std::unordered_map<std::string, int> atypeHash;
std::vector<double> coords;
std::string moleculeName;
standardBitset allEdges;
volatile bool interruptFlag = false;
clock_t startTime = 0;
unsigned long long runTimeMax = ULLONG_MAX;
unsigned int totalBonds = 0;
std::vector<edgeL> removedEdges;
std::vector<edgeL> originalEdgeList, univEdgeList;
std::unordered_map<standardBitset, pii> bitsetHashTable;
bool isPathway = 1, removeHydrogens = 1, disjointCompensation = 0;