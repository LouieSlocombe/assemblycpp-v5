#include "ioflag.h"
#include <stdlib.h>           // for atoll
#include <iosfwd>             // for std
#include <string>             // for basic_string, string, stoi, allocator
#include <unordered_map>      // for unordered_map
#include <vector>             // for vector
#include "globalPrimitives.h" // for ENUM_MAX, disjointCompensation, isPathway

using namespace std;

void owTime(string &_runtimeMax)
{
    runTimeMax = atoll(_runtimeMax.c_str());
}

void owEnumMax(string &_enum_max)
{
    ENUM_MAX = stoi(_enum_max);
}

void owPathway(string &_isPathway)
{
    isPathway = stoi(_isPathway);
}

void owRemoveHydrogens(string &_removeHydrogens)
{
    removeHydrogens = stoi(_removeHydrogens);
}

void owDisjointCompensate(string &_removeHydrogens)
{
    disjointCompensation = stoi(_removeHydrogens);
}

std::unordered_map<string, void (*)(string &)> fptrTable;

void fillFptrTable()
{
    void (*f)(string &);
    f = &owTime;
    fptrTable[string("runTime")] = f;
    f = &owEnumMax;
    fptrTable[string("enumMax")] = f;
    f = &owPathway;
    fptrTable[string("pathway")] = f;
    f = &owRemoveHydrogens;
    fptrTable[string("removeHydrogens")] = f;
    f = &owDisjointCompensate;
    fptrTable[string("compensateDisjoint")] = f;
}

void flagParser(int argc, char **argv)
{
    fillFptrTable();
    vector<string> args(argc);
    for (int i = 0; i < argc; i++)
    {
        args[i] = argv[i];
        if (args[i][0] == '-')
        {
            string flag = args[i].substr(1, args[i].find("=") - 1);
            string par = args[i].substr(args[i].find("=") + 1, args[i].length() - flag.length() + 2);
            if (fptrTable.count(flag))
            {
                void (*f)(string &) = fptrTable[flag];
                f(par);
            }
        }
    }
}