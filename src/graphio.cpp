#include "graphio.h"
#include <stddef.h>           // for size_t
#include <fstream>            // for basic_istream::operator>>, basic_istream
#include <iostream>           // for cout
#include <string>             // for getline, char_traits, operator>>, oper...
#include <utility>            // for pair
#include <vector>             // for vector
#include <sstream>            // std::istringstream
#include "globalPrimitives.h" // for moleculeName, pii
#include "molGraph.h"         // for molGraph

using namespace std;

void graphio(ifstream &ifs, molGraph &mg)
{
    string s;
    int graphSize, a, b;
    vector<pii> edgeList;
    getline(ifs, moleculeName);
    istringstream issm(moleculeName);
    issm >> moleculeName;
    cout << "Name of graph is: " << moleculeName << '\n';
    getline(ifs, s);
    istringstream iss1(s);
    iss1 >> graphSize;
    getline(ifs, s);
    istringstream iss(s);
    while (iss >> a >> b)
    {
        edgeList.push_back(pii(a, b));
    }
    getline(ifs, s);
    istringstream iss2(s);
    for (size_t i = 0; i < graphSize; i++)
    {
        iss2 >> s;
        mg.addAtom(s);
    }
    getline(ifs, s);
    istringstream iss3(s);
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        iss3 >> a;
        mg.addBond(edgeList[i].first - 1, edgeList[i].second - 1, a);
    }
    mg.printToCout();
}