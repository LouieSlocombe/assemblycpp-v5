
#include "assemblyCalculator.h"
#include <time.h>             // for clock, clock_t
#include <iostream>           // for operator<<, basic_ostream, cout, ifstream
#include <string>             // for char_traits, allocator, operator+, ope...
#include <vector>             // for vector
#include "globalPrimitives.h" // for moleculeName
#include "graphio.h"          // for graphio
#include "improvedBnB.h"      // for improvedBnB
#include "molGraph.h"         // for molGraph
#include "molfileParser.h"    // for molfileParser

using namespace std;

void assemblyCalculator(string &fileName)
{
    cout << "opening " << fileName << '\n';
    string in = fileName + ".mol";
    string out = fileName + "Out";
    ifstream molfile(in.c_str());
    if (molfile.is_open())
    {
        ofstream outputFile(out.c_str());
        moleculeName = fileName + "Pathway";
        molGraph mol_graph;
        vector<double> coords;
        molfileParser(molfile, mol_graph);
        clock_t t1 = clock();
        outputFile << fileName << " has assembly index: ";
        improvedBnB(mol_graph, outputFile);
        clock_t t2 = clock();
        outputFile << "time to completion: " << t2 - t1 << '\n';
    }
    else
    {
        ifstream graphFile(fileName.c_str());
        if (graphFile.is_open())
        {
            ofstream outputFile(out.c_str());
            molGraph mol_graph;
            graphio(graphFile, mol_graph);
            moleculeName = fileName + "Pathway";
            clock_t startTime = clock();
            outputFile << fileName << " has assembly index: ";
            improvedBnB(mol_graph, outputFile);
            outputFile << "time to completion: " << clock() - startTime << '\n';
        }
        else
            cout << "No file found\n";
    }
}