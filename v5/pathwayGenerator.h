/**
 * @brief Struct representing a recovered molecular state in a reconstruction pathway
 */
struct reconstructedState
{
    vector<molGraph> list;
    int assemblyIndex;
};

/**
 * @brief Struct representing a pair of edge masks used during pathway reconstruction
 */
struct reconstructedEdgelist
{
    vector<standardBitset> list;
    int assemblyIndex;
};

vector<reconstructedState> minAssemblyPathway;

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing edges to output
 * @param ofs Output file stream
 */
void printMaskAsEdgeList(standardBitset mask, ofstream &ofs)
{
    int msb;
    for (msb = univEdgeList.size() - 1; msb >= 0; msb--) if (mask[msb]) break;
    ofs << "[";
    for (size_t i = 0; i < univEdgeList.size(); i++)
    {
        if (mask[i])
        {
            ofs << "[" << univEdgeList[i].a << "," << univEdgeList[i].b << "]";
            if (i != msb) ofs << ',';
        }
    }
    ofs << "]";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param ofs Output file stream
 */
void printMaskAsEdgeList(ofstream &ofs)
{
    ofs << "[";
    for (size_t i = 0; i < originalEdgeList.size(); i++)
    {
        ofs << "[" << originalEdgeList[i].a << "," << originalEdgeList[i].b << "]";
        if (i < originalEdgeList.size() - 1) ofs << ',';
    }
    ofs << "]";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param maskList List of two bitsets (matching pair)
 * @param ofs Output file stream
 */
void printMatching(vector<standardBitset> &maskList, ofstream &ofs)
{
    ofs << "{\"Left\":";
    printMaskAsEdgeList(maskList[0], ofs);
    ofs << ",\"Right\":";
    printMaskAsEdgeList(maskList[1], ofs);
    ofs << "}";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param ofs Output file stream
 */
void printOriginalGraph(ofstream &ofs)
{
    ofs << "\"Vertices\": [";
    for (size_t i = 0; i < originalMolecule.mg.size(); i++)
    {
            ofs << i;
            if (i < originalMolecule.mg.size() - 1) ofs << ',';
    }
    ofs << "],\n";
    ofs << "\"Edges\": ";
    printMaskAsEdgeList(ofs);
    ofs << ",\n";
    ofs << "\"VertexColours\": [";
    for (size_t i = 0; i < originalMolecule.mg.size(); i++)
    {
        ofs << "\"" << originalMolecule.mg[i].type << "\"";
        if (i < originalMolecule.mg.size() - 1) ofs << ',';
    }
    ofs << "],\n";
    ofs << "\"EdgeColours\": [";
    for (size_t i = 0; i < originalEdgeList.size(); i++)
    {
        int x = originalMolecule.btypeS(originalEdgeList[i].a, originalEdgeList[i].c);
        if (x < 4)
        {
            switch (x)
            {
                case 1: ofs << "\"single\"";
                    break;
                case 2: ofs << "\"double\"";
                    break;
                case 3: ofs << "\"triple\"";
                    break;
                case 4: ofs << "\"quadruple\"";
                    break;
                case 5: ofs << "\"quintuple\"";
                    break;
                case 0: ofs << "error";
                    break;
            }
        }
        else
        {
            ofs << "\"" + to_string(x) + "\"";
        }
        if (i < originalEdgeList.size() - 1) ofs << ',';
    }
    ofs << "]\n";
}


/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing edges to remove from original molecule
 * @param ofs Output file stream
 */
void printRemnantGraph(standardBitset mask, ofstream &ofs)
{
    standardBitset dual = allEdges ^ mask, remnantAtoms = 0;
    int msb, msb2;
    for (msb = univEdgeList.size() - 1; msb >= 0; msb--) if (dual[msb]) break;
    for (size_t i = 0; i < univEdgeList.size(); i++)
    {
        if (dual[i])
        {
            standardBitset b1, b2; b1.set(univEdgeList[i].a), b2.set(univEdgeList[i].b);
            remnantAtoms |= (b1 | b2);
        }
    }
    for (msb2 = targetMolecule.mg.size() - 1; msb2 >= 0; msb2--) if (remnantAtoms[msb2]) break;
    ofs << "\"Vertices\": [";
    for (size_t i = 0; i < targetMolecule.mg.size(); i++)
    {
        if (remnantAtoms[i])
        {
            ofs << i;
            if (i != msb2) ofs << ',';
        }
    }
    ofs << "],\n";
    ofs << "\"Edges\": ";
    printMaskAsEdgeList(dual, ofs);
    ofs << ",\n";
    ofs << "\"VertexColours\": [";
    for (size_t i = 0; i < targetMolecule.mg.size(); i++)
    {
        if (remnantAtoms[i])
        {
            ofs << "\"" << targetMolecule.mg[i].type << "\"";
            if (i != msb2) ofs << ',';
        }
    }
    ofs << "],\n";
    ofs << "\"EdgeColours\": [";
    for (size_t i = 0; i < univEdgeList.size(); i++)
    {
        if (dual[i])
        {
            int x = targetMolecule.btypeS(univEdgeList[i].a, univEdgeList[i].c);
            if (x < 4)
            {
                switch (x)
                {
                    case 1: ofs << "\"single\"";
                        break;
                    case 2: ofs << "\"double\"";
                        break;
                    case 3: ofs << "\"triple\"";
                        break;
                    case 4: ofs << "\"quadruple\"";
                        break;
                    case 5: ofs << "\"quintuple\"";
                        break;
                    case 0: ofs << "error";
                        break;
                }
            }
            else
            {
                ofs << "\"" + to_string(x) + "\"";
            }
            if (i != msb) ofs << ',';
        }
    }
    ofs << "]\n";
}

/**
 * @brief Pathway reconstruction function, which outputs the original graph, remnants and duplicates to a file
 *
 * Outputs the pathway to a file whose name is stored in globals::moleculeName, which is the molecule name appended with "Pathway", e.g. "aspirinPathway"
 *
 * @param removedEdges Edges that were removed during assembly
 */
void recoverPathway2(vector<edgeL> &removedEdges)
{
    minAssemblyPathway.clear();
    assemblyPath *curr, *prev;
    stack<assemblyPath*> sap; sap.push(minAssemblyPath);
    vector<assemblyPath*> minPath;
    curr = sap.top();
    while (curr != nullptr)
    {
        prev = curr->parent;
        sap.push(prev);
        curr = prev;
    }
    sap.pop();
    while (sap.size())
    {
        minPath.push_back(sap.top());
        sap.pop();
    }
    vector<vector<standardBitset> > maskList(graphHashMap.size());
    for (auto it = graphHashMap.begin(); it != graphHashMap.end(); ++it)
    {
        maskList[it->second.first].resize(it->second.second + 1);
    }
    for (auto it = bitsetHashTable.begin(); it != bitsetHashTable.end(); ++it)
    {
        maskList[it->second.first][it->second.second] = it->first;
    }
    standardBitset allTakenEdges = 0;
    vector<reconstructedEdgelist> v;
    for (size_t i = 1; i < minPath.size(); i++)
    {
        reconstructedEdgelist r;
        standardBitset mask = maskList[minPath[i]->key[0]][minPath[i]->match], 
        duplicate = maskList[minPath[i]->key[0]][minPath[i]->duplicate];
        allTakenEdges |= duplicate;
        r.list.push_back(mask);
        r.list.push_back(duplicate);
        v.push_back(r);
    }
    ofstream ofs(moleculeName.c_str());
    ofs << "{\n";
    ofs << "\"file_graph\":[\n";
    ofs << "{\n";
    printOriginalGraph(ofs);
    ofs << "}\n";
    ofs << "],\n";
    ofs << "\"remnant\":[\n";
    ofs << "{\n";
    printRemnantGraph(allTakenEdges, ofs);
    ofs << "}\n";
    ofs << "],\n";
    ofs << "\"duplicates\":[\n";
    for (size_t i = 0; i < v.size(); i++)
    {
        printMatching(v[i].list, ofs);
        if (i < v.size() - 1) ofs << ",\n";
    }
    ofs << "\n],\n";
    ofs << "\"removed_edges\":[";    
    for (size_t i = 0; i < removedEdges.size(); i++)
    {
        ofs << "[" << removedEdges[i].a << "," << removedEdges[i].b << "]";
        if (i < removedEdges.size() - 1) ofs << ',';
    }
    ofs << "]\n";
    ofs << "}\n";
}