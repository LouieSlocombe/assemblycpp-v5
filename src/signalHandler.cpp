#include <iostream>           // std::cout
#include <cstdlib>            // std::exit
#include "globalPrimitives.h" // minAIfound, interruptFlag, removedEdges
#include "pathwayGenerator.h" // recoverPathway2 (if declared separately)
#include "signalHandler.h"

using namespace std;

#ifdef _WIN32
BOOL CtrlHandler(DWORD fdwCtrlType)
{
    switch (fdwCtrlType)
    {
    case CTRL_C_EVENT:
        cout << "min AI found so far: " << minAIfound << '\n';
        interruptFlag = true;
        return TRUE;

    default:
        return FALSE;
    }
}
#else
void signalHandler(int signum)
{
    cout << "Interrupt signal received. Lowest AI found: " << minAIfound << '\n';
    recoverPathway2(removedEdges);
    exit(signum);
}
#endif