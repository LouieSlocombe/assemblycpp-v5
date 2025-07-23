
#include <string>
#include <iostream>

std::string helpstring = R"(The following are possible flags for assemblyCpp:

-runTime=x : sets runTime to x miliseconds in windows and x microseconds in linux, default is unlimited

-enumMax=x: sets the number of subgraphs the algorithm is allowed to enumerate to x, default is 50,000,000. Setting this too large will cause memory issues for big graphs

-pathwayFolder=x: sets the pathway folder for the string assembly algorithm to x

-pathway=x: if x is 1 (the default), a pathway file is generated, else if x is 0 no pathway will be generated

-removeHydrogens=x: removes explicit hydrogens if greater than 0 (default), else leaves in explicit hydrogens. Only applies to molfile inputs.)";

void help()
{
    std::cout << helpstring << std::endl;
}