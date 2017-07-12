#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <omp.h>

int main(){

#pragma omp parallel for schedule(dynamic)
    for(int i=1;i<20;i++){

        std::string commandLine;

        std::stringstream ss;

        ss << i;

        ss >> commandLine;

        commandLine = "./LinearOpticalSimulation " + commandLine;

        system( commandLine.c_str() );

    }

    return 0;

}
