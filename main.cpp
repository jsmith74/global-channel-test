
#include "BFGS_Optimization.h"

void generateFilename(std::string& filename,int fileNumber);

int main( int argc, char *argv[] ){

    if(argc != 2){

        std::cout << "./LinearOpticalSimulation [number of encoding states]" << std::endl;

        return 0;

    }

    int encodingStates = std::atoi(argv[1]);

    BFGS_Optimization optimizer(4e-6,20.0,encodingStates);

    for(int i=0;i<60;i++) optimizer.minimize();

    std::string filename;

    generateFilename(filename,encodingStates);

    double globalResult,quantumEntropy,upperBound;

    std::ifstream infile( filename.c_str() );

    infile >> globalResult;

    infile >> quantumEntropy;

    infile >> upperBound;

    infile.close();

    std::ofstream outfile( "globalResults.dat",std::ofstream::app );

    outfile << encodingStates + 1 << "\t" << std::setprecision(16) << globalResult << "\t" << quantumEntropy << "\t" << upperBound << std::endl;

    outfile.close();

    return 0;

}

void generateFilename(std::string& filename,int fileNumber){

    std::stringstream ss;

    ss << fileNumber;

    ss >> filename;

    filename = "globalMonitor_" + filename + ".dat";

    return;

}
