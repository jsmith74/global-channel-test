#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

int main(){

    for(int X=1;X<20;X++){

        std::string filename;

        std::stringstream ss;

        ss << X;
        ss >> filename;

        filename = "globalMonitor_" + filename + ".dat";

        std::ifstream infile( filename.c_str() );

        if( !infile.is_open() ) continue;

        double mutualEntropy,vonNeumannEntropy,upperBound;

        infile >> mutualEntropy;
        infile >> vonNeumannEntropy;
        infile >> upperBound;

        infile.close();

        std::ofstream outfile("globalResultsProc.dat",std::ofstream::app);

        outfile << X + 1 << std::setprecision(16) << "\t" << mutualEntropy << "\t" << vonNeumannEntropy << "\t" << upperBound << std::endl;

        outfile.close();

    }

    return 0;

}
