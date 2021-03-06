#ifndef MERITFUNCTION_H_INCLUDED
#define MERITFUNCTION_H_INCLUDED

#include <omp.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unsupported/Eigen/KroneckerProduct>
#include "LinearOpticalTransform.h"

class MeritFunction{

    public:

        MeritFunction();
        void setMeritFunction(int intParam);
        double f(Eigen::VectorXd& position);
        int funcDimension;
        void printReport(Eigen::VectorXd& position);
        Eigen::VectorXd setInitialPosition();

    private:

        int ancillaModes,fileNumber;

        double upperBound;

        LinearOpticalTransform LOCircuit;

        void setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv);
        void setInBasis(Eigen::MatrixXi& ancillaBasis,int photons,int modes,Eigen::MatrixXi& inBasis);
        void printArr( int arr[] ,int size );

        double shannonEntropy(Eigen::VectorXd& position);
        double conditionalEntropy(Eigen::MatrixXcd& psiPrime,Eigen::VectorXd& position);
        double vonNeumannEntropy(Eigen::MatrixXcd& psiPrime,Eigen::VectorXd& position);

        Eigen::MatrixXcd genUnitary(Eigen::VectorXd a);
        Eigen::MatrixXcd matrixLog(Eigen::MatrixXcd X);
        Eigen::MatrixXcd matrixExp(Eigen::MatrixXcd X);
        Eigen::MatrixXcd genHermitian(Eigen::VectorXd& a);
        Eigen::VectorXd convertHermittoA(Eigen::MatrixXcd& H);

        void setAncillaVec(Eigen::VectorXd& position,int startPoint);
        void setCompVec(Eigen::VectorXd& position,int startPoint);
        void generateFilename(std::string& filename);

        std::vector<Eigen::MatrixXcd> UAlice;
        Eigen::MatrixXcd UAliceFull,UBob;
        Eigen::VectorXcd psiC,psiA,psi;

        inline int g(const int& n,const int& m);
        inline double doublefactorial(int x);
        int f_ds(int Na,int Ma,int Nb,int Mb);
        int setUpperBound(int N,int Ma,int Mb);
};




inline int MeritFunction::g(const int& n,const int& m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return (int)(doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1))+0.5);
    }
}

inline double MeritFunction::doublefactorial(int x){

    assert(x < 171);

    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;
}

#endif // MERITFUNCTION_H_INCLUDED
