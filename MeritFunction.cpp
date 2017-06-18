#include "MeritFunction.h"

#define PI 3.141592653589

#define INITIAL_CONDITION_RANDOM_DEGREE 2000

void MeritFunction::setMeritFunction(int intParam){

    int photons = 3;
    int modesAlice = 2;
    int modesBob = 2;

    int ancillaPhotons = 1;
    int ancillaModes = 2;

    int numberOfStates = 8;

    int totalPhotons = photons + ancillaPhotons;
    int totalModes = modesAlice + modesBob + ancillaModes;

    UAlice.resize(numberOfStates);

    for(int i=0;i<numberOfStates;i++) UAlice.at(i) = Eigen::MatrixXcd::Identity(modesAlice,modesAlice);

    UBob = Eigen::MatrixXcd::Identity( totalModes,totalModes );

    int initialStateDim = g(photons,modesAlice + modesBob);

    int ancillaStateDim = g(ancillaPhotons,ancillaModes);

    funcDimension = 2 * initialStateDim + 2 * ancillaStateDim + numberOfStates * modesAlice * modesAlice + totalModes * totalModes;

    Eigen::MatrixXi inBasis,outBasis,ancillaBasis;

    setToFullHilbertSpace(totalPhotons,totalModes,outBasis);

    setToFullHilbertSpace(ancillaPhotons,ancillaModes,ancillaBasis);

    setInBasis(ancillaBasis,photons,modesAlice+modesBob,inBasis);

    std::cout << "InBasis:\n" << inBasis << std::endl << std::endl;
    std::cout << "OutBasis: \n" << outBasis << std::endl << std::endl;

    LOCircuit.initializeCircuit(inBasis,outBasis);

    return;

}



Eigen::VectorXd MeritFunction::setInitialPosition(){

    Eigen::VectorXd output = PI * Eigen::VectorXd::Random(funcDimension);

    int aSize;

    std::complex<double> I(0.0,1.0);

    for(int i=0;i<UAlice.size();i++){

        for(int j=0;j<INITIAL_CONDITION_RANDOM_DEGREE;j++){

            Eigen::VectorXd a = Eigen::VectorXd::Random( UAlice[i].cols() * UAlice[i].cols() );
            a *= 2000 * PI;
            Eigen::MatrixXcd Utemp( UAlice[i].rows() , UAlice[i].cols() );
            Utemp = genUnitary(a);
            UAlice[i] *= Utemp;

        }

        Eigen::MatrixXcd H( UAlice[i].rows() , UAlice[i].cols() );

        H = matrixLog(UAlice[i]) / I;

        Eigen::VectorXd a = convertHermittoA(H);

        output.segment( i * ( a.size() ), a.size() ) = a;

        aSize = a.size();

    }

    for(int j=0;j<INITIAL_CONDITION_RANDOM_DEGREE;j++){

        Eigen::VectorXd a = Eigen::VectorXd::Random( UBob.cols() * UBob.cols() );
        a *= 2000 * PI;
        Eigen::MatrixXcd Utemp( UBob.rows() , UBob.cols() );
        Utemp = genUnitary(a);
        UBob *= Utemp;

    }

    Eigen::MatrixXcd H( UBob.rows() , UBob.cols() );

    H = matrixLog(UBob) / I;

    Eigen::VectorXd a = convertHermittoA(H);

    output.segment( aSize * UAlice.size() , a.size() ) = a;

    std::cout << output << std::endl << std::endl;

    assert(false);

    return output;

}



void MeritFunction::printArr( int arr[] ,int size ){

    for(int i=0;i<size;i++) std::cout << arr[i] << "\t";

    std::cout << std::endl;

    return;

}

double MeritFunction::f(Eigen::VectorXd& position){

    return 2.0;

}

double MeritFunction::conditionalEntropy(Eigen::MatrixXcd& psiPrime){

    double output = 0.0;

    for(int y=0;y<psiPrime.rows();y++){

        double sumTemp =0;

        for(int xx=0;xx<psiPrime.cols();xx++) sumTemp += std::norm( psiPrime(y,xx) );

        for(int x=0;x<psiPrime.cols();x++){

            double pyx = std::norm( psiPrime(y,x) );

            if(pyx != 0) output += pyx * log2( sumTemp / pyx );

        }

    }

    return output;

}


void MeritFunction::printReport(Eigen::VectorXd& position){


    return;

}


MeritFunction::MeritFunction(){



}


void MeritFunction::setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv){
    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );
    return;;
}


void MeritFunction::setInBasis(Eigen::MatrixXi& ancillaBasis,int photons,int modes,Eigen::MatrixXi& inBasis){

    Eigen::MatrixXi compBasis;

    setToFullHilbertSpace(photons,modes,compBasis);

    inBasis.resize( compBasis.rows() * ancillaBasis.rows(), compBasis.cols() + ancillaBasis.cols() );

    for(int i=0;i<ancillaBasis.rows();i++){

        inBasis.block(i*compBasis.rows(),ancillaBasis.cols(),compBasis.rows(),compBasis.cols()) = compBasis;

        for(int j=i*compBasis.rows();j<(i+1)*compBasis.rows();j++){

            inBasis.block( j,0,1,ancillaBasis.cols() ) = ancillaBasis.block(i,0,1,ancillaBasis.cols() );

        }

    }

    return;

}

Eigen::MatrixXcd MeritFunction::genUnitary(Eigen::VectorXd& a){

    return matrixExp(genHermitian(a));

}

Eigen::MatrixXcd MeritFunction::matrixExp(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();
    std::complex<double> I(0.0,1.0);

                                                                //THIS NEEDS TO BE AS EFFICIENT AS POSSIBLE
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::MatrixXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::MatrixXcd result(matrixSize,matrixSize);
    result = exp(I*evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result=result+exp(I*evalues(j))*sylvester[j];
    }

    return result;
}


Eigen::VectorXd MeritFunction::convertHermittoA(Eigen::MatrixXcd& H){

    int Hsize = H.rows();

    Eigen::VectorXd output(Hsize*Hsize);

    int rowIndex = 0;

    int outputIndex =0;

    for(int i=0;i<Hsize;i++){

        output(outputIndex) = real(H(rowIndex,rowIndex));
        outputIndex++;

        for(int j=rowIndex+1;j<Hsize;j++){
            output(outputIndex) = sqrt(norm(H(rowIndex,j)));
            outputIndex++;
            output(outputIndex) = arg(H(rowIndex,j));
            outputIndex++;
        }

        rowIndex++;
    }

    return output;

}

Eigen::MatrixXcd MeritFunction::genHermitian(Eigen::VectorXd& a){

    std::complex<double> I(0.0,1.0);
    int Hsize = sqrt(a.size());
    Eigen::MatrixXcd m(Hsize,Hsize);
    int extractIndex=0;                                     //REWRITE THIS FUNCTION IT NEEDS TO BE EFFICIENT- EIGEN SHOULD HAVE A STANDARD ONE

    for(int i=0;i<Hsize;i++){

        m(i,i)=a(extractIndex);
        extractIndex++;

        for(int j=i;j<Hsize;j++){

            if(i!=j){

                m(i,j) = a(extractIndex) * exp(I*a(extractIndex+1));
                m(j,i) = a(extractIndex) * exp(-I*a(extractIndex+1));
                extractIndex++;
                extractIndex++;

            }

        }

    }

    return m;
}


Eigen::MatrixXcd MeritFunction::matrixLog(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXcd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::MatrixXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::MatrixXcd result(matrixSize,matrixSize);
    result = log(evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result = result + log(evalues(j))*sylvester[j];
    }

    return result;
}

