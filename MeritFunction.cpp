#include "MeritFunction.h"

#define PI 3.141592653589

#define INITIAL_CONDITION_RANDOM_DEGREE 2000


/** ===== Fix Uniform Probabilities ================================ */

//#define USE_UNIFORM_PROBABILITIES

#define ALLOW_PROBABILITIES_TO_MOVE

/** ================================================================ */


void MeritFunction::setMeritFunction(int intParam){

    int photons = 2;
    int modesAlice = 2;
    int modesBob = 2;

    int ancillaPhotons = 0;
        ancillaModes = 0;

    int numberOfStates = intParam;

    upperBound = log2( setUpperBound(photons,modesAlice,modesBob) );

    fileNumber = numberOfStates;

    int totalPhotons = photons + ancillaPhotons;
    int totalModes = modesAlice + modesBob + ancillaModes;

    UAlice.resize(numberOfStates);

    for(int i=0;i<numberOfStates;i++) UAlice.at(i) = Eigen::MatrixXcd::Identity(modesAlice,modesAlice);

    UAliceFull = Eigen::MatrixXcd::Identity( totalModes,totalModes );

    UBob = Eigen::MatrixXcd::Identity( totalModes,totalModes );

    int initialStateDim = g(photons,modesAlice + modesBob);

    int ancillaStateDim = g(ancillaPhotons,ancillaModes);

    #ifdef USE_UNIFORM_PROBABILITIES

        funcDimension = 2 * initialStateDim + 2 * ancillaStateDim + numberOfStates * modesAlice * modesAlice + totalModes * totalModes;

    #endif // USE_UNIFORM_PROBABILITIES

    #ifdef ALLOW_PROBABILITIES_TO_MOVE

        funcDimension = 2 * initialStateDim + 2 * ancillaStateDim + numberOfStates * modesAlice * modesAlice + totalModes * totalModes + numberOfStates + 1;

    #endif // ALLOW_PROBABILITIES_TO_MOVE

    Eigen::MatrixXi inBasis,outBasis,ancillaBasis;

    setToFullHilbertSpace(totalPhotons,totalModes,outBasis);

    setToFullHilbertSpace(ancillaPhotons,ancillaModes,ancillaBasis);

    setInBasis(ancillaBasis,photons,modesAlice+modesBob,inBasis);

    psiC.resize( initialStateDim );
    psiA.resize( ancillaStateDim );

    LOCircuit.initializeCircuit(inBasis,outBasis);

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    int vecSize = UAlice[0].rows() * UAlice[0].rows();

    for(int i=0;i<UAlice.size();i++){

        UAlice[i] = genUnitary( position.segment(i*vecSize,vecSize) );

    }

    UBob = genUnitary( position.segment( UAlice.size() * vecSize , UBob.rows() * UBob.rows() ) );

    setAncillaVec( position, UAlice.size() * vecSize + UBob.rows() * UBob.rows() );

    setCompVec( position, UAlice.size() * vecSize + UBob.rows() * UBob.rows() + 2 * psiA.size() );

    if( psiA.size() > 0 ) psi = Eigen::kroneckerProduct(psiA,psiC);

    else psi = psiC;

    Eigen::MatrixXcd psiPrime( LOCircuit.omega.rows(),UAlice.size()+1 );

    for(int i=0;i<UAlice.size();i++){

        UAliceFull.block( ancillaModes , ancillaModes ,UAlice[i].rows(),UAlice[i].cols()) = UAlice[i];

        Eigen::MatrixXcd U = UAliceFull * UBob;

        LOCircuit.setOmega(U);

        psiPrime.col(i) = LOCircuit.omega * psi;

    }

    LOCircuit.setOmega(UBob);

    psiPrime.col( UAlice.size() ) = LOCircuit.omega * psi;

    return conditionalEntropy(psiPrime,position) - shannonEntropy(position);
    //return -vonNeumannEntropy(psiPrime,position);

}

double MeritFunction::vonNeumannEntropy(Eigen::MatrixXcd& psiPrime,Eigen::VectorXd& position){

    #ifdef USE_UNIFORM_PROBABILITIES

    Eigen::MatrixXcd rho = ( 1.0/psiPrime.cols() ) * psiPrime.col(0) * psiPrime.col(0).conjugate().transpose();

    for(int i=1;i<psiPrime.cols();i++) rho += (1.0/psiPrime.cols()) * psiPrime.col(i) * psiPrime.col(i).conjugate().transpose();

    #endif // USE_UNIFORM_PROBABILITIES

    #ifdef ALLOW_PROBABILITIES_TO_MOVE

    int numberOfStates = UAlice.size() + 1;

    Eigen::VectorXd probabilityVector = position.segment( funcDimension - numberOfStates,numberOfStates );

    for(int i=0;i<probabilityVector.size();i++) probabilityVector(i) *= probabilityVector(i);

    probabilityVector.normalize();

    Eigen::MatrixXcd rho = probabilityVector(0) * psiPrime.col(0) * psiPrime.col(0).conjugate().transpose();

    for(int i=1;i<psiPrime.cols();i++) rho += probabilityVector(i) * psiPrime.col(i) * psiPrime.col(i).conjugate().transpose();

    #endif // ALLOW_PROBABILITIES_TO_MOVE

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;

    ces.compute(rho);

    double output = 0.0;

    for(int i=0;i<ces.eigenvalues().size();i++) if(ces.eigenvalues()(i) > 0) output -= ces.eigenvalues()(i) * log2( ces.eigenvalues()(i) );

    return output;
}


double MeritFunction::shannonEntropy(Eigen::VectorXd& position){

    int numberOfStates = UAlice.size() + 1;

    #ifdef USE_UNIFORM_PROBABILITIES

    return log2( numberOfStates );

    #endif // USE_UNIFORM_PROBABILITIES

    #ifdef ALLOW_PROBABILITIES_TO_MOVE

    Eigen::VectorXd probabilityVector = position.segment( funcDimension - numberOfStates,numberOfStates );

    for(int i=0;i<probabilityVector.size();i++) probabilityVector(i) *= probabilityVector(i);

    probabilityVector.normalize();

    double output = 0.0;

    for(int x=0;x<probabilityVector.size();x++) if(probabilityVector(x) != 0) output -= probabilityVector(x) * log2( probabilityVector(x) );

    return output;

    #endif // ALLOW_PROBABILITIES_TO_MOVE

}

void MeritFunction::printReport(Eigen::VectorXd& position){

    int vecSize = UAlice[0].rows() * UAlice[0].rows();

    for(int i=0;i<UAlice.size();i++){

        UAlice[i] = genUnitary( position.segment(i*vecSize,vecSize) );

    }

    UBob = genUnitary( position.segment( UAlice.size() * vecSize , UBob.rows() * UBob.rows() ) );

    setAncillaVec( position, UAlice.size() * vecSize + UBob.rows() * UBob.rows() );

    setCompVec( position, UAlice.size() * vecSize + UBob.rows() * UBob.rows() + 2 * psiA.size() );

    if( psiA.size() > 0 ) psi = Eigen::kroneckerProduct(psiA,psiC);

    else psi = psiC;

    Eigen::MatrixXcd psiPrime( LOCircuit.omega.rows(),UAlice.size()+1 );

    for(int i=0;i<UAlice.size();i++){

        UAliceFull.block( ancillaModes , ancillaModes ,UAlice[i].rows(),UAlice[i].cols()) = UAlice[i];

        Eigen::MatrixXcd U = UAliceFull * UBob;

        LOCircuit.setOmega(U);

        psiPrime.col(i) = LOCircuit.omega * psi;

    }

    LOCircuit.setOmega(UBob);

    psiPrime.col( UAlice.size() ) = LOCircuit.omega * psi;

    double mutualEntropy = shannonEntropy(position) - conditionalEntropy( psiPrime,position );

    double quantumEntropy = vonNeumannEntropy( psiPrime,position );

    std::cout << "H(X:Y): " << mutualEntropy << std::endl;

    std::cout << "S(rho): " << quantumEntropy << std::endl;

    std::string filename = "";

    generateFilename( filename );

    std::ifstream infile( filename.c_str() );

    if(!infile.is_open()){

        std::ofstream outfile( filename.c_str() );

        outfile << std::setprecision(16) << mutualEntropy << "\t" << quantumEntropy << "\t" << upperBound;

        outfile.close();

    }

    else{

        double globalTest;

        infile >> globalTest;

        infile.close();

        if(mutualEntropy >= globalTest){

            std::ofstream outfile( filename.c_str() );

            outfile << std::setprecision(16) <<  mutualEntropy << "\t" << quantumEntropy << "\t" << upperBound;

            outfile.close();

        }

    }

    return;

}

void MeritFunction::generateFilename(std::string& filename){

    std::stringstream ss;

    ss << fileNumber;

    ss >> filename;

    filename = "globalMonitor_" + filename + ".dat";

    return;

}

void MeritFunction::setCompVec(Eigen::VectorXd& position,int startPoint){

    std::complex<double> I(0.0,1.0);

    for(int i=0;i<psiC.size();i++){

        psiC(i) = position( 2 * i + startPoint );

        psiC(i) *= std::exp( I * position(2 * i + startPoint + 1) );

    }

    psiC.normalize();

    return;

}

void MeritFunction::setAncillaVec(Eigen::VectorXd& position,int startPoint){

    std::complex<double> I(0.0,1.0);

    for(int i=0;i<psiA.size();i++){

        psiA(i) = position( 2 * i + startPoint );

        psiA(i) *= std::exp( I * position( 2 * i + startPoint + 1) );

    }

    psiA.normalize();

    return;

}


Eigen::VectorXd MeritFunction::setInitialPosition(){

    for(int i=0;i<UAlice.size();i++) UAlice.at(i) = Eigen::MatrixXcd::Identity(UAlice.at(i).rows(),UAlice.at(i).cols());

    UAliceFull = Eigen::MatrixXcd::Identity( UAliceFull.rows(),UAliceFull.cols() );

    UBob = Eigen::MatrixXcd::Identity( UBob.rows(),UBob.cols() );

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

    H = matrixLog( UBob ) / I;

    Eigen::VectorXd a = convertHermittoA(H);

    output.segment( aSize * UAlice.size() , a.size() ) = a;

//    #ifdef ALLOW_PROBABILITIES_TO_MOVE
//
//        int numberOfStates = UAlice.size() + 1;
//
//        Eigen::VectorXd probabilityVector = Eigen::VectorXd::Random(numberOfStates);
//
//        probabilityVector.normalize();
//
//        output.segment( funcDimension - numberOfStates,numberOfStates ) = probabilityVector;
//
//    #endif // ALLOW_PROBABILITIES_TO_MOVE

    return output;

}


double MeritFunction::conditionalEntropy(Eigen::MatrixXcd& psiPrime,Eigen::VectorXd& position){

    double output = 0.0;

    #ifdef USE_UNIFORM_PROBABILITIES

    for(int y=0;y<psiPrime.rows();y++){

        double sumTemp =0;

        for(int xx=0;xx<psiPrime.cols();xx++) sumTemp += std::norm( psiPrime(y,xx) );

        for(int x=0;x<psiPrime.cols();x++){

            double pyx = std::norm( psiPrime(y,x) );

            if(pyx != 0) output += pyx * log2( sumTemp / pyx );

        }

    }

    output *= 1.0 / ( UAlice.size() + 1 );

    #endif // USE_UNIFORM_PROBABILITIES

    #ifdef ALLOW_PROBABILITIES_TO_MOVE

    int numberOfStates = UAlice.size() + 1;

    Eigen::VectorXd probabilityVector = position.segment( funcDimension - numberOfStates,numberOfStates );

    for(int i=0;i<probabilityVector.size();i++) probabilityVector(i) *= probabilityVector(i);

    probabilityVector.normalize();

    for(int y=0;y<psiPrime.rows();y++){

        double sumTemp =0;

        for(int xx=0;xx<psiPrime.cols();xx++) sumTemp += probabilityVector(xx) * std::norm( psiPrime(y,xx) );

        for(int x=0;x<psiPrime.cols();x++){

            double pyx = probabilityVector(x) * std::norm( psiPrime(y,x) );

            if(pyx != 0) output += pyx * log2( sumTemp / pyx );

        }

    }

    #endif // ALLOW_PROBABILITIES_TO_MOVE

    return output;

}


MeritFunction::MeritFunction(){



}


void MeritFunction::setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv){

    if(subPhotons==0 && subModes == 0){

        nv.resize(0,0);

        return;

    }

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

    if( ancillaBasis.rows() == 0 ){

        inBasis = compBasis;

        return;

    }

    inBasis.resize( compBasis.rows() * ancillaBasis.rows(), compBasis.cols() + ancillaBasis.cols() );

    for(int i=0;i<ancillaBasis.rows();i++){

        inBasis.block(i*compBasis.rows(),ancillaBasis.cols(),compBasis.rows(),compBasis.cols()) = compBasis;

        for(int j=i*compBasis.rows();j<(i+1)*compBasis.rows();j++){

            inBasis.block( j,0,1,ancillaBasis.cols() ) = ancillaBasis.block(i,0,1,ancillaBasis.cols() );

        }

    }

    return;

}

Eigen::MatrixXcd MeritFunction::genUnitary(Eigen::VectorXd a){

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


int MeritFunction::f_ds(int Na,int Ma,int Nb,int Mb){

    return g(Na,Ma) * std::min( g(Na,Ma),g(Nb,Mb) );

}


int MeritFunction::setUpperBound(int N,int Ma,int Mb){

    int output = 0;

    for(int n=0;n<=N;n++) output += f_ds(n,Ma,N-n,Mb);

    return output;

}


void MeritFunction::printArr( int arr[] ,int size ){

    for(int i=0;i<size;i++) std::cout << arr[i] << "\t";

    std::cout << std::endl;

    return;

}
