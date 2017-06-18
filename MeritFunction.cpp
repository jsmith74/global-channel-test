#include "MeritFunction.h"

#define PI 3.141592653589

#define INITIAL_CONDITION_RANDOM_DEGREE 2000

void MeritFunction::setMeritFunction(int intParam){

    int wordLength = 1;

    int ancillaPhotons = 0;
    int ancillaModes = 0;

    int photons = 2 * wordLength + ancillaPhotons;
    int modes = 4 * wordLength + ancillaModes;

    Eigen::MatrixXi inBasis, outBasis;

    setToFullHilbertSpace(photons,modes,outBasis);

    setInbasisAndPsi(inBasis,wordLength,ancillaPhotons,ancillaModes);

    LOCircuit.initializeCircuit(inBasis,outBasis);

    funcDimension = modes * modes;

    U = Eigen::MatrixXcd::Identity(modes,modes);

    return;

}

void MeritFunction::setInbasisAndPsi(Eigen::MatrixXi& inBasis,int wordLength,int ancillaPhotons,int ancillaModes){

    int photons = 2 * wordLength + ancillaPhotons;
    int modes = 4 * wordLength + ancillaModes;

    Eigen::VectorXcd psiBell[4];

    Eigen::VectorXi psiBellBasis[4];

    for(int i=0;i<4;i++) psiBell[i].resize(4);
    for(int i=0;i<4;i++) psiBellBasis[i].resize(4);

    psiBellBasis[0] << 1,0,1,0;
    psiBellBasis[1] << 0,1,0,1;
    psiBellBasis[2] << 1,0,0,1;
    psiBellBasis[3] << 0,1,1,0;

    psiBell[0] << 1.0/sqrt(2.0),  1.0/sqrt(2.0),                 0,    0;
    psiBell[1] << 1.0/sqrt(2.0), -1.0/sqrt(2.0),                 0,    0;
    psiBell[2] << 0,                          0,      1.0/sqrt(2.0),   1.0/sqrt(2.0);
    psiBell[3] << 0,                          0,      1.0/sqrt(2.0),  -1.0/sqrt(2.0);

    Eigen::MatrixXcd nonOrthPsiGen(4,4);
    std::complex<double> I(0.0,1.0);

    nonOrthPsiGen <<    0.977274,	-0.074431+0.198483*I,	0,	0,
                        0.074431 +0.198483*I,	0.977274,	0,	0,
                        0,	0,	1.,	0,
                        0,	0,	0,	1.;





    psiBell[0] = nonOrthPsiGen * psiBell[0];

//    nonOrthPsiGen << 0.168868 +0.30214 *I,	-0.2906-0.363023 *I,	-0.512209-0.4431 *I,	0.434652 -0.127838 *I,
//                    -0.516986-0.166408 *I,	-0.0274105+0.627316 *I,	-0.151741-0.120282 *I,	0.522702 +0.00679282 *I,
//                    -0.679094+0.315506 *I,	-0.411681-0.277595 *I,	-0.11694+0.333938 *I,	-0.218589+0.140632 *I,
//                    0.154368 +0.026212 *I,	0.375968 -0.0399698 *I,	-0.473875+0.39251 *I,	0.173005 +0.651139 *I;
//
//    psiBell[1] = nonOrthPsiGen * psiBell[1];
//
//    nonOrthPsiGen << 0.891312 -0.102086 *I,	0.15643 -0.152356 *I,	-0.0825417-0.366368 *I,	0.0111417 +0.0793449 *I,
//                    -0.103664+0.0599398 *I,	0.640907 -0.697416 *I,	-0.0142119+0.29628 *I,	0.0203012 -0.0106868 *I,
//                    -0.0598877-0.418335 *I,	0.210226 +0.0828867 *I,	0.842687 -0.140316 *I,	-0.154087+0.12958 *I,
//                     0.0227513 +0.0411567 *I,	-0.0497852+0.0402866 *I,	0.0168252 +0.202197 *I,	0.46339 +0.858947 *I;
//
//    psiBell[2] = nonOrthPsiGen * psiBell[2];
//
//    nonOrthPsiGen << 0.195441 +0.861182 *I,	0.115601 -0.292968 *I,	0.285417 +0.153207 *I,	-0.118407-0.0449431 *I,
//                     0.260656 +0.155068 *I,	0.145408 +0.900742 *I,	0.129887 +0.152903 *I,	0.0509187 -0.180804 *I,
//                     0.216913 -0.25571 *I,	-0.0587541-0.158796 *I,	-0.11961+0.857281 *I,	-0.275643-0.183511 *I,
//                     0.0603394 +0.110008 *I,	0.134624 +0.146748 *I,	-0.285373+0.15537 *I,	-0.161075+0.901707 *I;
//
//    psiBell[3] = nonOrthPsiGen * psiBell[3];

    psi.resize( std::pow(4,wordLength) , std::pow(4,wordLength) );

    inBasis.resize( std::pow(4,wordLength) , 4 * wordLength + ancillaModes);

    for(int i=0;i<inBasis.rows();i++) for(int j=0;j<ancillaPhotons;j++){

        inBasis(i,j) = 1;

    }

    int ii[wordLength];

    for(int i=0;i<wordLength;i++) ii[i] = 0;

    int k = 0;

    do{

        Eigen::VectorXcd psiColTemp = psiBell[ ii[0] ];

        for(int i=1;i<wordLength;i++){

            psiColTemp = Eigen::kroneckerProduct( psiColTemp,psiBell[ ii[i] ] ).eval();

        }

        psi.col(k) = psiColTemp;

        Eigen::VectorXi basisVectorRowTemp( 4 * wordLength );

        for(int i=0;i<wordLength;i++) basisVectorRowTemp.segment(4*i,4) = psiBellBasis[ ii[i] ];

        inBasis.block(k,ancillaPhotons,1,4*wordLength) = basisVectorRowTemp.transpose();

        k++;

    } while( iterate( ii,wordLength ) );

    for(int i=0;i<std::pow(4,wordLength);i++){

            psi.col(i).normalize();

            //std::cout << "Norm Check: " << std::setprecision(16) << psi.col(i).norm() << std::endl;

    }

    for(int i=0;i<psi.cols();i++) for(int j=i;j<psi.cols();j++){

        //std::cout << i << "\t" << j << "\t" << std::setprecision(16) << (psi.col(i).conjugate().transpose() * psi.col(j)).norm() << std::endl;

    }

    printVonNeumannEntropy(psi);

    return;

}

void MeritFunction::printVonNeumannEntropy(Eigen::MatrixXcd& psi){

    Eigen::MatrixXcd rho = Eigen::MatrixXcd::Zero(psi.rows(),psi.rows());

    for(int i=0;i<psi.cols();i++) rho += (1.0 / psi.cols() ) * (psi.col(i) * psi.col(i).conjugate().transpose());

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;

    ces.compute(rho);

    double S = 0.0;

    for(int i=0;i<rho.rows();i++) S -= ces.eigenvalues()(i) * log2( ces.eigenvalues()(i) );

    std::ofstream outfile("results.dat",std::ofstream::app);

    outfile << "S(rho): " << S << std::endl << std::endl;

    outfile.close();

    return;

}

bool MeritFunction::iterate(int ii[],int wordLength){

    ii[ wordLength-1 ]++;

    for(int i=wordLength-1;i>=1;i--){

        if( ii[i] >= 4 ){

            ii[i] = 0;

            ii[ i-1 ]++;

        }

    }

    if( ii[0] >= 4 ) return false;

    return true;

}

void MeritFunction::printArr( int arr[] ,int size ){

    for(int i=0;i<size;i++) std::cout << arr[i] << "\t";

    std::cout << std::endl;

    return;

}

double MeritFunction::f(Eigen::VectorXd& position){

    U = genUnitary(position);

    LOCircuit.setOmega(U);

    psiPrime = LOCircuit.omega * psi;

    return conditionalEntropy(psiPrime);

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

    U = genUnitary(position);

    LOCircuit.setOmega(U);

    psiPrime = LOCircuit.omega * psi;

    double output = conditionalEntropy(psiPrime);

    std::cout << "H(X:Y): ";

    std::cout << log2( psi.cols() ) - ( 1.0 / psi.cols() ) * output << std::endl;

    std::ofstream outfile("results.dat",std::ofstream::app);

    outfile << "H(X:Y): " << log2( psi.cols() ) - ( 1.0 / psi.cols() ) * output << std::endl << std::endl;

    outfile << U << std::endl << std::endl << std::endl;

    outfile.close();

    return;

}



Eigen::VectorXd MeritFunction::setInitialPosition(){

    std::complex<double> I(0.0,1.0);

    for(int i=0;i<INITIAL_CONDITION_RANDOM_DEGREE;i++){
        Eigen::VectorXd a = Eigen::VectorXd::Random( funcDimension );
        a *= 2000 * PI;
        Eigen::MatrixXcd Utemp( U.rows() , U.cols() );
        Utemp = genUnitary(a);
        U *= Utemp;
    }

    Eigen::MatrixXcd H( U.rows(),U.cols() );

    H = matrixLog(U) / I;

    Eigen::VectorXd a = convertHermittoA(H);

    return a;

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

