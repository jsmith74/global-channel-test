
#include "BFGS_Optimization.h"


int main(){

    BFGS_Optimization optimizer(4e-4,20.0,1);

    for(int i=0;i<300;i++) optimizer.minimize();

    return 0;

}

