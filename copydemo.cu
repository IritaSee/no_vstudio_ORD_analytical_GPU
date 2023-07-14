#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

using namespace std;

double devAr[4];
double *d_devAr;
__global__ void foo(double *devAr)
{
    // printf("\n");
    // printf("dev before init: %d \n", devAr[0]);
    devAr[0] = 0.77;
    // printf("dev after init: %d \n", devAr[0]);
    // printf("\n");
}

int main()
{
    double test[4];
    devAr[0] = 0.0;
    devAr[1] = 0.33;
    devAr[2] = 0.66;
    devAr[3] = 0.99;
    cout << "host initially: " << devAr[0] << endl;
    devAr[0] = 7.0;
    cout << "host after init: " << devAr[0] << endl;
    printf("malloc: %d\n",cudaMalloc((void**)&d_devAr, 4*sizeof(double))); 
    // cudaMalloc((void**)&d_devAr, 5*sizeof(double)); 
    printf("\nstatus: %d\n",cudaMemcpy(d_devAr, devAr, 4*sizeof(double), cudaMemcpyHostToDevice));
    //foo << <1, 1 >> >(*d_devAr);
    cudaDeviceSynchronize();
    printf("\nstatus: %d\n",cudaMemcpy(test, d_devAr, 4*sizeof(double), cudaMemcpyDeviceToHost));
    for (int a=0; a<4;a++){
        cout << "after foo: " << test[a] << endl;
    }
    
}