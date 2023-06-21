#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

using namespace std;

__global__ int devAr[1];

__global__ void foo()
{
    printf("\n");
    printf("dev before init: %d \n", devAr[0]);
    devAr[0] = 77;
    printf("dev after init: %d \n", devAr[0]);
    printf("\n");
    // printf("cek: %d \n", cek);
}

int main()
{
    int test = 10;
    cout << "host initially: " << devAr[0] << endl;
    devAr[0] = 4;
    cout << "host after init: " << devAr[0] << endl;
    foo << <1, 1 >> >() ;
    cudaDeviceSynchronize();
    cout << "host after foo: " << devAr[0] << endl;
}