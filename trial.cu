#include <stdio.h>
#include <cuda.h>
__global__ void dkernel(){
    printf("Hello Cuda");
}
 // compile like usual c code (nvcc trial.cu -o trial) will create trial.exe
int main(){
    dkernel<<<1,1>>>();
    return 0;
}