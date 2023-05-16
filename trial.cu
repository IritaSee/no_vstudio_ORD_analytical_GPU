#include <stdio.h>
#include <cuda.h>
__global__ void dkernel(){
    printf("Hello Cuda\n");
}
 // compile like usual c code (nvcc trial.cu -o trial) will create trial.exe
int main(){
    //syntax:
    //[global or host void name]<<<kernel invocation, kernel threads>>>();
    dkernel<<<100,1000>>>(); //these are called kernels -> when called, are executed 32 times in parallel by 32 different CUDA threads, as opposed to only once like regular C functions.
    /*
    A kernel is defined using the __global__ declaration specifier and the number of CUDA threads 
    that execute that kernel for a given kernel call is specified using a new <<<...>>> execution
    configuration syntax (see C Language Extensions). Each thread that executes the kernel is given a unique 
    thread ID that is accessible within the kernel through the built-in threadIdx variable.
    */
    cudaThreadSynchronize();
    return 0;
}