#include <stdio.h>
#include <cuda.h>
__global__ void dkernel(){
    int t = threadIdx.x;
    int T = blockIdx.x;
    printf(" tread index: %d  block index: %d\n",t,T);
}
 // compile like usual c code (nvcc trial.cu -o trial) will create trial.exe
int main(){
    //syntax:
    //[global or host void name]<<<block, thread per block>>>();
    dkernel<<<3,10>>>(); //these are called kernels -> when called, are executed 32 times in parallel by 32 different CUDA threads, as opposed to only once like regular C functions.
    /*
    A kernel is defined using the __global__ declaration specifier and the number of CUDA threads 
    that execute that kernel for a given kernel call is specified using a new <<<...>>> execution
    configuration syntax (see C Language Extensions). Each thread that executes the kernel is given a unique 
    thread ID that is accessible within the kernel through the built-in threadIdx variable.
    */
    cudaDeviceSynchronize();
    /*
    The CUDA code you have provided is syntactically correct. 
    However, it will not print anything to the console. 
    This is because the kernel is not actually running on the GPU. 
    To run the kernel on the GPU, you need to call the cudaDeviceSynchronize() function. 
    This function will ensure that the kernel has finished executing before the main() function returns.
    */
    //cudaThreadSynchronize(); -> deprecated
    return 0;
/*
There are 4 variables each thread can access which contains information 
about the organization of threads and current thread. They are:

threadIdx : Id of the current thread.
blockIdx : Id of the current block.
blockDim : Size of each dimension of the current block.
gridDim : Size of each dimension of the current grid.

All of these are dim3 structure. we can use dot notation to access variable 
x,y,z which contains the information of the corresponding dimension. Example: threadIdx.x
*/

}