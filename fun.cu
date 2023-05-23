#include <stdio.h>
#define N 10

__global__ void fun(){
    printf("%d\n", threadIdx.x);
    }

int main(){
    fun<<<1,N>>>();
    cudaThreadSynchronize();
    return 0;
}