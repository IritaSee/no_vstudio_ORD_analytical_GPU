#include <iostream> 
#include "cuda_runtime.h"

int main(){
    int count;
    cudaGetDeviceCount(&count);
    printf("cuda device count: %d",count);
}