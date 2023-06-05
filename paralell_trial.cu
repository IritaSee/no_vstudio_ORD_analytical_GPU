#include <iostream>

 __global__ void child_kernel() {printf("hello\n");}

 __global__ void addKernel(int *c, const int *a, const int *b)
 {
     child_kernel << <1, 1 >> > ();
     int i = threadIdx.x;
     c[i] = a[i] + b[i];
 }

 int main (){
    int c[5];
    int a[5] = {1,2,3,4,5};
    int b[5] = {2,3,4,5,6};
    addKernel<<<1,1>>>(c,a,b);
    cudaDeviceSynchronize();
 }