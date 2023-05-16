#include <stdio.h>
__global__ void multiply(float *a, float *b, float *c, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    c[i] = a[i] * b[i];
  }
}


int main(){
   
    float A[5] = { 1, 2, 3, 4, 5 };
    float B[5] = { 10, 20, 30, 40, 50 };
    float C[5];
    multiply<<<1,5>>>(A, B, C, 5);
    int i;
    for(i=0;i<5;i++){
        printf("%f\n",C[i]);
    }
    return 0;
}