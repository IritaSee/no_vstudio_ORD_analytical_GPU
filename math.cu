#include <stdio.h>

__global__ void multiply(float *a, float *b, float *c, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  printf("%d\n",i);
  if (i < n) {
    c[i] = a[i] * b[i];
  }
}

int main() {
  float A[6] = { 1, 2, 3, 4, 5, 6 };
  float B[6] = { 10, 20, 30, 40, 50, 60 };
  float C[6] = { 0 };
  multiply<<<1,6>>>(A, B, C, 6);
  cudaDeviceSynchronize();
  int i;
  for(i=0;i<6;i++){
    printf("%f * %f = %f\n",A[i],B[i],C[i]);
  }
  return 0;
}