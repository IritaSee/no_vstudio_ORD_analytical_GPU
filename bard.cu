#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

__global__ void multiply(int *a, int *b, int *c, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    c[i] = a[i] * b[i];
  }
}

int main() {
  // Declare the input and output arrays.
  int *a = (int *)malloc(sizeof(int) * 10);
  int *b = (int *)malloc(sizeof(int) * 10);
  int *c = (int *)malloc(sizeof(int) * 10);

  // Initialize the input arrays.
  for (int i = 0; i < 10; i++) {
    a[i] = i;
    b[i] = i + 1;
  }

  // Check for errors when copying the input arrays to the GPU.
  cudaError_t cudaStatus = cudaMemcpy(a, a, sizeof(int) * 10, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "CUDA memcpy failed: %s\n", cudaGetErrorString(cudaStatus));
    return 1;
  }

  // Check for errors when copying the input arrays to the GPU.
  cudaStatus = cudaMemcpy(b, b, sizeof(int) * 10, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "CUDA memcpy failed: %s\n", cudaGetErrorString(cudaStatus));
    return 1;
  }

  // Execute the kernel function.
  dim3 blockDim(1024, 1);
  dim3 gridDim(1, 1);
  multiply<<<gridDim, blockDim>>>(a, b, c, 10);

  // Check for errors when copying the output array from the GPU.
  cudaStatus = cudaMemcpy(c, c, sizeof(int) * 10, cudaMemcpyDeviceToHost);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "CUDA memcpy failed: %s\n", cudaGetErrorString(cudaStatus));
    return 1;
  }

  // Print the output array.
  for (int i = 0; i < 10; i++) {
    printf("%d * %d = %d\n",a[i],b[i],c[i]);
  }

  // Free the memory allocated for the input and output arrays.
  free(a);
  free(b);
  free(c);

  return 0;
}
