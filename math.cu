#include <cuda_runtime.h>
#include <iostream>

using namespace std;
// reference github.com/colina118/Multin-cleo

__global__ void multiplyElements(int *a, int *b, int *c) {
  // Get the thread ID.
  int tid = threadIdx.x;

  // Calculate the element index.
  int index = blockIdx.x * blockDim.x + tid;

  // Multiply the elements and store the result.
  c[index] = a[index] * b[index];
}

int main() {
  // Declare the arrays on the host.
  int size = 10;
  int a[size], b[size], c[size];

  // Initialize the arrays.
  for (int i = 0; i < size; i++) {
    a[i] = i;
    b[i] = i * 2;
  }

  // Allocate memory on the device.
  int *d_a, *d_b, *d_c;
  cudaMalloc(&d_a, size * sizeof(int));
  cudaMalloc(&d_b, size * sizeof(int));
  cudaMalloc(&d_c, size * sizeof(int));

  // Copy the arrays to the device.
  cudaMemcpy(d_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, size * sizeof(int), cudaMemcpyHostToDevice);

  // Launch the kernel.
  dim3 grid(1, 1); //grid size
  dim3 block(1024, 1); //block size: 1 is the dimension, each block has 1024 threads
  multiplyElements<<<grid, block>>>(d_a, d_b, d_c);

  // Copy the results back to the host.
  cudaMemcpy(c, d_c, size * sizeof(int), cudaMemcpyDeviceToHost);

  // Print the results.
  for (int i = 0; i < size; i++) {
    cout <<a[i]<<"x"<<b[i]<<"="<< c[i] << endl;
  }

  // Free the memory on the device.
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);

  return 0;
}

