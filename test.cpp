#include <hip/hip_runtime.h>
#include <iostream>
__global__ void test(float* A){
  int i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
  A[i] = A[i]* A[i];
}

int main() {
     float *A;
     hipMalloc(&A , sizeof(float)*1);
     hipLaunchKernelGGL(test, dim3(1) , dim3(1), 0, 0, A);
     hipDeviceSynchronize();
     std::cout << "Done\n";
     hipFree(A);
     return 0;
}
