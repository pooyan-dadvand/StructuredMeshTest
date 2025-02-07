
#include <iostream>
#include <math.h>
// Kernel function to add the elements of two arrays
__global__
void LHSMult(int n, float length, float *x, float *y, float *conductivities, float *densities)
{
    const int block_first_element_index = blockIdx.x * blockDim.x * blockDim.y * blockDim.z;
    const int block_first_node_index = blockIdx.x * (blockDim.x + 1) * (blockDim.y + 1) * (blockDim.z + 1);
    const int thread_index = threadIdx.x + blockDim.x * threadIdx.y + blockDim.x * blockDim.y * threadIdx.z;
    const int element_index = block_first_element_index + thread_index;
    const int node_0_index = block_first_node_index + threadIdx.x + threadIdx.y * ( blockDim.x + 1) + threadIdx.z * ( blockDim.x + 1) * ( blockDim.x + 1);
    const int node_1_index=node_0_index+1;
    const int node_2_index=node_0_index+blockDim.x + 2;
    const int node_3_index=node_0_index+blockDim.x + 1;
    const int node_4_index=node_0_index+(blockDim.x + 1) * (blockDim.x + 1);
    const int node_5_index=node_0_index+1+(blockDim.x + 1) * (blockDim.x + 1);
    const int node_6_index=node_0_index+blockDim.x + 2+(blockDim.x + 1) * (blockDim.x + 1);
    const int node_7_index=node_0_index+blockDim.x + 1+(blockDim.x + 1) * (blockDim.x + 1);

    float d = conductivities[element_index];
    float rho = densities[element_index];

    const float X_0 = x[node_0_index];
    const float X_1 = x[node_1_index];
    const float X_2 = x[node_2_index];
    const float X_3 = x[node_3_index];
    const float X_4 = x[node_4_index];
    const float X_5 = x[node_5_index];
    const float X_6 = x[node_6_index];
    const float X_7 = x[node_7_index];

    const float x0 = d*length;
    const float l3 = length*length*length;
    const float x2 = rho*l3;
    const float x3 = (1.0/3.0)*x0 + (1.0/27.0)*x2;
    const float x4 = (1.0/12.0)*x0;
    const float x5 = (1.0/216.0)*rho*l3 - x4;
    const float x6 = (1.0/108.0)*rho*l3 - x4;
    const float x7 = X_5*x6;
    const float x8 = X_7*x6;
    const float x9 = (1.0/54.0)*x2;
    const float x10 = X_1*x9;
    const float x11 = X_3*x9;
    const float x12 = x10 + x11 + x7 + x8;
    const float x13 = X_2*x6 + X_4*x9;
    const float x14 = X_4*x6;
    const float x15 = X_6*x6;
    const float x16 = X_0*x9;
    const float x17 = X_2*x9;
    const float x18 = x14 + x15 + x16 + x17;
    const float x19 = X_3*x6 + X_5*x9;
    const float x20 = X_0*x6 + X_6*x9;
    const float x21 = X_1*x6 + X_7*x9;
    const float x22 = x19 + x21;
    const float x23 = x13 + x20;

    y[node_0_index]+= X_0*x3 + X_6*x5 + x12 + x13 ;
    y[node_1_index]+= X_1*x3 + X_7*x5 + x18 + x19 ;
    y[node_2_index]+= X_2*x3 + X_4*x5 + x12 + x20 ;
    y[node_3_index]+= X_3*x3 + X_5*x5 + x18 + x21 ;
    y[node_4_index]+= X_2*x5 + X_4*x3 + x15 + x16 + x22 ;
    y[node_5_index]+= X_3*x5 + X_5*x3 + x10 + x23 + x8 ;
    y[node_6_index]+= X_0*x5 + X_6*x3 + x14 + x17 + x22 ;
    y[node_7_index]+= X_1*x5 + X_7*x3 + x11 + x23 + x7 ;

}

__global__
void SetValue(float* rVector, const int rSize, const float rValue) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

    for (int i = index; i < rSize; i+=stride) {
        rVector[i] = rValue;
    }
}

template <class T>
double norm_2(const T* rVector, const int rSize) {
    double result = 0.0;
    for (int i = 0; i < rSize; i++) {
        result += rVector[i] * rVector[i];
    }
    return sqrt(result);
}

float* CreateArray(int size, float value) {
  float* array;
  cudaError_t error = cudaMallocManaged(&array, size * sizeof(float));

  if (error != cudaSuccess) {
    std::cerr << "CUDA error: " << cudaGetErrorString(error) << std::endl;
    exit(EXIT_FAILURE);
  }
  int blockSize = 256;
  int numBlocks = (size + blockSize - 1) / blockSize;
  SetValue<<<numBlocks, blockSize>>>(array, size, value);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
  }

  return array;
}


template<int N>
static double MultiplyBenchmark() {
  std::cout << N << "\t";
  // Adding a timer to measure the time taken for initialization
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start);

  // Setting the model parameters
  const int size = 256 / N;
  // double element_size = 0.216; // 0.1 is the size of the domain
  const int num_nodes_per_block = (N + 1) * (N + 1) * (N + 1);
  const int num_elements_per_block = N * N * N;
  const int num_blocks = size * size * size;
  const int num_nodes = num_blocks * num_nodes_per_block;
  const int num_elements = num_blocks * num_elements_per_block;

  std::cout << int(num_nodes/1e6) << "M\t\t";

  // allocate memory for the arrays
  float *x = CreateArray(num_nodes, 1.0f);
  float *y = CreateArray(num_nodes, 0.0f);
  float *conductivity = CreateArray(num_elements, 100.0f);
  float *density = CreateArray(num_elements, 1000.0f);

  // print the time taken for initialization
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  std::cout << milliseconds << " ms" << "\t";

  dim3 gridDim(num_blocks, 1, 1);
  dim3 blockDim(N, N, N);

  // Adding a timer to measure the time taken for the multiplication
  cudaEvent_t startMult, stopMult;
  cudaEventCreate(&startMult);
  cudaEventCreate(&stopMult);
  cudaEventRecord(startMult);

  const float length = 0.1;
  for (int i = 0; i < 100; i++){
    LHSMult<<<gridDim, blockDim>>>(num_elements, length, x, y, conductivity, density);
    // cudaDeviceSynchronize();
  }

  cudaEventRecord(stopMult);
  cudaEventSynchronize(stopMult);
  float millisecondsMult = 0;
  cudaEventElapsedTime(&millisecondsMult, startMult, stopMult);
  std::cout << millisecondsMult << " ms\t\t"; // returning the norm of the y array
  std::cout << norm_2(y, num_nodes) << std::endl;
  return norm_2(y, num_nodes);
}

int main(void)
{
  std::cout << "N\tN. Nodes\tBuild\t\tMultiply\tNorm2" << std::endl;

  MultiplyBenchmark<4>();
  MultiplyBenchmark<8>();
  // MultiplyBenchmark<16>();
  // MultiplyBenchmark<32>();
  // MultiplyBenchmark<64>();
  // int N = 1<<20;
  // float *x, *y,*conductivity, *density;

  // // Allocate Unified Memory – accessible from CPU or GPU
  // cudaMallocManaged(&x, N*sizeof(float));
  // cudaMallocManaged(&y, N*sizeof(float));


  // // initialize x and y arrays on the host
  // for (int i = 0; i < N; i++) {
  //   x[i] = 1.0f;
  //   y[i] = 2.0f;
  // }

  // // Run kernel on 1M elements on the GPU
  // int blockSize = 256;
  // int numBlocks = (N + blockSize - 1) / blockSize;
  // LHSMult<<<numBlocks, blockSize>>>(N, x, y,conductivity,density);

  // // Wait for GPU to finish before accessing on host
  // cudaDeviceSynchronize();

  // // Check for errors (all values should be 3.0f)
  // float maxError = 0.0f;
  // for (int i = 0; i < N; i++)
  //   maxError = fmax(maxError, fabs(y[i]-3.0f));
  // std::cout << "Max error: " << maxError << std::endl;

  // // Free memory
  // cudaFree(x);
  // cudaFree(y);
  
  // return 0;
}