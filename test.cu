#include <stdio.h>

__global__ void fillArray(int* arr, int size) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x + 1; // Compute the unique thread ID
    if (tid <= size) {
        arr[tid - 1] = tid; // Fill the array with numbers from 1 to 1024
    }
}

int main() {
    const int arraySize = 1024;
    const int blockSize = 256;
    const int gridSize = (arraySize + blockSize - 1) / blockSize;

    int* hostArray = new int[arraySize];
    int* deviceArray;

    // Allocate memory on the device
    cudaMalloc((void**)&deviceArray, arraySize * sizeof(int));

    // Launch the kernel to fill the array
    fillArray<<<gridSize, blockSize>>>(deviceArray, arraySize);

    // Copy the result back to the host
    cudaMemcpy(hostArray, deviceArray, arraySize * sizeof(int), cudaMemcpyDeviceToHost);

    // Print the filled array
    printf("Filled Array: ");
    for (int i = 0; i < arraySize; ++i) {
        printf("%d ", hostArray[i]);
    }
    printf("\n");

    // Free allocated memory
    delete[] hostArray;
    cudaFree(deviceArray);

    return 0;
}

