#include <stdio.h>

// Use nvcc compilier to run *.cu


__global__ void hello_cuda(){
    printf("Hello CUDA!\n");
}


int main(int argc, char ** argv){
    //  One block containig one thread.
    hello_cuda<<<1, 1>>>();
    //  To see output.
    cudaDeviceSynchronize();

    return 0;
}