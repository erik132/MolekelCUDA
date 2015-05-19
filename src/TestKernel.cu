#include "TestKernel.h"
#include "ESLogger.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>


__global__ void addKernel(int *result, int adder){
	
	atomicAdd(result,adder);
}

void activateCuda(){
	ESLogger esl("cudaLog.txt");
	int result=0,*d_result;
	char buffer[50];
	esl.logMessage("cuda activated");

	cudaMalloc(&d_result,sizeof(int));
	cudaMemcpy(d_result,&result,sizeof(int),cudaMemcpyHostToDevice);
	addKernel<<<1,128>>>(d_result,2);
	cudaMemcpy(&result,d_result,sizeof(int),cudaMemcpyDeviceToHost);
	cudaFree(d_result);
	
	sprintf(buffer, "cuda result was %d", result);
	esl.logMessage(buffer);
	

}

