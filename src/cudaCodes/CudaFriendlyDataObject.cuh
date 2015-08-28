#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class CudaFriendlyDataObject{
public:
	virtual cudaError_t cpyInternalPointers() = 0;
	virtual void clearCudaData() = 0;


};