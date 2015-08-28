#pragma once

#include <vector>
#include "../../old/molekeltypes.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../CudaFriendlyDataObject.cuh"



class CudaShell: public CudaFriendlyDataObject{
public:
	short int nBase;
	Gauss *gaussians;
	Gauss *deviceGaussians;
	int gaussiansSize;

	~CudaShell();
	void setProperties(Shell *shell);
	cudaError_t cpyInternalPointers();
	void clearCudaData();

};