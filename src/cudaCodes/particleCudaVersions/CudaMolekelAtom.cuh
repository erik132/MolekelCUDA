#pragma once

#include <vector>

#include "../../old/molekeltypes.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../CudaFriendlyDataObject.cuh"
#include "CudaShell.cuh"

class CudaMolekelAtom: public CudaFriendlyDataObject{
public:
	float coord[3];
	CudaShell *shells;
	CudaShell *deviceShells;
	int shellsSize;

	~CudaMolekelAtom();
	void setProperties(MolekelAtom *atom);
	cudaError_t cpyInternalPointers();
	void clearCudaData();

};