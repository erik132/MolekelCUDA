#pragma once

#include "../../old/molekeltypes.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../CudaFriendlyDataObject.cuh"

class CudaMolecularOrbital: public CudaFriendlyDataObject{
public:
	double *coefficients;
	double *deviceCoefficients;
	int coefficientsSize;

	~CudaMolecularOrbital();
	void setProperties(MolecularOrbital *molOrb, Molecule *mol);
	cudaError_t cpyInternalPointers();
	void clearCudaData();
	
};