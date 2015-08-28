#pragma once

#include <vector>

#include "../../old/molekeltypes.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../CudaFriendlyDataObject.cuh"
#include "CudaMolekelAtom.cuh"

class CudaMolecule: public CudaFriendlyDataObject{
public:
	CudaMolekelAtom *atoms;
	CudaMolekelAtom *deviceAtoms;
	int atomsSize;
	int nBasisFunctions;

	~CudaMolecule();
	void setProperties(Molecule *molecule);
	cudaError_t cpyInternalPointers();
	void clearCudaData();
	
};