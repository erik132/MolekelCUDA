#pragma once

#include <math.h>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "../old/molekeltypes.h"
#include "../old/constant.h"

#include "CalcDensDataPack.h"
#include "particleCudaVersions/CudaMolecule.cuh"
#include "particleCudaVersions/CudaMolecularOrbital.cuh"
#include "CalcDensInternalData.h"
#include "../ESLogger.h"

//#define BLOCK_DIM 5

class CalcDensCudaFunction{
private:
	//variables
	

	//methods
	int getSingleGridSize(int elements, int blockSize);

protected:
	//variables
	
	CalcDensInternalData calcData;
	CudaMolecule cudaMolecule;
	CudaMolecule *deviceMolecule;
	CudaMolecularOrbital cudaOrbital;
	CudaMolecularOrbital *deviceOrbital;
	int BLOCK_DIM;
	
	//methods
	dim3 getGridSize();
	vtkImageData* initImageData();
	


public:
	virtual vtkImageData* calcImageData() = 0;
	CalcDensCudaFunction(CalcDensDataPack *data);

	cudaError_t moleculeToDevice();
	cudaError_t orbitalToDevice();
	void deleteDeviceMoleculeData();
	void deleteDeviceOrbitalData();

};