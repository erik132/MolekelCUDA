#pragma once

#include <math.h>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "../old/molekeltypes.h"
#include "../old/constant.h"

#include "CalcDensDataPack.h"
#include "particleCudaVersions/CudaMolecule.cuh"
#include "../ESLogger.h"

#define BLOCK_DIM 5

class CalcDensCudaFunction{
protected:
	vtkImageData* initImageData();
	int ncub[3];
	float x, y, z, dx, dy, dz;
	CalcDensDataPack *calcData;
	CudaMolecule cudaMolecule;
	CudaMolecule *deviceMolecule;

public:
	virtual vtkImageData* calcImageData() = 0;
	CalcDensCudaFunction(CalcDensDataPack *data);
	CalcDensDataPack *getDataPack();

	cudaError_t moleculeToDevice();
	void deleteDeviceData();

};