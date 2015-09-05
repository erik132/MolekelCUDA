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

#define BLOCK_DIM 5

class CalcDensCudaFunction{
private:
	cudaError_t cpyOrbital();


protected:
	vtkImageData* initImageData();
	CalcDensInternalData calcData;
	CudaMolecule cudaMolecule;
	CudaMolecule *deviceMolecule;
	CudaMolecularOrbital cudaOrbital;
	CudaMolecularOrbital *deviceOrbital;


public:
	virtual vtkImageData* calcImageData() = 0;
	CalcDensCudaFunction(CalcDensDataPack *data);

	cudaError_t moleculeToDevice();
	void deleteDeviceData();

};