#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "CalcDensCudaFunction.cuh"
#include "CalcDensDataPack.h"
#include "CalcDensInternalData.h"

#include "particleCudaVersions/CudaMolecule.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../ESLogger.h"



class CalcDensCalcPoint: public CalcDensCudaFunction{


public:
	vtkImageData* calcImageData();
	CalcDensCalcPoint(CalcDensDataPack *data);

};