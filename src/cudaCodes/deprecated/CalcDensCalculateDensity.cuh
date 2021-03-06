#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "CalcDensCudaFunction.cuh"
#include "CalcDensDataPack.h"
#include "CalcDensInternalData.h"

#include "particleCudaVersions/CudaMolecule.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../ESLogger.h"



class CalcDensCalculateDensity: public CalcDensCudaFunction{


public:
	vtkImageData* calcImageData() override;
	CalcDensCalculateDensity(CalcDensDataPack *data);

};