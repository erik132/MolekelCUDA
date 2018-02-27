#include "CalcDensCalculateDensity.cuh"
#include "CalcDensDataPack.h"
#include "CalcDensInternalData.h"

#include "particleCudaVersions/CudaMolecule.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


class CalcDensCalculateDensityUnrolled: public CalcDensCalculateDensity{
protected:
	vtkImageData* runComputation() override;

public:
	CalcDensCalculateDensityUnrolled(CalcDensDataPack *data);

};