#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "../old/constant.h"

#include "CalcDensCudaFunction.cuh"
#include "CalcDensDataPack.h"
#include "CalcDensCalculateDensity.cuh"

#include "../ESLogger.h"



class CalcDensCudaController{
private:
	CalcDensCudaFunction *calcFunction;
	void getOrbitalFunction(CalcDensDataPack *data);
	void getElectroDensityFunction(CalcDensDataPack *data);
	void getSpinDensityFunction(CalcDensDataPack *data);
public:
	vtkImageData* vtkProcessCalc(CalcDensDataPack *data);
	CalcDensCudaController();
	~CalcDensCudaController();

};