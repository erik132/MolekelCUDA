#include "CalcDensCudaService.h"

#include "CalcDensCudaController.cuh"

vtkImageData* CalcDensCudaService::vtkProcessCalc(CalcDensDataPack *data){

	CalcDensCudaController controller;

	return controller.vtkProcessCalc(data);

}