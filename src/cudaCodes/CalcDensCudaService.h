#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "CalcDensDataPack.h"

class CalcDensCudaService{

public:
	vtkImageData* vtkProcessCalc(CalcDensDataPack *data);

};