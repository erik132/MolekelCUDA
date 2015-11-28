#include "CalcDensCalculateDensity.cuh"

vtkImageData* CalcDensCalculateDensity::calcImageData(){
	ESLogger esl("CalcDensCalculateDensity.txt");
	esl.logMessage("function started");

	const int resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	cudaError_t status;
	char buffer[1000];
	
	sprintf(buffer, "resultsLength is %d" , resultsLength);
	esl.logMessage(buffer);

	status = moleculeToDevice();
	
	if(status == cudaSuccess){
		esl.logMessage("Molecule copied with success");
	}else{
		esl.logMessage("Molecule copy failed");
	}

	deleteDeviceMoleculeData();
	return NULL;
}

CalcDensCalculateDensity::CalcDensCalculateDensity(CalcDensDataPack *data): CalcDensCudaFunction(data){

}