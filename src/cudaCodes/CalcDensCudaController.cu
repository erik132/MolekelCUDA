#include "CalcDensCudaController.cuh"
#include "gputimer.h"
#include "CalcDensCalcPoint.cuh"

void CalcDensCudaController::getOrbitalFunction(CalcDensDataPack *data){

	switch(data->mol->alphaOrbital[0].flag) {
		  case GAMESS_ORB :
		  case HONDO_ORB  :
		  case GAUSS_ORB  : 
			  CalcDensCudaController::calcFunction = new CalcDensCalcPoint(data); break;

		  case MOS_ORB   :
		  case ZINDO_ORB  :
		  case PRDDO_ORB  : 
			  /*funct = calc_prddo_point;*/ break;
		  case MLD_SLATER_ORB  : 
			  /*funct = calc_sltr_point;*/ break;
	}
}

void CalcDensCudaController::getElectroDensityFunction(CalcDensDataPack *data){
	
	switch(data->mol->alphaOrbital[0].flag) {
	  case GAMESS_ORB :
	  case HONDO_ORB  :
	  case GAUSS_ORB  :
	   /*if (!generate_density_matrix(mol, key)) {
		fprintf(stderr, "Can't generate the density matrix!\n");
		strcpy(timestring, "Can't generate the density matrix!");
		CalcDensCudaController::calcFunction = NULL;
	   }
	   else {
		printf("density matrix generated...\n");
		esl.logMessage("function calculate_density activated");
		funct = calculate_density;
	   }*/
	   break;


	  case MOS_ORB   :
	  case ZINDO_ORB  :
	  case PRDDO_ORB  : /*funct = calc_prddo_density;*/ break;
	  case MLD_SLATER_ORB  : /*funct = calc_sltr_density;*/ break;
	}

}

/*
	TODO: when density matrix generation is impossible the program is terminated with a null pointer error. Find a more polite wait to end it.
*/

void CalcDensCudaController::getSpinDensityFunction(CalcDensDataPack *data){

	switch(data->mol->alphaOrbital[0].flag) {
	  case GAMESS_ORB :
	  case HONDO_ORB  :
	  case GAUSS_ORB  :
	   /*if (!mol->alphaBeta) {
		funct = calculateSomo;
	   }
	   else if(!generate_density_matrix(mol, key)) {
		fprintf(stderr, "Can't generate the density matrix!\n");
		strcpy(timestring, "Can't generate the density matrix!");
		CalcDensCudaController::calcFunction = NULL;
	   }
	   else {
		printf("density matrix generated...\n");
		funct = calculate_density;
	   }
	   break;*/
	  case MOS_ORB   :
	  case ZINDO_ORB  :
	  case PRDDO_ORB  : /*funct = calc_prddo_spindensity;*/ break;
	  case MLD_SLATER_ORB  : /*funct = calc_sltr_spindensity;*/ break;
	}
}

vtkImageData* CalcDensCudaController::vtkProcessCalc(CalcDensDataPack *data){
	
	ESLogger esl("CalcDensCudaController.txt");
	GpuTimer gputimer;
	char buffer[100];
	vtkImageData *returnData;


	esl.logMessage("function started");

	switch(data->key) {
	   case CALC_ORB  :
			CalcDensCudaController::getOrbitalFunction(data);
		break;

	   case EL_DENS :
			CalcDensCudaController::getElectroDensityFunction(data);
		break;

	   case SPIN_DENS :
		   CalcDensCudaController::getSpinDensityFunction(data);
		break;

	   case MEP :
		/*funct = calc_mep;*/
		break;
	  }

	/*if(calcFunction->getDataPack()->mol == NULL){
		esl.logMessage("data pack molecule was NULL");
	}else{
		esl.logMessage("data pac molecule was NOT NULL");
	}*/
	if (calcFunction != NULL){
		gputimer.Start();
		returnData = calcFunction->calcImageData();
		gputimer.Stop();
		sprintf(buffer,"total gpu elapsed time: %f",gputimer.Elapsed());
		esl.logMessage(buffer);
		return returnData;
		
	}else{
		return NULL;
	}
}

CalcDensCudaController::~CalcDensCudaController(){
	delete calcFunction;
}