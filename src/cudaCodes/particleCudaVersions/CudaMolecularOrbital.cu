#include "CudaMolecularOrbital.cuh"


void CudaMolecularOrbital::setProperties(MolecularOrbital *molOrb, Molecule *mol){

	coefficients = molOrb->coefficient;
	coefficientsSize = mol->nBasisFunctions;
}


cudaError_t CudaMolecularOrbital::cpyInternalPointers(void){
	cudaError_t status;

	status = cudaMalloc((void**)&deviceCoefficients,sizeof(double)*coefficientsSize);
	if(status!=cudaSuccess){
		return status;
	}

	status = cudaMemcpy(deviceCoefficients, coefficients, sizeof(double)*coefficientsSize, cudaMemcpyHostToDevice);
	return status;
}

void CudaMolecularOrbital::clearCudaData(void){
	cudaFree(deviceCoefficients);

}

CudaMolecularOrbital::~CudaMolecularOrbital(void){


}