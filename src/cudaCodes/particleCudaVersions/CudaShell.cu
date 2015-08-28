#include "CudaShell.cuh"



void CudaShell::setProperties(Shell *shell){
	int i=0;
	
	gaussiansSize = shell->gaussians.size();
	nBase = shell->n_base;
	gaussians = new Gauss[gaussiansSize];

	for(i=0;i<gaussiansSize;i++){
		gaussians[i]=shell->gaussians[i];
	}

}

CudaShell::~CudaShell(void){
	delete[] gaussians;
};

cudaError_t CudaShell::cpyInternalPointers(void){

	cudaError_t status;
	
	status = cudaMalloc((void**)&deviceGaussians, sizeof(Gauss)*gaussiansSize);
	if(status != cudaSuccess){
		return status;
	}
	status = cudaMemcpy(deviceGaussians, gaussians, sizeof(Gauss)*gaussiansSize, cudaMemcpyHostToDevice);

	return status;
}

void CudaShell::clearCudaData(void){
	cudaFree(deviceGaussians);
	//CudaShell::deviceGaussians = NULL;
}