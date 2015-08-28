#include "CudaMolekelAtom.cuh"



void CudaMolekelAtom::setProperties(MolekelAtom *atom){
	int i=0;

	for(int i=0;i<3;i++){
		coord[i] = atom->coord[i];
	}
	shellsSize = atom->Shells.size();
	shells = new CudaShell[shellsSize];

	for(i=0;i<shellsSize;i++){
		
		shells[i].setProperties(&atom->Shells[i]);
	}

}

CudaMolekelAtom::~CudaMolekelAtom(void){
	delete[] shells;

}

cudaError_t CudaMolekelAtom::cpyInternalPointers(void){
	
	cudaError_t status;
	int i;

	for(i=0; i<shellsSize; i++){
		status = shells[i].cpyInternalPointers();
		if(status != cudaSuccess){
			return status;
		}
	}
	
	status = cudaMalloc((void**)&deviceShells,sizeof(CudaShell)*shellsSize);
	if(status != cudaSuccess){
		return status;
	}
	status = cudaMemcpy(deviceShells, shells, sizeof(CudaShell)*shellsSize, cudaMemcpyHostToDevice);

	return status;
}

void CudaMolekelAtom::clearCudaData(void){
	int i;

	for(i=0; i<shellsSize; i++){
		shells[i].clearCudaData();
	}

	cudaFree(deviceShells);
	//CudaMolekelAtom::deviceShells = NULL;
}