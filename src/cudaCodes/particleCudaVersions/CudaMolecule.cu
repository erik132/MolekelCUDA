#include "CudaMolecule.cuh"


void CudaMolecule::setProperties(Molecule *molecule){
	int i=0;

	nBasisFunctions = molecule->nBasisFunctions;
	atomsSize = molecule->Atoms.size();
	atoms = new CudaMolekelAtom[atomsSize];

	for(i=0; i<CudaMolecule::atomsSize; i++){
		atoms[i].setProperties(&molecule->Atoms[i]);
	}

}


CudaMolecule::~CudaMolecule(void){
	delete[] atoms;
}

cudaError_t CudaMolecule::cpyInternalPointers(void){

	cudaError_t status;
	int i;

	for(i=0; i<atomsSize; i++){
		status = atoms[i].cpyInternalPointers();
		if(status != cudaSuccess){
			return status;
		}
	}
	
	status = cudaMalloc((void**)&deviceAtoms,sizeof(CudaMolekelAtom)*atomsSize);
	if(status != cudaSuccess){
		return status;
	}
	status = cudaMemcpy(deviceAtoms, atoms, sizeof(CudaMolekelAtom)*atomsSize, cudaMemcpyHostToDevice);

	return status;
}

void CudaMolecule::clearCudaData(void){
	int i;

	for(i=0; i<atomsSize; i++){
		atoms[i].clearCudaData();
	}

	cudaFree(deviceAtoms);
}