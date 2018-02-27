#include "CalcDensCalculateDensityUnrolled.cuh"

__global__ void calculateDensityUnrolled(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, float *densities, double *results, int offsetx){
	double result = 0;
	int indexZ = threadIdx.z + (blockDim.z*blockIdx.z);
	int	indexY = threadIdx.y + (blockDim.y*blockIdx.y);
	int	indexX = threadIdx.x + (blockDim.x*(blockIdx.x + offsetx));
	float x,y,z;

	if(indexX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
		
		x = internalData.dim0 + indexX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim4 + indexZ*internalData.dz;
		
		//result = calcChiCalculateDensity(densities,orbital, molecule, x, y, z, internalData.densityLength);
		results[indexX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] = result;
	}
}

CalcDensCalculateDensityUnrolled::CalcDensCalculateDensityUnrolled(CalcDensDataPack *data): CalcDensCalculateDensity(data){
	
}
