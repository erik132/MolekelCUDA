#pragma once

struct CalcDensInternalData{
	int ncub0, ncub1, ncub2;
	float dim0, dim1, dim2, dim3, dim4, dim5;
	float dx, dy, dz;
	double minValue, maxValue;
	int datasource;
	int densityLength; //nr of elements in density matrix
	int offsetx, offsety, offsetz;

};