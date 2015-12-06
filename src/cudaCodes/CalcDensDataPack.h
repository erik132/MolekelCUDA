#pragma once
#include "../old/molekelTypes.h"

struct CalcDensDataPack{
	Molecule *mol;
	MolecularOrbital *orbital;
	float *dim;
	int *ncubes;
	int key;
	double minValue, maxValue;
	int datasource;

};