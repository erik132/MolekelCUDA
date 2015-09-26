#pragma once

#include "../particleCudaVersions/CudaMolecule.cuh"
#include "../particleCudaVersions/CudaMolecularOrbital.cuh"

#include "../../old/constant.h"
#include "../../old/molekeltypes.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


static __device__ double calcChiCalcPoint(CudaMolecularOrbital *orbital, CudaMolecule *molecule, float x, float y, float z){

	double radial_part;
	float xa, ya, za, ra2;  /* atomic units !! */
	int atom, shell, gauss, i;

	const int nBasisFunctions = molecule->nBasisFunctions;
	const int atomsSize = molecule->atomsSize;
	int shellsSize, gaussSize, count=0;
	Gauss tempGauss;

	double cp[10], result = 0; //temporary cells will be used to gather gaussian calculations in them to later add them to the result.
	
	
	for (atom=0; atom<atomsSize; atom++) {
		xa = (x - molecule->deviceAtoms[atom].coord[0]) * _1_BOHR;
		ya = (y - molecule->deviceAtoms[atom].coord[1]) * _1_BOHR;
		za = (z - molecule->deviceAtoms[atom].coord[2]) * _1_BOHR;

		ra2 = xa*xa + ya*ya + za*za;    /* cutoff-distance ? */
		shellsSize=molecule->deviceAtoms[atom].shellsSize;

		for (shell=0; shell<shellsSize; shell++) {
			gaussSize = molecule->deviceAtoms[atom].deviceShells[shell].gaussiansSize;
			switch(molecule->deviceAtoms[atom].deviceShells[shell].nBase){
			case 1  :        /*** S-orbital ***/
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = exp(-ra2*tempGauss.exponent);
					*cp += tempGauss.coeff * radial_part;
				}
				result += cp[0]*orbital->deviceCoefficients[count];
				count++;
				break;

			case 4 :        /*** SP-orbital ***/
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = exp(-ra2*tempGauss.exponent);
					*cp     += tempGauss.coeff * radial_part;
					*(cp+1) += tempGauss.coeff2 * xa * radial_part;
					*(cp+2) += tempGauss.coeff2 * ya * radial_part;
					*(cp+3) += tempGauss.coeff2 * za * radial_part;
				}
				for(i=0; i<4; i++)
					result += cp[i]* orbital->deviceCoefficients[count + i];
				count += 4;
				break;

			case 3  :        /*** P-orbital ***/
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = tempGauss.coeff * exp(-ra2*tempGauss.exponent);
					*cp     += xa * radial_part;
					*(cp+1) += ya * radial_part;
					*(cp+2) += za * radial_part;
					}
				for(i=0; i<3; i++)
					result += cp[i]* orbital->deviceCoefficients[count + i];
				count += 3;
				break;

			case 5  :        /*** D-orbital (5) ***/
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = tempGauss.coeff * exp(-ra2*tempGauss.exponent);
					*cp    += 0.288675135 *
						   (2*za*za - xa*xa - ya*ya) * radial_part;
					*(cp+3) += 0.5 * (xa*xa - ya*ya) * radial_part;
					*(cp+4) += xa * ya * radial_part;
					*(cp+1) += xa * za * radial_part;
					*(cp+2) += ya * za * radial_part;
				}
				for(i=0; i<5; i++)
					result += cp[i]* orbital->deviceCoefficients[count + i];
				count += 5;
				break;

			case 6  :        /*** D-orbital (6) ***/
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = tempGauss.coeff * exp(-ra2*tempGauss.exponent);
					*cp     += radial_part * xa * xa * 0.57735027;
					*(cp+1) += radial_part * ya * ya * 0.57735027;
					*(cp+2) += radial_part * za * za * 0.57735027;
					*(cp+3) += radial_part * xa * ya;
					*(cp+4) += radial_part * xa * za;
					*(cp+5) += radial_part * ya * za;
				}
				for(i=0; i<6; i++)
					result += cp[i]* orbital->deviceCoefficients[count + i];
				count += 6;
				break;

			case 7  :        /*** F-orbital (7) ***/
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = exp(-ra2*tempGauss.exponent) * tempGauss.coeff;
					*cp     += radial_part * za * (5. * za * za - 3. * ra2)/* * k */;
					*(cp+1) += radial_part * xa * (5. * za * za - ra2)/* * k */;
					*(cp+2) += radial_part * ya * (5. * za * za - ra2)/* * k */;
					*(cp+3) += radial_part * za * (xa * xa - ya * ya)/* * k */;
					*(cp+4) += radial_part * xa * ya * za;
					*(cp+5) += radial_part * (xa * xa * xa - 3. * xa * ya * ya)/* * k */;
					*(cp+6) += radial_part * (3. * xa * xa * ya - ya * ya * ya)/* * k */;
				}
				for(i=0; i<7; i++)
					result += cp[i]* orbital->deviceCoefficients[count + i];
				count += 7;
				break;

			case 10 :        /*** F-orbital (10) ***/
				   /* correct order ??? */
				for (gauss=0; gauss<gaussSize; gauss++) {
					tempGauss = molecule->deviceAtoms[atom].deviceShells[shell].deviceGaussians[gauss];
					radial_part = tempGauss.coeff * exp(-ra2*tempGauss.exponent);
					*cp     += radial_part * xa * xa * xa * .25819889;
					*(cp+1) += radial_part * ya * ya * ya * .25819889;
					*(cp+2) += radial_part * za * za * za * .25819889;
					*(cp+3) += radial_part * xa * xa * ya * .57735027;
					*(cp+4) += radial_part * xa * xa * za * .57735027;
					*(cp+5) += radial_part * xa * ya * ya * .57735027;
					*(cp+6) += radial_part * ya * ya * za * .57735027;
					*(cp+7) += radial_part * xa * za * za * .57735027;
					*(cp+8) += radial_part * ya * za * za * .57735027;
					*(cp+9) += radial_part * xa * ya * za;
				}
				for(i=0; i<10; i++)
					result += cp[i]* orbital->deviceCoefficients[count + i];
				count += 10;
				break;

			} /* end of switch */
			for(i=0; i<10; i++)
				cp[i] = 0;
		} /* end of loop over the shells (for(sp...) */
	} /* end of loop over the atoms (for(ap...)*/
	return result;
}