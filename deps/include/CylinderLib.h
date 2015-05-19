/*
 * CylinderLib.h
 *
 *     Cylinder Library functions
 *
 * Copyright 1996, 1997, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 *
 * UNPUBLISHED -- Rights reserved under the copyright laws of the United
 * States.   Use of a copyright notice is precautionary only and does not
 * imply publication or disclosure.
 *
 * U.S. GOVERNMENT RESTRICTED RIGHTS LEGEND:
 * Use, duplication or disclosure by the Government is subject to restrictions
 * as set forth in FAR 52.227.19(c)(2) or subparagraph (c)(1)(ii) of the Rights
 * in Technical Data and Computer Software clause at DFARS 252.227-7013 and/or
 * in similar or successor clauses in the FAR, or the DOD or NASA FAR
 * Supplement.  Contractor/manufacturer is Silicon Graphics, Inc.,
 * 2011 N. Shoreline Blvd. Mountain View, CA 94039-7311.
 *
 * THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY
 * INFORMATION OF SILICON GRAPHICS, INC. ANY DUPLICATION, MODIFICATION,
 * DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY
 * PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF SILICON
 * GRAPHICS, INC.
 */
/**************************************************************************\
 *
 * OpenMOIV - C++ library for molecular visualization using Inventor.
 * Copyright (C) 2001-2003 Universitat Pompeu Fabra - Barcelona (Spain)
 * 
 * Developers: Interactive Technology Group (GTI)
 * Team: Josep Blat, Eduard Gonzalez, Sergi Gonzalez,
 *       Daniel Soto, Alejandro Ramirez, Oscar Civit.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details (see the file 
 * LICENSE.LGPL at the root directory).
 *
 * REMARK: This library is a derived product.
 *         You need also to accept all other applicable licenses.
 * 
 * Homepage: http://www.tecn.upf.es/openMOIV/
 * Contact:  openmoiv@upf.es
 *
\**************************************************************************/

#ident "$Revision: 1.3 $"

#ifndef __CYLINDER_LIB_H__
#define __CYLINDER_LIB_H__

#include <Inventor/SbColor.h>  
#include <Inventor/SbLinear.h>  

void cylinderLibUnitCylinder(SbBool doCaps, int numSections, int numSides,
    int &cylinderNormCount, int &cylinderVertCount,
	SbVec2f *ringCoords,
    SbVec3f *&cylinderNormals, SbVec3f *&cylinderVertices);
                
void cylinderLibUnitCylinderCaps(int numSections, int numSides,
    int &cylinderCapNormCount, int &cylinderCapVertCount,  
    int &cylinderBottomStartNorm, int &cylinderBottomStartVert,
	SbVec2f *ringCoords,
    SbVec3f *&cylinderCapNormals, SbVec3f *&cylinderCapVertices);
        
void cylinderLibGenerateCylinderDisplayList(const SbMatrix &theMatrix,
    const SbBool doCaps, const SbBool doNormals, const SbBool doTextures,
    const int numSections, const int numSides,
    const SbVec3f *cylinderNormals, const SbVec3f *cylinderVertices);

void cylinderLibGenerateColoredCylinderDisplayList(const SbMatrix &theMatrix,
    const SbBool doCaps, const SbBool doNormals, const SbBool doTextures,
    const int numSections, const int numSides, 
    const SbColor color1, const SbColor color2,
    const SbVec3f *cylinderNormals, const SbVec3f *cylinderVertices);
                
void cylinderLibGenerateCylinderCapDisplayList(const SbMatrix &theMatrix,
    const SbBool bottom, const SbBool doNormals, const SbBool doTextures,  
    const int numSections, const int numSides,
    const int cylinderBottomStartNorm, const int cylinderBottomStartVert,
    const SbVec3f *cylinderCapNormals, const SbVec3f *cylinderCapVertices);
                
void cylinderLibComputeDepth(float complexity,
    int &numSections, int &numSides);
         
void cylinderLibComputeRingCoords(int numSides, SbVec2f *&ringCoords);

#endif /* !__CYLINDER_LIB_H__ */
