/*
 * ChemCondenser.h
 *
 *     ChemCondenser class: takes a scene graph produced by flattening a graph
 *     and condenses vertices and material indices, if possible, by
 *     removing duplicates.
 *     This assumes that flattening has produced an
 *     SoIndexedTriangleStripSet with lots of 1-triangle "strips", one
 *     for each triangle generated during flattening.
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

/////////////////////////////////////////////////////////////////////////////
//
// ChemCondenser class: takes a scene graph produced by flattening a graph
// and condenses vertices and material indices, if possible, by
// removing duplicates.
//
// This assumes that flattening has produced an
// SoIndexedTriangleStripSet with lots of 1-triangle "strips", one
// for each triangle generated during flattening.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef  __CHEM_CONDENSER_H__
#define  __CHEM_CONDENSER_H__

#include <Inventor/SbBasic.h>

class SbDict;
class SoMFInt32;
class SoIndexedLineSet;

class ChemCondenser {

  public:
    ChemCondenser();
    ~ChemCondenser();

    void condense(SoIndexedLineSet *_lineSet, SbBool doNormals, 
                  SbBool doTexCoords);

    void condenseCoordinates(SoIndexedLineSet *_lineSet);

  private:
    SoIndexedLineSet *lineSet;

    // A vertex
    struct StripVertex {
        int             coordIndex;     // Coordinate index
        int             normIndex;      // Normal index
        int             texCoordIndex;  // Texture coordinate index
        int             mtlIndex;       // Material index
        int             uniqueID;       // Unique ID for edge hashing
        int             curLine;
        SbPList         lines;          // lines the vertex belongs to
        StripVertex     *next;          // For making lists
    };

    // A line
    struct StripLine {
        StripVertex     *v[2];          // Vertices of line
        SbBool          isUsed;         // TRUE when used in strip
    };

    int         numLines;
    int         numVertices;
    int         *vertexMap;
    StripVertex *vertices;
    StripLine   *lines;
    SbPList     *vertexPtrList;

    SbBool      haveMaterials;
    SbBool      haveNormals;
    SbBool      haveTexCoords;

    // Removes duplicate coordinates, normals, or texture coordinates,
    // updating the indices in the lineSet.
    void condenseCoordinates();
    void condenseNormals();
    void condenseTextureCoordinates();

    // Condenses the material indices if possible
    void condenseMaterials();

    // Returns TRUE if the two sets of indices are the same
    SbBool sameIndices(const SoMFInt32 *indexField1,
                       const SoMFInt32 *indexField2);

    void createVertexList();
    void createLineList();
    void createPolyLines();
    void adjustIndices();


};

#endif /* __CHEM_CONDENSER_H__ */