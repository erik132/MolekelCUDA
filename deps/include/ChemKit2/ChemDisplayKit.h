/*
 * ChemDisplayKit.h
 *
 *	This file defines the ChemDisplayKit class.
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

#ident "$Revision: 1.5 $"

#ifndef __CHEM_DISPLAY_KIT_H__
#define __CHEM_DISPLAY_KIT_H__

#include <Inventor/SbLinear.h>
#include <Inventor/nodekits/SoBaseKit.h>
#include <ChemKit2/ChemkitBasic.h>

////////////////////////////////////////////////////////////////////
//    Class: ChemDisplayKit
// 
// NOTE TO DEVELOPERS:
//     For info about the structure of ChemDisplayKit:
//     [1] compile: /usr/share/src/Inventor/samples/ivNodeKitStructure
//     [2] type:    ivNodeKitStructure ChemKit.
//     [3] The program prints a diagram of the scene graph and a table with
//         information about each part.
//    
//      A parent node that manages a collection of child nodes
//      into a unit with the following structure:
//
//                            this
//                              |
//           ----------------------------------------------------------
//           |               |               |            |           |
//     "callbackList" "chemDisplayParam" "chemRadii" "chemColor" "chemDisplay"
//
////////////////////////////////////////////////////////////////////


class CHEMKIT_DLL_API ChemDisplayKit : public SoBaseKit {

	SO_KIT_HEADER(ChemDisplayKit);

	SO_KIT_CATALOG_ENTRY_HEADER(chemDisplayParam);
	SO_KIT_CATALOG_ENTRY_HEADER(chemRadii);
	SO_KIT_CATALOG_ENTRY_HEADER(chemColor);
	SO_KIT_CATALOG_ENTRY_HEADER(chemDisplay);

  public:
	// Constructor
	ChemDisplayKit();

  SoINTERNAL public:
	static void initClass();	// initialize the class

  protected:
    /* 
       ~ChemDisplayKit(); was moved from private mode by fabien
       fontaine the 13/12/2000 to avoid g++ compiler warning. 
    */
    // Destructor
    ~ChemDisplayKit();
    
    /* 
       Moved to protected mode by fabien fontaine the 13/12/2000 to avoid 
       g++ compiler warning
  private:
	// Destructor
	~ChemDisplayKit();
    */
};

#endif /* !__CHEM_DISPLAY_KIT_H__ */

