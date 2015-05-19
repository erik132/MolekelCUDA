/*
 * ChemDisplayParamElement.h
 *
 *	This file defines the ChemDisplayParamElement class.
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

#ifndef  __CHEM_DISPLAY_PARAM_ELEMENT_H__
#define  __CHEM_DISPLAY_PARAM_ELEMENT_H__

#include <ChemKit2/ChemDisplayParam.h>
#include <Inventor/elements/SoReplacedElement.h>
#include <ChemKit2/ChemkitBasic.h>

//////////////////////////////////////////////////////////////////////////////
//
//  Class: ChemDisplayParamElement
//
//  Element that stores the current chemical system information.
//
//////////////////////////////////////////////////////////////////////////////

SoEXTENDER class CHEMKIT_DLL_API ChemDisplayParamElement : public SoReplacedElement {

    SO_ELEMENT_HEADER(ChemDisplayParamElement);

  public:
    virtual void    init(SoState *state);

    // Set
    static void     set(SoState *state, SoNode *node, ChemDisplayParam *cdp);

    // Get
    static ChemDisplayParam *get(SoState *state);

    // Returns the top (current) instance of the element in the state
    static const ChemDisplayParamElement * getInstance(SoState *state);

    virtual SbBool matches(const SoElement *elt) const;

    virtual SoElement *copyMatchInfo() const;

    virtual void print(FILE *fp) const;

  SoINTERNAL public:
    // Initializes the ChemDisplayParamElement class
    static void     initClass();

  protected:
    ChemDisplayParam *data;
    virtual ~ChemDisplayParamElement() {}
};

#endif /* !__CHEM_DISPLAY_PARAM_ELEMENT_H__ */
