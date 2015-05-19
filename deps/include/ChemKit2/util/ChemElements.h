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

#ifndef __ELEMENTS__
#define __ELEMENTS__

#define AXIAL_CUTOFF 4.000f

static char *theElements[] = {
   " ",
   "H","He",
   "Li","Be","B","C","N","O","F","Ne",
   "Na","Mg","Al","Si","P","S","Cl","Ar",
   "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
   "Zn","Ga","Ge","As","Se","Br","Kr",
   "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
   "Cd","In","Sn","Sb","Te","I","Xe",
   "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
   "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
   "Au","Hg","Tl","Pb","Bi","Po","At","Rn",
   "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
   "Fm","Md","No","Lr"
};
                                                                                                                   
static char *theElementTable[] = {
   " ",
   "H","HE",
   "LI","BE","B","C","N","O","F","NE",
   "NA","MG","AL","SI","P","S","CL","AR",
   "K","CA","SC","TI","V","CR","MN","FE","CO","NI","CU",
   "ZN","GA","GE","AS","SE","BR","KR",
   "RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG",
   "CD","IN","SN","SB","TE","I","XE",
   "CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY",
   "HO","ER","TM","YB","LU","HF","TA","W","RE","OS","IR","PT",
   "AU","HG","Tl","PB","BI","PO","AT","RN",
   "FR","RA","AC","TH","PA","U","NP","PU","AM","CM","BK","CF","ES",
   "FM","MD","NO","LR"
};

static float covRadius[] = {
   0.00f,
   0.32f,1.60f,
   0.68f,0.35f,0.82f,0.77f,0.75f,0.73f,0.72f,1.12f,
   1.53f,1.43f,1.18f,1.10f,1.06f,1.02f,1.06f,1.54f,
   1.75f,1.81f,1.80f,1.64f,1.48f,1.50f,1.49f,1.47f,1.40f,1.44f,1.35f,
   1.38f,1.37f,1.40f,1.41f,1.37f,1.14f,1.37f,
   1.44f,2.09f,1.78f,1.82f,1.66f,1.61f,1.52f,1.46f,1.41f,1.38f,1.76f,
   1.70f,1.68f,1.66f,1.54f,1.50f,1.37f,1.42f,
   2.57f,2.24f,1.97f,1.98f,1.98f,1.97f,1.96f,1.94f,2.22f,1.93f,1.91f,1.91f,
   1.90f,1.88f,1.87f,2.09f,1.87f,1.75f,1.60f,1.57f,1.54f,1.51f,1.43f,1.43f,
   1.50f,1.54f,1.85f,2.00f,1.60f,1.71f,1.45f,1.74f,
   1.68f,1.68f,1.68f,1.98f,1.98f,1.76f,1.70f,1.70f,1.70f,1.70f,1.70f,1.70f,1.70f,
   1.70f,1.70f,1.70f,1.70f
};

enum stdNames{
   _GLY, _ALA, _VAL, _PHE, _PRO, _MET,
   _ILE, _LEU, _ASP, _GLU, _LYS, _ARG,
   _SER, _THR, _TYR, _HIS, _CYS, _ASN,
   _GLN, _TRP, _AAA, _TTT, _CCC, _GGG,
   _RESNUM
};

static const char* stdName3[_RESNUM] = { "GLY",
          "ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU",
          "ASP", "GLU", "LYS", "ARG",
          "SER", "THR", "TYR", "HIS", "CYS", "ASN", "GLN", "TRP",
          "  A", "  T", "  C", "  G", };
                                                                                                                     
static const char* stdName1[_RESNUM]= { "G",
          "A", "V", "F", "P", "M", "I", "L",
          "D", "E", "K", "R",
          "S", "T", "Y", "H", "C", "N", "Q", "W",
          "A", "T", "C", "G" };

static short aminoCol[_RESNUM][3] = {
{235,235,235},    //_GLY,
{200,200,200},    //_ALA,
{ 15,130, 15},    //_VAL, 
{ 50, 50,170},    //_PHE,
{220,150,130},    //_PRO,
{230,230,  0},    //_MET, 
{ 15,130, 15},    //_ILE,
{ 15,130, 15},    //_LEU,
{230, 10, 10},    //_ASP,
{230, 10, 10},    //_GLU,
{ 20, 90,255},    //_LYS,
{ 20, 90,255},    //_ARG,
{250,150,  0},    //_SER, 
{250,150,  0},    //_THR, 
{ 50, 50,170},    //_TYR, 
{130,130,210},    //_HIS, 
{230,230,  0},    //_CYS, 
{  0,230,230},    //_ASN, 
{  0,230,230},    //_GLN, 
{180, 90,180},    //_TRP,
{ 50, 50, 50},    //_AAA,
{ 50, 50, 50},    //_TTT,
{ 50, 50, 50},    //_CCC,
{ 50, 50, 50}     //_GGG,
};

// based on hydrophobicity scale by Abraham D.J. and Leo A.J.
// Proteins 2:130-152(1987)

static float hydroCol[_RESNUM][3] = {
{1.00f,1.00f,1.00f},    //_GLY
{1.00f,0.83f,0.83f},    //_ALA
{1.00f,0.32f,0.32f},    //_VAL
{1.00f,0.01f,0.01f},    //_PHE
{1.00f,0.50f,0.50f},    //_PRO
{1.00f,0.57f,0.57f},    //_MET
{1.00f,0.04f,0.04f},    //_ILE
{1.00f,0.04f,0.04f},    //_LEU
{0.87f,0.87f,1.00f},    //_ASP
{0.71f,0.71f,1.00f},    //_GLU
{0.00f,0.00f,1.00f},    //_LYS
{0.01f,0.01f,1.00f},    //_ARG
{0.66f,0.66f,1.00f},    //_SER
{0.83f,0.83f,1.00f},    //_THR
{1.00f,0.36f,0.36f},    //_TYR
{1.00f,1.00f,1.00f},    //_HIS
{1.00f,0.77f,0.77f},    //_CYS
{0.46f,0.46f,1.00f},    //_ASN
{0.71f,0.71f,1.00f},    //_GLN
{1.00f,0.00f,0.00f},    //_TRP
{0.50f,0.50f,0.50f},    //_AAA
{0.50f,0.50f,0.50f},    //_TTT
{0.50f,0.50f,0.50f},    //_CCC
{0.50f,0.50f,0.50f}     //_GGG
};


#endif
