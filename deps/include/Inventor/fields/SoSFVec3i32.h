#ifndef COIN_SOSFVEC3I32_H
#define COIN_SOSFVEC3I32_H

/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) 1998-2008 by Kongsberg SIM.  All rights reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Kongsberg SIM about acquiring
 *  a Coin Professional Edition License.
 *
 *  See http://www.coin3d.org/ for more information.
 *
 *  Kongsberg SIM, Postboks 1283, Pirsenteret, 7462 Trondheim, NORWAY.
 *  http://www.sim.no/  sales@sim.no  coin-support@coin3d.org
 *
\**************************************************************************/

#include <Inventor/fields/SoSField.h>
#include <Inventor/fields/SoSubField.h>
#include <Inventor/SbVec3i32.h>

class COIN_DLL_API SoSFVec3i32 : public SoSField {
  typedef SoSField inherited;

  SO_SFIELD_HEADER(SoSFVec3i32, SbVec3i32, const SbVec3i32 &);

public:
  static void initClass(void);

  void setValue(int32_t x, int32_t y, int32_t z);
  void setValue(const int32_t xyz[3]);

}; // SoSFVec3i32

#endif // !COIN_SOSFVEC3I32_H
