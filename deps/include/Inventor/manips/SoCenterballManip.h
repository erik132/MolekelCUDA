#ifndef COIN_SOCENTERBALLMANIP_H
#define COIN_SOCENTERBALLMANIP_H

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

#include <Inventor/manips/SoTransformManip.h>
#include <Inventor/tools/SbLazyPimplPtr.h>

class SoCenterballManipP;

class COIN_DLL_API SoCenterballManip : public SoTransformManip {
  typedef SoTransformManip inherited;

  SO_NODE_HEADER(SoCenterballManip);

public:
  static void initClass(void);
  SoCenterballManip(void);

protected:
  virtual ~SoCenterballManip(void);

  virtual void setDragger(SoDragger * newDragger);

  static void fieldSensorCB(void * f, SoSensor * s);
  static void valueChangedCB(void * f, SoDragger * d);

private:
  SbLazyPimplPtr<SoCenterballManipP> pimpl;

  // NOT IMPLEMENTED:
  SoCenterballManip(const SoCenterballManip & rhs);
  SoCenterballManip & operator = (const SoCenterballManip & rhs);
}; // SoCenterBallManip

#endif // !COIN_SOCENTERBALLMANIP_H
