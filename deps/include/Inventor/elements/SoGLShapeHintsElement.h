#ifndef COIN_SOGLSHAPEHINTSELEMENT_H
#define COIN_SOGLSHAPEHINTSELEMENT_H

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

#include <Inventor/elements/SoShapeHintsElement.h>

class COIN_DLL_API SoGLShapeHintsElement : public SoShapeHintsElement {
  typedef SoShapeHintsElement inherited;

  SO_ELEMENT_HEADER(SoGLShapeHintsElement);
public:
  static void initClass(void);
protected:
  virtual ~SoGLShapeHintsElement();

public:
  virtual void init(SoState * state);

  virtual void push(SoState * state);
  virtual void pop(SoState * state,
                   const SoElement * prevTopElement);

  static void forceSend(SoState * const state, const SbBool twoside);
  static void forceSend(SoState * const state,
                        const SbBool ccw, const SbBool cull);
  static void forceSend(SoState * const state, const SbBool ccw,
                        const SbBool cull, const SbBool twoside);

protected:
  virtual void setElt(VertexOrdering vertexOrdering,
                      ShapeType shapeType, FaceType faceType);
private:
  SoState * state;
};

#endif // !COIN_SOGLSHAPEHINTSELEMENT_H
