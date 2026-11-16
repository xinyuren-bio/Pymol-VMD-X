
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#ifndef _H_FontType
#define _H_FontType

#include"Font.h"
#include"TypeFace.h"

struct CFontType : public CFont {
  CTypeFace* TypeFace;

  ~CFontType() override;

  CFontType(PyMOLGlobals* G, unsigned char* dat, unsigned int len);
  /* 先建空字体，再由调用方填入从文件读出的字形。 */
  explicit CFontType(PyMOLGlobals* G) : CFont(G), TypeFace(nullptr) {}
  const char* RenderOpenGL(const RenderInfo* info, const char* text, float size,
      const float* rpos, bool needSize, short relativeMode, bool shouldRender,
      CGO* shaderCGO) override;
  const char* RenderOpenGLFlat(const RenderInfo* info, const char* text, float size,
      const float* rpos, bool needSize, short relativeMode, bool shouldRender,
      CGO* shaderCGO) override;
  const char* RenderRay(CRay* ray, const char* text, float size, const float* rpos,
      bool needSize, short relativeMode) override;
};

CFont *FontTypeNew(PyMOLGlobals * G, unsigned char *dat, unsigned int len);
/* 用字体文件创建界面字体。文件打不开时返回空。 */
CFont *FontTypeNewFromFile(PyMOLGlobals * G, const char *path);

#endif
