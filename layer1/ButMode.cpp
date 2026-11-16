
/* 
   A* -------------------------------------------------------------------
   B* This file contains source code for the PyMOL computer program
   C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
   D* -------------------------------------------------------------------
   E* It is unlawful to modify or remove this copyright notice.
   F* -------------------------------------------------------------------
   G* Please see the accompanying LICENSE file for further information. 
   H* -------------------------------------------------------------------
   I* Additional authors of this source file include:
   -* 
   -* 
   -*
   Z* -------------------------------------------------------------------
*/
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"
#include"MemoryDebug.h"
#include "main.h"
#include "Base.h"
#include "ButMode.h"
#include "Scene.h"
#include "Util.h"
#include "Ortho.h"
#include "Setting.h"
#include "P.h"
#include "Text.h"
#include "Menu.h"
#include "CGO.h"
#include "Movie.h"

/* 侧栏图例：18 DIP 行距、8 DIP 左右边距，适配默认约 220 DIP 宽。 */
#define cButModeLineHeight DIP2PIXEL(18)
#define cButModeLeftMargin DIP2PIXEL(8)
#define cButModeTopMargin DIP2PIXEL(4)
#define cButModeBottomMargin DIP2PIXEL(4)

struct CButMode : public Block {
  CodeType Code[cButModeCount + 1] {};
  int NCode {};
  int Mode[cButModeInputCount] {};
  int NBut {};
  float Rate { 0.0f };
  float RateShown { 0.0f };
  float Samples { 0.0f }, Delay { 0.0f };
  /* 列/行标签：哑光灰 */
  float TextColor1[3] = { 0.55f, 0.55f, 0.58f };
  /* 动作名：柔和白 */
  float TextColor2[3] = { 0.88f, 0.88f, 0.88f };
  /* 当前模式与选择：浅青强调色 */
  float TextColor3[3] = { 0.45f, 0.82f, 0.85f };
  int DeferCnt { 0 };
  float DeferTime { 0.0f };

  CButMode(PyMOLGlobals * G) : Block(G) {}

  int click(int button, int x, int y, int mod) override;
  void draw(CGO *orthoCGO) override;
  bool fastDraw(CGO* orthoCGO) override;
};


/*========================================================================*/
Block *ButModeGetBlock(PyMOLGlobals * G)
{
  CButMode *I = G->ButMode;
  {
    return (I);
  }
}

int ButModeGetHeight(PyMOLGlobals * G)
{
  /* 展开：约 10 行×18 + 呼吸边距；紧凑：模式/选择/状态三行。 */
  if(SettingGetGlobal_b(G, cSetting_mouse_grid))
    return DIP2PIXEL(196);
  else
    return DIP2PIXEL(60);
}


/*========================================================================*/
int ButModeGet(PyMOLGlobals * G, int button)
{
  CButMode *I = G->ButMode;
  if((button >= 0) && (button < I->NBut)) {
    return I->Mode[button];
  }
  return 0;
}

void ButModeSet(PyMOLGlobals * G, int button, int action)
{
  CButMode *I = G->ButMode;
  if((button >= 0) && (button < I->NBut) && (action >= 0) && (action < I->NCode)) {
    I->Mode[button] = action;
    OrthoDirty(G);
  }
}


/*========================================================================*/
void ButModeSetRate(PyMOLGlobals * G, float interval)
{
  CButMode *I = G->ButMode;

  if(interval >= 0.001F) {      /* sub-millisecond, defer... */
    if(I->DeferCnt) {
      interval = (interval + I->DeferTime) / (I->DeferCnt + 1);
      I->DeferCnt = 0;
      I->DeferTime = 0.0F;
    }
    I->Delay -= interval;
    if(interval < 1.0F) {
      I->Samples *= 0.95 * (1.0F - interval);
      I->Rate *= 0.95 * (1.0F - interval);
    } else {
      I->Samples = 0.0F;
      I->Rate = 0.0F;
    }

    I->Samples++;
    I->Rate += 1.0F / interval;
  } else {
    I->DeferCnt++;
    I->DeferTime += interval;
  }
}


/*========================================================================*/
void ButModeResetRate(PyMOLGlobals * G)
{
  CButMode *I = G->ButMode;
  I->Samples = 0.0;
  I->Rate = 0.0;
  I->RateShown = 0;
  I->Delay = 0;
}


/*========================================================================*/
void ButModeFree(PyMOLGlobals * G)
{
  DeleteP(G->ButMode);
}


/*========================================================================*/
int CButMode::click(int button, int x, int y, int mod)
{
  /* 自底向上：第 0 行状态、第 1 行选择（随 18 DIP 行高缩放）；之上循环鼠标模式。 */
  int dy = (y - rect.bottom) / cButModeLineHeight;
  //  int dx = (x - block->rect.left);
  int forward = 1;
  // TAKEN OUT : BB 12/11 : Mouse position should not have
  // an effect on whether the mouse ring goes forwards or
  // backwards.
  // forward = (dx > ((block->rect.right - block->rect.left) / 2));
  //
  /*  register CButMode *I=block->G->ButMode; */
  if(button == P_GLUT_BUTTON_SCROLL_BACKWARD || button == P_GLUT_RIGHT_BUTTON)
    forward = !forward;
  if(mod == cOrthoSHIFT)
    forward = !forward;
  if(dy < 2) {
    if(ButModeTranslate(m_G, P_GLUT_SINGLE_LEFT, 0) != cButModePickAtom) {
      if(!forward) {
        PLog(m_G, "cmd.mouse('select_backward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse select_backward,quiet=1");
      } else {
        PLog(m_G, "cmd.mouse('select_forward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse select_forward,quiet=1");
      }
    }
  } else {
    if(button == P_GLUT_RIGHT_BUTTON) {
      MenuActivate0Arg(m_G,x,y,x,y,false,"mouse_config");
    } else {
      if(!forward) {
        PLog(m_G, "cmd.mouse('backward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse backward,quiet=1");
      } else {
        PLog(m_G, "cmd.mouse('forward')", cPLog_pym);
        OrthoCommandIn(m_G, "mouse forward,quiet=1");
      }
    }
  }
  return (1);
}

static bool ButModeDrawFastImpl(Block * block, short definitely , CGO *orthoCGO);

/*========================================================================*/
/* 鼠标模式的内部名字仍是英文，这里只改屏幕上看到的字。 */
static const char *zh_mouse_mode(const char *s)
{
  if (!s)
    return "";
  if (!strcmp(s, "3-Button Viewing")) return "三键查看";
  if (!strcmp(s, "3-Button Editing")) return "三键编辑";
  if (!strcmp(s, "3-Button Motions")) return "三键动画";
  if (!strcmp(s, "3-Button Lights")) return "三键光照";
  if (!strcmp(s, "3-Button Maestro")) return "三键大师";
  if (!strcmp(s, "3-Button All Modes")) return "三键全部";
  if (!strcmp(s, "2-Button Viewing")) return "两键查看";
  if (!strcmp(s, "2-Btn. Selecting")) return "两键选择";
  if (!strcmp(s, "2-Button Editing")) return "两键编辑";
  if (!strcmp(s, "2-Button Lights")) return "两键光照";
  if (!strcmp(s, "1-Button Viewing")) return "单键查看";
  return s;
}

/* 按 Mode 下标取动作文案，越界或未映射时返回占位。 */
static const char *butmode_action(const CButMode *I, int idx)
{
  if(!I || idx < 0 || idx >= I->NBut)
    return "-";
  int mode = I->Mode[idx];
  if(mode < 0 || mode >= I->NCode)
    return "-";
  return mode == cButModeNone ? "-" : I->Code[mode];
}

/* 选择粒度中文标签。 */
static const char *zh_sele_mode(int m)
{
  switch(m) {
  case 0: return "原子";
  case 1: return "残基";
  case 2: return "链";
  case 3: return "区段";
  case 4: return "对象";
  case 5: return "分子";
  case 6: return "α碳";
  default: return "";
  }
}

/* 中性分隔线，兼容 CGO 与立即模式 GL。 */
static void butmode_draw_hline(PyMOLGlobals *G, int x0, int x1, int y,
                               const float *col, CGO *orthoCGO)
{
  if(!(G && G->HaveGUI && G->ValidContext))
    return;
  if(orthoCGO) {
    CGOColorv(orthoCGO, col);
    CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
    CGOVertex(orthoCGO, (float) x0, (float) y, 0.f);
    CGOVertex(orthoCGO, (float) x1, (float) y, 0.f);
    CGOVertex(orthoCGO, (float) x0, (float) (y + 1), 0.f);
    CGOVertex(orthoCGO, (float) x1, (float) (y + 1), 0.f);
    CGOEnd(orthoCGO);
  } else {
#ifndef PURE_OPENGL_ES_2
    glColor3fv(col);
    glBegin(GL_LINES);
    glVertex2i(x0, y);
    glVertex2i(x1, y);
    glEnd();
#endif
  }
}

/* 画一行四个鼠标列动作（左/中/右/滚轮），全部 TextDrawStrAt。 */
static void butmode_draw_action_row(PyMOLGlobals *G, const CButMode *I,
                                    const int *idx4, const int *col_x, int y,
                                    const float *act_col, CGO *orthoCGO)
{
  TextSetColor(G, act_col);
  for(int c = 0; c < 4; c++) {
    int idx = idx4[c];
    if(idx < 0)
      continue;
    TextDrawStrAt(G, butmode_action(I, idx), col_x[c], y, orthoCGO);
  }
}

void CButMode::draw(CGO* orthoCGO)
{
  const float *labelColor = TextColor1;
  const float *actionColor = TextColor2;
  const float *accentColor = TextColor3;
  const float *sectionColor = TextColor;
  /* 分隔线：略亮于更深背景的中性灰 */
  const float dividerColor[3] = { 0.30f, 0.30f, 0.32f };
  CButMode *I = this;

  if(m_G->HaveGUI && m_G->ValidContext && ((rect.right - rect.left) > 6)) {
    /* 这一块说明用中文字体，对象名字仍走原来的点阵字 */
    TextOverlayFontScope zh_font(m_G);
    if(SettingGet<InternalGUIMode>(m_G, cSetting_internal_gui_mode) == InternalGUIMode::Default) {
      if (orthoCGO)
        CGOColorv(orthoCGO, BackColor);
#ifndef PURE_OPENGL_ES_2
      else
        glColor3fv(BackColor);
#endif
      fill(orthoCGO);
      drawLeftEdge(orthoCGO);
    } else {
      drawLeftEdge(orthoCGO);
      if (orthoCGO)
        CGOColor(orthoCGO, .5f, .5f, .5f);
#ifndef PURE_OPENGL_ES_2
      else
        glColor3f(0.5, 0.5, 0.5);
#endif
      drawTopEdge();
      /* 非 Default 背景时跟 overlay 色，避免发灰发糊 */
      actionColor = OrthoGetOverlayColor(m_G);
      sectionColor = actionColor;
      labelColor = actionColor;
    }

    int x = rect.left + cButModeLeftMargin;
    int xRight = rect.right - cButModeLeftMargin;
    int content = xRight - x;
    if(content < 1)
      content = 1;

    /* 底两行固定：状态 / 选择；click 仍用 dy<2 命中这两行。 */
    int y_sele = rect.bottom + cButModeBottomMargin + cButModeLineHeight;
    int y = (rect.top - cButModeLineHeight) - cButModeTopMargin;

    /* 首列宽度随侧栏自适应，保证 >=176 DIP 时四列仍可读 */
    int col0 = DIP2PIXEL(36);
    if(content >= DIP2PIXEL(200))
      col0 = DIP2PIXEL(44);
    else if(content < DIP2PIXEL(160))
      col0 = DIP2PIXEL(30);
    if(col0 > content / 3)
      col0 = content / 3;
    int rest = content - col0;
    if(rest < 4)
      rest = 4;
    int colw = rest / 4;
    int col_x[4] = {
      x + col0,
      x + col0 + colw,
      x + col0 + 2 * colw,
      x + col0 + 3 * colw
    };

    /* 行 0：模式名 */
    TextSetColor(m_G, sectionColor);
    TextDrawStrAt(m_G, "鼠标", x, y, orthoCGO);
    TextSetColor(m_G, accentColor);
    TextDrawStrAt(m_G, zh_mouse_mode(SettingGetGlobal_s(m_G, cSetting_button_mode_name)),
                 x + col0, y, orthoCGO);
    y -= cButModeLineHeight;

    if(SettingGetGlobal_b(m_G, cSetting_mouse_grid)) {
      /* 模式与表格之间一条中性线，不占整行 */
      butmode_draw_hline(m_G, x, xRight, y + cButModeLineHeight - DIP2PIXEL(2),
                         dividerColor, orthoCGO);

      /* 表头：四列独立绘制 */
      TextSetColor(m_G, labelColor);
      TextDrawStrAt(m_G, "左键", col_x[0], y, orthoCGO);
      TextDrawStrAt(m_G, "中键", col_x[1], y, orthoCGO);
      TextDrawStrAt(m_G, "右键", col_x[2], y, orthoCGO);
      TextDrawStrAt(m_G, "滚轮", col_x[3], y, orthoCGO);
      y -= cButModeLineHeight;

      /* 修饰键行：动作一律来自 Mode[]→Code[]，不写死查看模式文案 */
      const char *sc_label = "S+C";
      struct {
        const char *label;
        int idx[4];
      } rows[] = {
        { "普通", { 0, 1, 2, 12 } },
        { "Shift", { 3, 4, 5, 13 } },
        { "Ctrl", { 6, 7, 8, 14 } },
        { sc_label, { 9, 10, 11, 15 } },
        { "单击", { 19, 20, 21, -1 } },
        { "双击", { 16, 17, 18, -1 } },
      };

      for(int r = 0; r < 6; r++) {
        if(y <= y_sele)
          break;
        TextSetColor(m_G, labelColor);
        TextDrawStrAt(m_G, rows[r].label, x, y, orthoCGO);
        butmode_draw_action_row(m_G, I, rows[r].idx, col_x, y, actionColor, orthoCGO);
        y -= cButModeLineHeight;
      }

      /* 表格与选择行之间细线 */
      butmode_draw_hline(m_G, x, xRight, y_sele + cButModeLineHeight - DIP2PIXEL(2),
                         dividerColor, orthoCGO);
    }

    /* 选择 / 拾取（青强调），落在倒数第二行 */
    {
      TextSetColor(m_G, sectionColor);
      if(ButModeTranslate(m_G, P_GLUT_SINGLE_LEFT, 0) == cButModePickAtom) {
        TextDrawStrAt(m_G, "拾取", x, y_sele, orthoCGO);
        TextSetColor(m_G, accentColor);
        TextDrawStrAt(m_G, "原子（及关节点）", x + col0, y_sele, orthoCGO);
      } else {
        TextDrawStrAt(m_G, "选择", x, y_sele, orthoCGO);
        TextSetColor(m_G, accentColor);
        TextDrawStrAt(m_G, zh_sele_mode(SettingGetGlobal_i(m_G, cSetting_mouse_selection_mode)),
                     x + col0, y_sele, orthoCGO);
      }
    }
  }
  if (!orthoCGO || !(SettingGetGlobal_b(m_G, cSetting_show_frame_rate) || MoviePlaying(m_G))) {
    ButModeDrawFastImpl(this, true, orthoCGO);
  }
}

bool CButMode::fastDraw(CGO* orthoCGO){
  return ButModeDrawFastImpl(this, false, orthoCGO);
}

static bool ButModeDrawFastImpl(Block * block, short definitely , CGO *orthoCGO)
{
  PyMOLGlobals *G = block->m_G;
  CButMode *I = block->m_G->ButMode;
  int x, y;
  const float labelColor[3] = { 0.55f, 0.55f, 0.58f };
  const float valueColor[3] = { 0.82f, 0.82f, 0.82f };

  if (!definitely && (!(SettingGetGlobal_b(G, cSetting_show_frame_rate) || MoviePlaying(G)))) {
    return false;
  }

  /* 底栏中文同样需要系统字体，否则「帧/状态/赫兹」会变问号 */
  TextOverlayFontScope zh_font(G);

  x = I->rect.left + cButModeLeftMargin;
  int xRight = I->rect.right - cButModeLeftMargin;
  y = I->rect.bottom + cButModeBottomMargin;

#ifndef PURE_OPENGL_ES_2
    {
        int buffer;
        /* TODO : Why do we only do this for the back right buffer,
         for performance? */
        glGetIntegerv(GL_DRAW_BUFFER, (GLint *) & buffer);
        if(buffer != GL_BACK_RIGHT) {
#else
    {
        {
#endif
            if(I->Delay <= 0.0F) {
                if(I->Samples > 0.0F)
                    I->RateShown = (I->Rate / I->Samples);
                else
                    I->RateShown = 0.0F;
                I->Delay = 0.2F;
            }
        }
    }

  {
    int has_movie = false;
    int frame_rate = SettingGetGlobal_b(G, cSetting_show_frame_rate);
    int nf;
    char rateStr[255];
    int col0 = DIP2PIXEL(36);
    int width = I->rect.right - I->rect.left;
    if(width >= DIP2PIXEL(200))
      col0 = DIP2PIXEL(44);
    else if(width < DIP2PIXEL(176))
      col0 = DIP2PIXEL(30);

    nf = SceneGetNFrame(G, &has_movie);
    if(nf == 0)
      nf = 1;
    TextSetColor(G, labelColor);
    if(has_movie) {
      TextDrawStrAt(G, "帧", x, y, orthoCGO);
    } else {
      TextDrawStrAt(G, "状态", x, y, orthoCGO);
    }
    TextSetColor(G, valueColor);
    snprintf(rateStr, sizeof(rateStr), "%d/%d", SceneGetFrame(G) + 1, nf);
    TextDrawStrAt(G, rateStr, x + col0, y, orthoCGO);

    if(frame_rate) {
      /* 右侧预留「赫兹」，避免单位溢出侧栏 */
      int unit_w = DIP2PIXEL(28);
      int rate_w = DIP2PIXEL(40);
      int unit_x = xRight - unit_w;
      int rate_x = unit_x - rate_w;
      if(rate_x < x + col0 + DIP2PIXEL(48))
        rate_x = x + col0 + DIP2PIXEL(48);
      snprintf(rateStr, sizeof(rateStr), "%5.1f", I->RateShown);
      TextSetColor(G, valueColor);
      TextDrawStrAt(G, rateStr, rate_x, y, orthoCGO);
      TextSetColor(G, labelColor);
      TextDrawStrAt(G, "赫兹", unit_x, y, orthoCGO);
    } else if(has_movie) {
      int st_x = x + col0 + DIP2PIXEL(72);
      if(st_x + DIP2PIXEL(48) > xRight)
        st_x = xRight - DIP2PIXEL(56);
      if(st_x < x + col0 + DIP2PIXEL(40))
        st_x = x + col0 + DIP2PIXEL(40);
      TextSetColor(G, labelColor);
      TextDrawStrAt(G, "状态", st_x, y, orthoCGO);
      TextSetColor(G, valueColor);
      snprintf(rateStr, sizeof(rateStr), "%d", SceneGetState(G) + 1);
      TextDrawStrAt(G, rateStr, st_x + DIP2PIXEL(32), y, orthoCGO);
    }
  }
  return true;
}


/*========================================================================*/
int ButModeInit(PyMOLGlobals * G)
{
  CButMode *I = nullptr;
  if((I = (G->ButMode = new CButMode(G)))) {

    int a;

    I->Rate = 0.0;
    I->Samples = 0.0;
    I->RateShown = 0.0;
    I->Delay = 0.0;
    I->DeferCnt = 0;
    I->DeferTime = 0.0F;
    I->NCode = cButModeCount;
    I->NBut = cButModeInputCount;

    for(a = 0; a < I->NBut; a++) {
      I->Mode[a] = -1;
    }

    strcpy(I->Code[cButModeRotXYZ], "旋转");
    strcpy(I->Code[cButModeRotZ], "转轴");
    strcpy(I->Code[cButModeTransXY], "平移");
    strcpy(I->Code[cButModeTransZ], "远近");
    strcpy(I->Code[cButModeClipNF], "裁剪");
    strcpy(I->Code[cButModeClipN], "近裁");
    strcpy(I->Code[cButModeClipF], "远裁");
    strcpy(I->Code[cButModePickAtom], "点选");
    strcpy(I->Code[cButModePickBond], "选键");
    strcpy(I->Code[cButModeTorFrag], "扭转");
    strcpy(I->Code[cButModeRotFrag], "转段");
    strcpy(I->Code[cButModeMovFrag], "移段");
    strcpy(I->Code[cButModeLB], "左键");
    strcpy(I->Code[cButModeMB], "中键");
    strcpy(I->Code[cButModeRB], "右键");
    strcpy(I->Code[cButModeAddToLB], "加左");
    strcpy(I->Code[cButModeAddToMB], "加中");
    strcpy(I->Code[cButModeAddToRB], "加右");
    strcpy(I->Code[cButModeOrigAt], "原点");
    strcpy(I->Code[cButModeRectAdd], "加框");
    strcpy(I->Code[cButModeRectSub], "减框");
    strcpy(I->Code[cButModeRect], "框选");
    strcpy(I->Code[cButModeNone], "  -  ");
    strcpy(I->Code[cButModeCent], "居中");
    strcpy(I->Code[cButModePkTorBnd], "选扭");
    strcpy(I->Code[cButModeScaleSlab], "厚度");
    strcpy(I->Code[cButModeMoveSlab], "移层");
    strcpy(I->Code[cButModePickAtom1], "选一");
    strcpy(I->Code[cButModeMoveAtom], "移原");
    strcpy(I->Code[cButModeMenu], "菜单");
    strcpy(I->Code[cButModeSeleSet], "选定");
    strcpy(I->Code[cButModeSeleToggle], "加减");
    strcpy(I->Code[cButModeSeleAddBox], "加框");
    strcpy(I->Code[cButModeSeleSubBox], "减框");
    strcpy(I->Code[cButModeMoveSlabAndZoom], "层放");
    strcpy(I->Code[cButModeSimpleClick], "单击");
    strcpy(I->Code[cButModeRotDrag], "拖转");
    strcpy(I->Code[cButModeMovDrag], "拖移");
    strcpy(I->Code[cButModeMovDragZ], "拖远");
    strcpy(I->Code[cButModeRotObj], "转物");
    strcpy(I->Code[cButModeMovObj], "移物");
    strcpy(I->Code[cButModeMovObjZ], "物远");
    strcpy(I->Code[cButModeMovFragZ], "段远");
    strcpy(I->Code[cButModeMoveAtomZ], "原远");
    strcpy(I->Code[cButModeDragMol], "拖分");
    strcpy(I->Code[cButModeRotView], "转视");
    strcpy(I->Code[cButModeMovView], "移视");
    strcpy(I->Code[cButModeMovViewZ], "视远");
    strcpy(I->Code[cButModeDragObj], "拖物");
    strcpy(I->Code[cButModeInvMoveSlabAndZoom], "反层");
    strcpy(I->Code[cButModeInvTransZ], "反远");
    strcpy(I->Code[cButModeSeleSetBox], "方框");
    strcpy(I->Code[cButModeInvRotZ], "反转");
    strcpy(I->Code[cButModeRotL], "转光");
    strcpy(I->Code[cButModeMovL], "移光");
    strcpy(I->Code[cButModeMvzL], "光远");

    I->active = true;

    /* 卡片背景略深于旧版 0.2 灰，分区标签用中性灰 */
    I->BackColor[0] = 0.12F;
    I->BackColor[1] = 0.12F;
    I->BackColor[2] = 0.13F;

    I->TextColor[0] = 0.62F;
    I->TextColor[1] = 0.62F;
    I->TextColor[2] = 0.64F;

    I->TextColor1[0] = 0.55F;
    I->TextColor1[1] = 0.55F;
    I->TextColor1[2] = 0.58F;

    I->TextColor2[0] = 0.88F;
    I->TextColor2[1] = 0.88F;
    I->TextColor2[2] = 0.88F;

    I->TextColor3[0] = 0.45F;
    I->TextColor3[1] = 0.82F;
    I->TextColor3[2] = 0.85F;

    OrthoAttach(G, I, cOrthoTool);
    return 1;
  } else
    return 0;
}


/*========================================================================*/
int ButModeCheckPossibleSingleClick(PyMOLGlobals * G, int button, int mod)
{
  int click_button = -1;
  switch (button) {
  case P_GLUT_LEFT_BUTTON:
    click_button = P_GLUT_SINGLE_LEFT;
    break;
  case P_GLUT_MIDDLE_BUTTON:
    click_button = P_GLUT_SINGLE_MIDDLE;
    break;
  case P_GLUT_RIGHT_BUTTON:
    click_button = P_GLUT_SINGLE_RIGHT;
    break;
  }
  if(click_button < 0)
    return false;
  else
    return (ButModeTranslate(G, click_button, mod) >= 0);
}

int ButModeTranslate(PyMOLGlobals * G, int button, int mod)
{
  int mode = cButModeNothing;
  CButMode *I = G->ButMode;
  switch (button) {
  case P_GLUT_LEFT_BUTTON:
    mode = 0;
    break;
  case P_GLUT_MIDDLE_BUTTON:
    mode = 1;
    break;
  case P_GLUT_RIGHT_BUTTON:
    mode = 2;
    break;
  case P_GLUT_BUTTON_SCROLL_FORWARD:
  case P_GLUT_BUTTON_SCROLL_BACKWARD:
    switch (mod) {
    case 0:
      mode = 12;
      break;
    case cOrthoSHIFT:
      mode = 13;
      break;
    case cOrthoCTRL:
      mode = 14;
      break;
    case (cOrthoCTRL + cOrthoSHIFT):
      mode = 15;
    }
    mod = 0;
    switch (I->Mode[mode]) {
    case cButModeScaleSlab:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeScaleSlabExpand;
      } else {
        return cButModeScaleSlabShrink;
      }
      break;
    case cButModeMoveSlab:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeMoveSlabForward;
      } else {
        return cButModeMoveSlabBackward;
      }
      break;
    case cButModeMoveSlabAndZoom:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeMoveSlabAndZoomForward;
      } else {
        return cButModeMoveSlabAndZoomBackward;
      }
      break;
    case cButModeInvMoveSlabAndZoom:
      if(button != P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeMoveSlabAndZoomForward;
      } else {
        return cButModeMoveSlabAndZoomBackward;
      }
      break;
    case cButModeTransZ:
      if(button == P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeZoomForward;
      } else {
        return cButModeZoomBackward;
      }
      break;
    case cButModeInvTransZ:
      if(button != P_GLUT_BUTTON_SCROLL_FORWARD) {
        return cButModeZoomForward;
      } else {
        return cButModeZoomBackward;
      }
      break;
    }
    return -1;
    break;
  case P_GLUT_DOUBLE_LEFT:
  case P_GLUT_DOUBLE_MIDDLE:
  case P_GLUT_DOUBLE_RIGHT:
  case P_GLUT_SINGLE_LEFT:
  case P_GLUT_SINGLE_MIDDLE:
  case P_GLUT_SINGLE_RIGHT:
    switch (button) {
    case P_GLUT_DOUBLE_LEFT:
      mode = 16;
      break;
    case P_GLUT_DOUBLE_MIDDLE:
      mode = 17;
      break;
    case P_GLUT_DOUBLE_RIGHT:
      mode = 18;
      break;
    case P_GLUT_SINGLE_LEFT:
      mode = 19;
      break;
    case P_GLUT_SINGLE_MIDDLE:
      mode = 20;
      break;
    case P_GLUT_SINGLE_RIGHT:
      mode = 21;
      break;
    }
    switch (mod) {
    case cOrthoSHIFT:
      mode += 6;
      break;
    case cOrthoCTRL:
      mode += 12;
      break;
    case (cOrthoCTRL + cOrthoSHIFT):
      mode += 18;
      break;
    case cOrthoALT:
      mode += 24;
      break;
    case (cOrthoALT + cOrthoSHIFT):
      mode += 30;
      break;
    case (cOrthoALT + cOrthoCTRL):
      mode += 36;
      break;
    case (cOrthoALT + cOrthoCTRL + cOrthoSHIFT):
      mode += 42;
      break;
    }
    mod = 0;
    break;
  }
  switch (mod) {
  case 0:
    break;
  case cOrthoSHIFT:
    mode += 3;
    break;
  case cOrthoCTRL:
    mode += 6;
    break;
  case (cOrthoCTRL + cOrthoSHIFT):
    mode += 9;
    break;
  case cOrthoALT:
    mode += 68;
    break;
  case (cOrthoALT + cOrthoSHIFT):
    mode += 71;
    break;
  case (cOrthoALT + cOrthoCTRL):
    mode += 74;
    break;
  case (cOrthoALT + cOrthoCTRL + cOrthoSHIFT):
    mode += 77;
    break;
  }
  return (I->Mode[mode]);
}
