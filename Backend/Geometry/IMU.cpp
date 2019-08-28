/******************************************************************************
 * Copyright 2017-2018 Baidu Robotic Vision Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/
#include "stdafx.h"
//#ifndef CFG_DEBUG
//#define CFG_DEBUG
//#endif
//#ifdef CFG_DEBUG
//#define CFG_DEBUG_EIGEN
//#endif
#include "IMU.h"

namespace IMU {

void Delta::Factor::Auxiliary::Global::Set(const Jacobian::Global &J, const Error &e,
                                           const float w, const Weight &W, const float Tpv) {
  J.GetTranspose(m_JT);
  W.GetScaled(w, &m_W);
  const xp128f _Tpv = xp128f::get(Tpv);

  for (int i = 0; i < 5; ++i) {
    const LA::AlignedMatrix3x3f *Wi = m_W[i];
    const LA::AlignedMatrix3x3f &Wir = Wi[0];
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrr1, Wir, m_JTW[1][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrbw1, Wir, m_JTW[4][i]);
    m_JTW[1][i].GetMinus(m_JTW[6][i]);
    const LA::AlignedMatrix3x3f &Wiv = Wi[1];
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvr1, Wiv, m_JTW[1][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv1, Wiv, m_JTW[2][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvba1, Wiv, m_JTW[3][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvbw1, Wiv, m_JTW[4][i]);
    m_JTW[2][i].GetMinus(m_JTW[7][i]);
    const LA::AlignedMatrix3x3f &Wip = Wi[2];
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jpp1, Wip, m_JTW[0][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr1, Wip, m_JTW[1][i]);
    //LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpv1, Wip, m_JTW[2][i]);
    LA::AlignedMatrix3x3f::AddsATo(_Tpv, m_JTW[0][i], m_JTW[2][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpba1, Wip, m_JTW[3][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpbw1, Wip, m_JTW[4][i]);
    m_JTW[0][i].GetMinus(m_JTW[5][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr2, Wip, m_JTW[6][i]);
    m_JTW[3][i] += m_W[3][i];
    m_W[3][i].GetMinus(m_JTW[8][i]);
    m_JTW[4][i] += m_W[4][i];
    m_W[4][i].GetMinus(m_JTW[9][i]);
  }

  LA::AlignedMatrix3x3f *A[10] = {&m_A[ 0], &m_A[ 9], &m_A[17], &m_A[24], &m_A[30],
                                  &m_A[35], &m_A[39], &m_A[42], &m_A[44], &m_A[45]};
  const LA::AlignedVector3f *_e = (LA::AlignedVector3f *) &e;
  LA::AlignedMatrix3x3f WJ;
  for (int i = 0; i < 10; ++i) {
    const LA::AlignedMatrix3x3f *JTWi = m_JTW[i];
    if (i >= 1) {
      const LA::AlignedMatrix3x3f &JTWir = JTWi[0];
      LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrr1, JTWir, A[1][i], i == 1);
      if (i >= 4) {
        LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrbw1, JTWir, A[4][i], i == 4);
        if (i >= 6) {
          A[1][i].GetMinus(A[6][i]);
        }
      }
      const LA::AlignedMatrix3x3f &JTWiv = JTWi[1];
      LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvr1, JTWiv, A[1][i], i == 1);
      if (i >= 2) {
        LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv1, JTWiv, A[2][i], i == 2);
        if (i >= 3) {
          LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvba1, JTWiv, A[3][i], i == 3);
          if (i >= 4) {
            LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvbw1, JTWiv, A[4][i], i == 4);
          }
          if (i >= 7) {
            A[2][i].GetMinus(A[7][i]);
          }
        }
      }
    }
    const LA::AlignedMatrix3x3f &JTWip = JTWi[2];
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jpp1, JTWip, A[0][i], i == 0);
    if (i >= 1) {
      LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr1, JTWip, A[1][i], i == 1);
      if (i >= 2) {
        //LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpv1, JTWip, A[2][i], i == 2);
        LA::AlignedMatrix3x3f::AddsATo(_Tpv, A[0][i], A[2][i], i == 2);
        if (i >= 3) {
          LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpba1, JTWip, A[3][i], i == 3);
          if (i >= 4) {
            LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpbw1, JTWip, A[4][i], i == 4);
            if (i >= 5) {
              A[0][i].GetMinus(A[5][i]);
              if (i >= 6) {
                LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr2, JTWip, A[6][i], i == 6);
              }
            }
          }
        }
      }
    }
    if (i >= 3) {
      JTWi[3].GetTranspose(WJ);
      A[3][i] += WJ;
      if (i == 3) {
        A[3][i].SetLowerFromUpper();
      } else if (i >= 8) {
        WJ.GetMinus(A[8][i]);
      }
      if (i >= 4) {
        JTWi[4].GetTranspose(WJ);
        A[4][i] += WJ;
        if (i == 4) {
          A[4][i].SetLowerFromUpper();
        } else if (i >= 9) {
          WJ.GetMinus(A[9][i]);
        }
      }
    }
    LA::AlignedVector3f &bi = m_b[i];
    bi.MakeZero();
    float *_bi = bi;
    for (int j = 0; j < 5; ++j) {
      LA::AlignedMatrix3x3f::AddAbTo(JTWi[j], _e[j], _bi);
    }
  }
}

void Delta::Factor::Auxiliary::RelativeLF::Set(const Jacobian::RelativeLF &J, const Error &e,
                                               const float w, const Weight &W, const float Tpv) {
  J.GetTranspose(m_JT);
  J.m_Jvr2.GetTranspose(m_Jvr2T);
  J.m_Jvv2.GetTranspose(m_Jvv2T);
  W.GetScaled(w, &m_W);
  const xp128f Jpv1 = xp128f::get(-Tpv);

  for (int i = 0; i < 5; ++i) {
    const LA::AlignedMatrix3x3f *Wi = m_W[i];
    const LA::AlignedMatrix3x3f &Wir = Wi[0];
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrr1, Wir, m_JTW[1][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrbw1, Wir, m_JTW[4][i]);
    m_JTW[1][i].GetMinus(m_JTW[6][i]);
    const LA::AlignedMatrix3x3f &Wiv = Wi[1];
    LA::AlignedMatrix2x3f::ABT(J.m_JvgT, Wiv, m_JTWg[i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvr1, Wiv, m_JTW[1][i]);
    //LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv1, Wiv, m_JTW[2][i]);
    m_W[1][i].GetMinus(m_JTW[2][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvba1, Wiv, m_JTW[3][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvbw1, Wiv, m_JTW[4][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_Jvr2T, Wiv, m_JTW[6][i]);
    LA::AlignedMatrix3x3f::ABT(m_Jvv2T, Wiv, m_JTW[7][i]);
    const LA::AlignedMatrix3x3f &Wip = Wi[2];
    LA::AlignedMatrix2x3f::AddABTTo(J.m_JpgT, Wip, m_JTWg[i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jpp1, Wip, m_JTW[0][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr1, Wip, m_JTW[1][i]);
    //LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpv1, Wip, m_JTW[2][i]);
    LA::AlignedMatrix3x3f::AddsATo(Jpv1, m_W[2][i], m_JTW[2][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpba1, Wip, m_JTW[3][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpbw1, Wip, m_JTW[4][i]);
    m_JTW[0][i].GetMinus(m_JTW[5][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr2, Wip, m_JTW[6][i]);
    m_JTW[3][i] += m_W[3][i];
    m_W[3][i].GetMinus(m_JTW[8][i]);
    m_JTW[4][i] += m_W[4][i];
    m_W[4][i].GetMinus(m_JTW[9][i]);
  }

  LA::AlignedMatrix2x3f::ABT(m_JTWg[1], J.m_JvgT, m_Agg);
  LA::AlignedMatrix2x3f::AddABTTo(m_JTWg[2], J.m_JpgT, m_Agg);
  const LA::AlignedVector3f *_e = (LA::AlignedVector3f *) &e;
  m_bg.MakeZero();
  for (int i = 0; i < 5; ++i) {
    LA::AlignedMatrix2x3f::AddAbTo(m_JTWg[i], _e[i], m_bg);
  }
  LA::AlignedMatrix3x3f *A[10] = {&m_A[ 0], &m_A[ 9], &m_A[17], &m_A[24], &m_A[30],
                                  &m_A[35], &m_A[39], &m_A[42], &m_A[44], &m_A[45]};
  LA::AlignedMatrix3x3f WJ;
  for (int i = 0; i < 10; ++i) {
    const LA::AlignedMatrix3x3f *JTWi = m_JTW[i], &JTWiv = JTWi[1], &JTWip = JTWi[2];
    LA::AlignedMatrix2x3f::ABT(J.m_JvgT, JTWiv, m_Agc[i]);
    LA::AlignedMatrix2x3f::AddABTTo(J.m_JpgT, JTWip, m_Agc[i]);
    if (i >= 1) {
      const LA::AlignedMatrix3x3f &JTWir = JTWi[0];
      LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrr1, JTWir, A[1][i], i == 1);
      if (i >= 4) {
        LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrbw1, JTWir, A[4][i], i == 4);
        if (i >= 6) {
          A[1][i].GetMinus(A[6][i]);
        }
      }
      LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvr1, JTWiv, A[1][i], i == 1);
      if (i >= 2) {
        //LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv1, JTWiv, A[2][i], i == 2);
        JTWiv.GetTranspose(A[2][i]);
        A[2][i].MakeMinus();
        if (i == 2) {
          A[2][i].SetLowerFromUpper();
        }
        if (i >= 3) {
          LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvba1, JTWiv, A[3][i], i == 3);
          if (i >= 4) {
            LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvbw1, JTWiv, A[4][i], i == 4);
          }
          if (i >= 6) {
            LA::AlignedMatrix3x3f::AddABTTo(m_Jvr2T, JTWiv, A[6][i], i == 6);
            if (i >= 7) {
              LA::AlignedMatrix3x3f::ABT(m_Jvv2T, JTWiv, A[7][i], i == 7);
            }
          }
        }
      }
    }
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jpp1, JTWip, A[0][i], i == 0);
    if (i >= 1) {
      LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr1, JTWip, A[1][i], i == 1);
      if (i >= 2) {
        //LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpv1, JTWip, A[2][i], i == 2);
        JTWip.GetTranspose(WJ);
        LA::AlignedMatrix3x3f::AddsATo(Jpv1, WJ, A[2][i], i == 2);
        if (i >= 3) {
          LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpba1, JTWip, A[3][i], i == 3);
          if (i >= 4) {
            LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpbw1, JTWip, A[4][i], i == 4);
            if (i >= 5) {
              A[0][i].GetMinus(A[5][i]);
              if (i >= 6) {
                LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr2, JTWip, A[6][i], i == 6);
              }
            }
          }
        }
      }
    }
    if (i >= 3) {
      JTWi[3].GetTranspose(WJ);
      A[3][i] += WJ;
      if (i == 3) {
        A[3][i].SetLowerFromUpper();
      } else if (i >= 8) {
        WJ.GetMinus(A[8][i]);
      }
      if (i >= 4) {
        JTWi[4].GetTranspose(WJ);
        A[4][i] += WJ;
        if (i == 4) {
          A[4][i].SetLowerFromUpper();
        } else if (i >= 9) {
          WJ.GetMinus(A[9][i]);
        }
      }
    }
    LA::AlignedVector3f &bi = m_b[i];
    bi.MakeZero();
    float *_bi = bi;
    for (int j = 0; j < 5; ++j) {
      LA::AlignedMatrix3x3f::AddAbTo(JTWi[j], _e[j], _bi);
    }
  }
}

void Delta::Factor::Auxiliary::RelativeKF::Set(const Jacobian::RelativeKF &J, const Error &e,
                                               const float w, const Weight &W, const float Tpv) {
  J.GetTranspose(m_JT);
  W.GetScaled(w, &m_W);
  const xp128f Jpv1 = xp128f::get(-Tpv);

  for (int i = 0; i < 5; ++i) {
    const LA::AlignedMatrix3x3f *Wi = m_W[i];
    const LA::AlignedMatrix3x3f &Wir = Wi[0];
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrbw1, Wir, m_JTWc[2][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrr2, Wir, m_JTWc[4][i]);
    const LA::AlignedMatrix3x3f &Wiv = Wi[1];
    LA::AlignedMatrix2x3f::ABT(J.m_JvgT, Wiv, m_JTWg[i]);
    //LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv1, Wiv, m_JTWc[0][i]);
    m_W[1][i].GetMinus(m_JTWc[0][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvba1, Wiv, m_JTWc[1][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvbw1, Wiv, m_JTWc[2][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvr2, Wiv, m_JTWc[4][i]);
    LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv2, Wiv, m_JTWc[5][i]);
    const LA::AlignedMatrix3x3f &Wip = Wi[2];
    LA::AlignedMatrix2x3f::AddABTTo(J.m_JpgT, Wip, m_JTWg[i]);
    //LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpv1, Wip, m_JTW[2][i]);
    LA::AlignedMatrix3x3f::AddsATo(Jpv1, m_W[2][i], m_JTWc[0][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpba1, Wip, m_JTWc[1][i]);
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpbw1, Wip, m_JTWc[2][i]);
    //LA::AlignedMatrix3x3f::ABT(m_JT.m_Jpp2, Wip, m_JTWc[3][i]);
    m_JTWc[3][i] = m_W[2][i];
    LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr2, Wip, m_JTWc[4][i]);
    m_JTWc[1][i] += m_W[3][i];
    m_W[3][i].GetMinus(m_JTWc[6][i]);
    m_JTWc[2][i] += m_W[4][i];
    m_W[4][i].GetMinus(m_JTWc[7][i]);
  }

  LA::AlignedMatrix2x3f::ABT(m_JTWg[1], J.m_JvgT, m_Agg);
  LA::AlignedMatrix2x3f::AddABTTo(m_JTWg[2], J.m_JpgT, m_Agg);
  const LA::AlignedVector3f *_e = (LA::AlignedVector3f *) &e;
  m_bg.MakeZero();
  for (int i = 0; i < 5; ++i) {
    LA::AlignedMatrix2x3f::AddAbTo(m_JTWg[i], _e[i], m_bg);
  }
  LA::AlignedMatrix3x3f *Ac[8] = {&m_Ac[ 0], &m_Ac[ 7], &m_Ac[13], &m_Ac[18],
                                  &m_Ac[22], &m_Ac[25], &m_Ac[27], &m_Ac[28]};
  LA::AlignedMatrix3x3f WJ;
  for (int i = 0; i < 8; ++i) {
    const LA::AlignedMatrix3x3f *JTWi = m_JTWc[i], &JTWiv = JTWi[1], &JTWip = JTWi[2];
    LA::AlignedMatrix2x3f::ABT(J.m_JvgT, JTWiv, m_Agc[i]);
    LA::AlignedMatrix2x3f::AddABTTo(J.m_JpgT, JTWip, m_Agc[i]);
    const LA::AlignedMatrix3x3f &JTWir = JTWi[0];
    if (i >= 2) {
      LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrbw1, JTWir, Ac[2][i], i == 2);
      if (i >= 4) {
        LA::AlignedMatrix3x3f::ABT(m_JT.m_Jrr2, JTWir, Ac[4][i], i == 4);
      }
    }
    //LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv1, JTWiv, Ac[0][i], i == 0);
    JTWiv.GetTranspose(Ac[0][i]);
    Ac[0][i].MakeMinus();
    if (i == 0) {
      Ac[0][i].SetLowerFromUpper();
    }
    if (i >= 1) {
      LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvba1, JTWiv, Ac[1][i], i == 1);
      if (i >= 2) {
        LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvbw1, JTWiv, Ac[2][i], i == 2);
      }
      if (i >= 4) {
        LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jvr2, JTWiv, Ac[4][i], i == 4);
        if (i >= 5) {
          LA::AlignedMatrix3x3f::ABT(m_JT.m_Jvv2, JTWiv, Ac[5][i], i == 5);
        }
      }
    }
    //LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpv1, JTWip, Ac[0][i], i == 0);
    JTWip.GetTranspose(WJ);
    LA::AlignedMatrix3x3f::AddsATo(Jpv1, WJ, Ac[0][i], i == 0);
    if (i >= 1) {
      LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpba1, JTWip, Ac[1][i], i == 1);
      if (i >= 2) {
        LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpbw1, JTWip, Ac[2][i], i == 2);
        if (i >= 3) {
          //LA::AlignedMatrix3x3f::ABT(m_JT.m_Jpp2, JTWip, Ac[3][i], i == 3);
          Ac[3][i] = WJ;
          if (i == 3) {
            Ac[3][i].SetLowerFromUpper();
          }
          if (i >= 4) {
            LA::AlignedMatrix3x3f::AddABTTo(m_JT.m_Jpr2, JTWip, Ac[4][i], i == 4);
          }
        }
      }
      JTWi[3].GetTranspose(WJ);
      Ac[1][i] += WJ;
      if (i == 1) {
        Ac[1][i].SetLowerFromUpper();
      } else if (i >= 6) {
        WJ.GetMinus(Ac[6][i]);
      }
      if (i >= 2) {
        JTWi[4].GetTranspose(WJ);
        Ac[2][i] += WJ;
        if (i == 2) {
          Ac[2][i].SetLowerFromUpper();
        } else if (i >= 7) {
          WJ.GetMinus(Ac[7][i]);
        }
      }
    }
    LA::AlignedVector3f &bi = m_bc[i];
    bi.MakeZero();
    float *_bi = bi;
    for (int j = 0; j < 5; ++j) {
      LA::AlignedMatrix3x3f::AddAbTo(JTWi[j], _e[j], _bi);
    }
  }
}

void InitializeCamera(const AlignedVector<Measurement> &us, Camera &C) {
  const int N = us.Size();
  if (N == 0 || IMU_GRAVITY_EXCLUDED) {
    C.MakeIdentity();
  } else {
    LA::AlignedVector3f g;
    g.MakeZero();
    for (int i = 0; i < N; ++i) {
      g += us[i].m_a;
    }
    g.Normalize();
    g.MakeMinus();
    C.MakeIdentity(&g);
  }
}

void PreIntegrate(const AlignedVector<Measurement> &us, const float t1, const float t2,
                  const Camera &C1, Delta *D, AlignedVector<float> *work, const bool jac,
                  const Measurement *u1, const Measurement *u2, const float eps) {
  PreIntegrate(us.Data(), us.Size(), t1, t2, C1, D, work, jac, u1, u2, eps);
}

void PreIntegrate(const Measurement *us, const int N, const float t1, const float t2,
                  const Camera &C1, Delta *D, AlignedVector<float> *work, const bool jac,
                  const Measurement *u1, const Measurement *u2, const float eps) {
  if (u1) {
    D->m_u1 = *u1;
  } else {
    D->m_u1.Invalidate();
  }
  if (u2) {
    D->m_u2 = *u2;
  } else {
    D->m_u2.Invalidate();
  }
  D->m_ba = C1.m_ba;
  D->m_bw = C1.m_bw;

  D->m_v.MakeZero();
  D->m_p.MakeZero();
  D->m_Jvba.MakeZero();
  D->m_Jvbw.MakeZero();
  D->m_Jpba.MakeZero();
  D->m_Jpbw.MakeZero();

  int r1 = 0, r2 = 1;
  Quaternion dq, q;
  Rotation3D dR, RT[2];
  LA::AlignedMatrix3x3f Jrbw[2];
  q.MakeIdentity();
  RT[r1].MakeIdentity();
  RT[r2].MakeIdentity();
  Jrbw[r1].MakeZero();
  Jrbw[r2].MakeZero();

  Delta::Transition F;
  Delta::Covariance P[2];
  Delta::Covariance::DD U;
  P[r1].MakeZero();
  P[r2].MakeZero();

  float Tpg = 0.0f;
  float _t1, _t2;
  LA::AlignedVector3f a, adt, w, wdt;
  LA::AlignedMatrix3x3f Jr[2], Jrdt;
  const xp128f s = xp128f::get(0.5f);
  for (int i1 = -1, i2 = 0; i1 < N; i1 = i2++) {
    if (i1 != -1) {
      const Measurement &_u1 = us[i1];
      _t1 = _u1.t();
      if (i2 != N) {
        const Measurement &_u2 = us[i2];
        _t2 = _u2.t();
        a.xyzr() = (_u1.m_a.xyzr() + _u2.m_a.xyzr()) * s;
        w.xyzr() = (_u1.m_w.xyzr() + _u2.m_w.xyzr()) * s;
      } else {
        _t2 = t2;
        if (u2) {
          Measurement::Interpolate(_u1, *u2, (_t1 + _t2) * 0.5f, a, w);
        } else {
          a = _u1.m_a;
          w = _u1.m_w;
        }
      }
    } else {
      const Measurement &_u2 = us[i2];
      _t1 = t1;
      _t2 = _u2.t();
      if (u1) {
        Measurement::Interpolate(*u1, _u2, (_t1 + _t2) * 0.5f, a, w);
      } else {
        a = _u2.m_a;
        w = _u2.m_w;
      }
    }
    if (_t1 == _t2) {
      continue;
    }
    a -= D->m_ba;
    w -= D->m_bw;
    const xp128f dt = xp128f::get(_t2 - _t1);
    const xp128f dT = xp128f::get(_t1 - t1);
    const xp128f dt_2 = xp128f::get(dt[0] * 0.5f);
    Tpg = (dt_2[0] + (_t1 - t1)) * dt[0] + Tpg;

    w.GetScaled(dt, wdt);
    //dR.SetRodrigues(wdt);
    //Rotation3D::ABT(dR, RT[r1], D->m_R);
    dq.SetRodrigues(wdt, eps);
    dR.SetQuaternion(dq);
    q = dq * q;
    RT[r2].SetQuaternion(q);
    RT[r2].Transpose();

    const LA::AlignedMatrix3x3f RTdt = (RT[r1] + RT[r2]) * dt_2;
    const LA::AlignedVector3f dv = RTdt * a;
    const LA::AlignedVector3f dp = D->m_v * dt + dv * dt_2;
    D->m_v += dv;
    D->m_p += dp;

    if (jac) {
      Rotation3D::GetRodriguesJacobian(wdt, Jrdt, eps);
      Jrdt *= dt;
      Jrbw[r2] = dR * Jrbw[r1] - Jrdt;
      const LA::AlignedMatrix3x3f dJvba = RTdt.GetMinus();
      const LA::AlignedMatrix3x3f dJpba = D->m_Jvba * dt + dJvba * dt_2;
      D->m_Jvba += dJvba;
      D->m_Jpba += dJpba;
      a.GetScaled(dt_2, adt);
      SkewSymmetricMatrix::ABT(RT[r1], adt, Jr[r1]);
      SkewSymmetricMatrix::ABT(RT[r2], adt, Jr[r2]);
      const LA::AlignedMatrix3x3f dJvbw = (Jr[r1] * Jrbw[r1]) + (Jr[r2] * Jrbw[r2]);
      const LA::AlignedMatrix3x3f dJpbw = D->m_Jvbw * dt + dJvbw * dt_2;
      D->m_Jvbw += dJvbw;
      D->m_Jpbw += dJpbw;

      if (i1 != -1) {
        RT[r2].GetScaled(dt, F.m_Fdb.m_Frbw);
        F.m_Fdb.m_Frbw.MakeMinus();
        dv.GetMinus(F.m_Fdd.m_Fvr);
        F.m_Fdb.m_Fvba = dJvba;
        F.m_Fdb.m_Fvbw = dJvbw;
        dp.GetMinus(F.m_Fdd.m_Fpr);
        F.m_Fdd.m_Fpv = dt;
        F.m_Fdb.m_Fpba = dJpba;
        F.m_Fdb.m_Fpbw = dJpbw;
        Delta::Covariance::FPFT(F, P[r1], &U, &P[r2]);
      }
      const float s2r = dt[0] * IMU_VARIANCE_GYROSCOPE_NOISE;
      const float s2v = dt[0] * IMU_VARIANCE_ACCELERATION_NOISE;
      const float s2p = dt[0] * dt_2[0] * s2v;
      const float s2ba = dt[0] * IMU_VARIANCE_ACCELERATION_BIAS_WALK;
      const float s2bw = dt[0] * IMU_VARIANCE_GYROSCOPE_BIAS_WALK;
      P[r2].IncreaseDiagonal(s2r, s2v, s2p, s2ba, s2bw);
    }

    r1 = r2;
    r2 = 1 - r2;
  }
  D->m_RT = RT[r1];
  D->m_Jrbw = Jrbw[r1];
  D->m_Tvg = D->m_Tpv = t2 - t1;
  D->m_Tpg = Tpg;
  //D->m_Tpg = 0.5f * D.m_Tvg * D.m_Tvg;
  if (jac) {
    P[r1].IncreaseDiagonal(IMU_VARIANCE_EPSILON_ROTATION,
                           IMU_VARIANCE_EPSILON_VELOCITY,
                           IMU_VARIANCE_EPSILON_POSITION,
                           IMU_VARIANCE_EPSILON_BIAS_ACCELERATION,
                           IMU_VARIANCE_EPSILON_BIAS_GYROSCOPE);
    D->m_W.Set(P[r1], work);
    const float w[5] = {IMU_WEIGHT_ROTATION, IMU_WEIGHT_VELOCITY, IMU_WEIGHT_POSITION,
                        IMU_WEIGHT_BIAS_ACCELERATION, IMU_WEIGHT_BIAS_GYROSCOPE};
    for (int i = 0; i < 5; ++i) {
      D->m_W[i][i] *= w[i];
      for (int j = i + 1; j < 5; ++j) {
        D->m_W[i][j] *= sqrtf(w[i] * w[j]);
        D->m_W[i][j].GetTranspose(D->m_W[j][i]);
      }
    }
  }
}

void Propagate(const Point3D &pu, const Delta &D, const Camera &C1, Camera &C2, const float eps) {
  const LA::AlignedVector3f dba = C1.m_ba - D.m_ba;
  const LA::AlignedVector3f dbw = C1.m_bw - D.m_bw;
  const Rotation3D R1T = C1.m_T.GetRotationTranspose();
  C2.m_T = D.GetRotationMeasurement(dbw, eps).GetTranspose() / R1T;
  const LA::AlignedVector3f dp = D.m_p + D.m_Jpba * dba + D.m_Jpbw * dbw;
  C2.m_p = (R1T - C2.m_T.GetRotationTranspose()) * pu + C1.m_p +
            C1.m_v * D.m_Tpv + R1T.GetApplied(dp);
  if (!IMU_GRAVITY_EXCLUDED) {
    C2.m_p.z() -= IMU_GRAVITY_MAGNITUDE * D.m_Tpg;
  }
  C2.m_T.SetPosition(C2.m_p);
  const LA::AlignedVector3f dv = D.m_v + D.m_Jvba * dba + D.m_Jvbw * dbw;
  C2.m_v = C1.m_v + R1T.GetApplied(dv);
  if (!IMU_GRAVITY_EXCLUDED) {
    C2.m_v.z() -= IMU_GRAVITY_MAGNITUDE * D.m_Tvg;
  }
  //////////////////////////////////////////////////////////////////////////
  //C2.m_p = C1.m_p;
  //C2.m_T.SetPosition(C2.m_p);
  //C2.m_v.MakeZero();
  //////////////////////////////////////////////////////////////////////////
  C2.m_ba = C1.m_ba;
  C2.m_bw = C1.m_bw;
}

}
