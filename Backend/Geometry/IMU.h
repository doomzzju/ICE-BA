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
#ifndef _IMU_H_
#define _IMU_H_

#include "Camera.h"
#include "Parameter.h"
#include "AlignedVector.h"

namespace IMU {
class Measurement {
 public:
  inline const float& t() const { return m_w.r(); }
  inline       float& t()       { return m_w.r(); }
  inline bool operator < (const float t) const { return this->t() < t; }
  inline bool Valid() const { return m_a.Valid(); }
  inline bool Invalid() const { return m_a.Invalid(); }
  inline void Invalidate() { m_a.Invalidate(); }
  static inline void Interpolate(const Measurement &u1, const Measurement &u2, const float t,
                                 LA::AlignedVector3f &a, LA::AlignedVector3f &w) {
    const xp128f w1 = xp128f::get((u2.t() - t) / (u2.t() - u1.t()));
    const xp128f w2 = xp128f::get(1.0f - w1[0]);
    a.xyzr() = (u1.m_a.xyzr() * w1) + (u2.m_a.xyzr() * w2);
    w.xyzr() = (u1.m_w.xyzr() * w1) + (u2.m_w.xyzr() * w2);
  }
  inline void Print(const bool e = false) const {
    if (e) {
      UT::Print("%e  %e %e %e  %e %e %e\n", t(), m_a.x(), m_a.y(), m_a.z(),
                                                 m_w.x(), m_w.y(), m_w.z());
    } else {
      UT::Print("%f  %f %f %f  %f %f %f\n", t(), m_a.x(), m_a.y(), m_a.z(),
                                                 m_w.x(), m_w.y(), m_w.z());
    }
  }
 public:
  LA::AlignedVector3f m_a, m_w;
};

class Delta {
 public:
  class Transition {
   public:
    class DD {
     public:
      //LA::AlignedMatrix3x3f m_Fvr, m_Fpr, m_Fpv;
      SkewSymmetricMatrix m_Fvr, m_Fpr;
      xp128f m_Fpv;
    };
    class DB {
     public:
      LA::AlignedMatrix3x3f m_Frbw, m_Fvba, m_Fvbw, m_Fpba, m_Fpbw;
    };
   public:
    DD m_Fdd;
    DB m_Fdb;
  };
  class Covariance {
   public:
    class DD {
     public:
      inline void GetUpper(float **P) const {
        memcpy(&P[0][0], &m_Prr.m00(), 12);
        memcpy(&P[0][3], &m_Prv.m00(), 12);
        memcpy(&P[0][6], &m_Prp.m00(), 12);
        memcpy(&P[1][1], &m_Prr.m11(), 8);
        memcpy(&P[1][3], &m_Prv.m10(), 12);
        memcpy(&P[1][6], &m_Prp.m10(), 12);
        P[2][2] = m_Prr.m22();
        memcpy(&P[2][3], &m_Prv.m20(), 12);
        memcpy(&P[2][6], &m_Prp.m20(), 12);
        memcpy(&P[3][3], &m_Pvv.m00(), 12);
        memcpy(&P[3][6], &m_Pvp.m00(), 12);
        memcpy(&P[4][4], &m_Pvv.m11(), 8);
        memcpy(&P[4][6], &m_Pvp.m10(), 12);
        P[5][5] = m_Pvv.m22();
        memcpy(&P[5][6], &m_Pvp.m20(), 12);
        memcpy(&P[6][6], &m_Ppp.m00(), 12);
        memcpy(&P[7][7], &m_Ppp.m11(), 8);
        P[8][8] = m_Ppp.m22();
      }
      inline void IncreaseDiagonal(const float s2r, const float s2v, const float s2p) {
        m_Prr.IncreaseDiagonal(s2r);
        m_Pvv.IncreaseDiagonal(s2v);
        m_Ppp.IncreaseDiagonal(s2p);
      }
      inline void SetLowerFromUpper() {
        m_Prr.SetLowerFromUpper();
        m_Prv.GetTranspose(m_Pvr);
        m_Prp.GetTranspose(m_Ppr);
        m_Pvv.SetLowerFromUpper();
        m_Pvp.GetTranspose(m_Ppv);
        m_Ppp.SetLowerFromUpper();
      }
     public:
      LA::AlignedMatrix3x3f m_Prr, m_Prv, m_Prp;
      LA::AlignedMatrix3x3f m_Pvr, m_Pvv, m_Pvp;
      LA::AlignedMatrix3x3f m_Ppr, m_Ppv, m_Ppp;
    };
    class BD {
     public:
      LA::AlignedMatrix3x3f m_Pbav, m_Pbap;
      LA::AlignedMatrix3x3f m_Pbwr, m_Pbwv, m_Pbwp;
    };
    class DB {
     public:
      inline void Get(float **P, const int j) const {
        const int jbw = j + 3;
        memset(P[0] + j, 0, 12);              memcpy(P[0] + jbw, &m_Prbw.m00(), 12);
        memset(P[1] + j, 0, 12);              memcpy(P[1] + jbw, &m_Prbw.m10(), 12);
        memset(P[2] + j, 0, 12);              memcpy(P[2] + jbw, &m_Prbw.m20(), 12);
        memcpy(P[3] + j, &m_Pvba.m00(), 12);  memcpy(P[3] + jbw, &m_Pvbw.m00(), 12);
        memcpy(P[4] + j, &m_Pvba.m10(), 12);  memcpy(P[4] + jbw, &m_Pvbw.m10(), 12);
        memcpy(P[5] + j, &m_Pvba.m20(), 12);  memcpy(P[5] + jbw, &m_Pvbw.m20(), 12);
        memcpy(P[6] + j, &m_Ppba.m00(), 12);  memcpy(P[6] + jbw, &m_Ppbw.m00(), 12);
        memcpy(P[7] + j, &m_Ppba.m10(), 12);  memcpy(P[7] + jbw, &m_Ppbw.m10(), 12);
        memcpy(P[8] + j, &m_Ppba.m20(), 12);  memcpy(P[8] + jbw, &m_Ppbw.m20(), 12);
      }
      inline void GetTranspose(BD *P) const {
        m_Prbw.GetTranspose(P->m_Pbwr);
        m_Pvba.GetTranspose(P->m_Pbav);
        m_Pvbw.GetTranspose(P->m_Pbwv);
        m_Ppba.GetTranspose(P->m_Pbap);
        m_Ppbw.GetTranspose(P->m_Pbwp);
      }
     public:
      LA::AlignedMatrix3x3f m_Prbw;
      LA::AlignedMatrix3x3f m_Pvba, m_Pvbw;
      LA::AlignedMatrix3x3f m_Ppba, m_Ppbw;
    };
    class BB {
     public:
      inline void GetUpper(float **P, const int i, const int j) const {
        P[i][j] = P[i + 1][j + 1] = P[i + 2][j + 2] = m_Pbaba;
        P[i + 3][j + 3] = P[i + 4][j + 4] = P[i + 5][j + 5] = m_Pbwbw;
        memset(P[i] + j + 1, 0, 20);
        memset(P[i + 1] + j + 2, 0, 16);
        memset(P[i + 2] + j + 3, 0, 12);
        memset(P[i + 3] + j + 4, 0, 8);
        P[i + 4][j + 5] = 0.0f;
      }
      inline void IncreaseDiagonal(const float s2ba, const float s2bw) {
        m_Pbaba += s2ba;
        m_Pbwbw += s2bw;
      }
     public:
      float m_Pbaba, m_Pbwbw;
    };
   public:
    inline void MakeZero() { memset(this, 0, sizeof(Covariance)); }
    inline void IncreaseDiagonal(const float s2r, const float s2v, const float s2p,
                                 const float s2ba, const float s2bw) {
      m_Pdd.IncreaseDiagonal(s2r, s2v, s2p);
      m_Pbb.IncreaseDiagonal(s2ba, s2bw);
    }
    inline void SetLowerFromUpper() {
      m_Pdd.SetLowerFromUpper();
      //m_Pdb.GetTranspose(&m_Pbd);
    }
   public:
    static inline void ABT(const Transition::DD &A, const DD &B, DD *ABT) {
      ABT->m_Prr = B.m_Prr;
      ABT->m_Prv = B.m_Prv;
      ABT->m_Prp = B.m_Prp;
      //LA::AlignedMatrix3x3f::ABT(A.m_Fvr, B.m_Prr, ABT->m_Pvr);
      SkewSymmetricMatrix::AB(A.m_Fvr, B.m_Prr, ABT->m_Pvr);
      ABT->m_Pvr += B.m_Pvr;
      //LA::AlignedMatrix3x3f::ABT(A.m_Fvr, B.m_Pvr, ABT->m_Pvv);
      SkewSymmetricMatrix::AB(A.m_Fvr, B.m_Prv, ABT->m_Pvv);
      ABT->m_Pvv += B.m_Pvv;
      //LA::AlignedMatrix3x3f::ABT(A.m_Fvr, B.m_Ppr, ABT->m_Pvp);
      SkewSymmetricMatrix::AB(A.m_Fvr, B.m_Prp, ABT->m_Pvp);
      ABT->m_Pvp += B.m_Pvp;
      //LA::AlignedMatrix3x3f::ABT(A.m_Fpr, B.m_Prr, ABT->m_Ppr);
      SkewSymmetricMatrix::AB(A.m_Fpr, B.m_Prr, ABT->m_Ppr);
      //LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpv, B.m_Prv, ABT->m_Ppr);
      LA::AlignedMatrix3x3f::AddsATo(A.m_Fpv, B.m_Pvr, ABT->m_Ppr);
      ABT->m_Ppr += B.m_Ppr;
      //LA::AlignedMatrix3x3f::ABT(A.m_Fpr, B.m_Pvr, ABT->m_Ppv);
      SkewSymmetricMatrix::AB(A.m_Fpr, B.m_Prv, ABT->m_Ppv);
      //LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpv, B.m_Pvv, ABT->m_Ppv);
      LA::AlignedMatrix3x3f::AddsATo(A.m_Fpv, B.m_Pvv, ABT->m_Ppv);
      ABT->m_Ppv += B.m_Ppv;
      //LA::AlignedMatrix3x3f::ABT(A.m_Fpr, B.m_Ppr, ABT->m_Ppp);
      SkewSymmetricMatrix::AB(A.m_Fpr, B.m_Prp, ABT->m_Ppp);
      //LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpv, B.m_Ppv, ABT->m_Ppp);
      LA::AlignedMatrix3x3f::AddsATo(A.m_Fpv, B.m_Pvp, ABT->m_Ppp);
      ABT->m_Ppp += B.m_Ppp;
    }
    //static inline void ABT(const Transition::DD &A, const BD &B, DB *ABT) {
    //  B.m_Pbwr.GetTranspose(ABT->m_Prbw);
    //  B.m_Pbav.GetTranspose(ABT->m_Pvba);
    //  B.m_Pbwv.GetTranspose(ABT->m_Pvbw);
    //  LA::AlignedMatrix3x3f::AddABTTo(A.m_Fvr, B.m_Pbwr, ABT->m_Pvbw);
    //  B.m_Pbap.GetTranspose(ABT->m_Ppba);
    //  LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpv, B.m_Pbav, ABT->m_Ppba);
    //  B.m_Pbwp.GetTranspose(ABT->m_Ppbw);
    //  LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpr, B.m_Pbwr, ABT->m_Ppbw);
    //  LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpv, B.m_Pbwv, ABT->m_Ppbw);
    //}
    static inline void AB(const Transition::DD &A, const DB &B, DB *AB) {
      AB->m_Prbw = B.m_Prbw;
      AB->m_Pvba = B.m_Pvba;
      AB->m_Pvbw = B.m_Pvbw;
      SkewSymmetricMatrix::AddABTo(A.m_Fvr, B.m_Prbw, AB->m_Pvbw);
      AB->m_Ppba = B.m_Ppba;
      LA::AlignedMatrix3x3f::AddsATo(A.m_Fpv, B.m_Pvba, AB->m_Ppba);
      AB->m_Ppbw = B.m_Ppbw;
      SkewSymmetricMatrix::AddABTo(A.m_Fpr, B.m_Prbw, AB->m_Ppbw);
      LA::AlignedMatrix3x3f::AddsATo(A.m_Fpv, B.m_Pvbw, AB->m_Ppbw);
    }
    static inline void AddABTTo(const Transition::DB &A, const DB &B, DD *ABT) {
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Frbw, B.m_Prbw, ABT->m_Prr);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Frbw, B.m_Pvbw, ABT->m_Prv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Frbw, B.m_Ppbw, ABT->m_Prp);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fvbw, B.m_Prbw, ABT->m_Pvr);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fvba, B.m_Pvba, ABT->m_Pvv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fvbw, B.m_Pvbw, ABT->m_Pvv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fvba, B.m_Ppba, ABT->m_Pvp);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fvbw, B.m_Ppbw, ABT->m_Pvp);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpbw, B.m_Prbw, ABT->m_Ppr);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpba, B.m_Pvba, ABT->m_Ppv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpbw, B.m_Pvbw, ABT->m_Ppv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpba, B.m_Ppba, ABT->m_Ppp);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Fpbw, B.m_Ppbw, ABT->m_Ppp);
    }
    static inline void AddABTTo(const Transition::DB &A, const BB &B, DB *ABT) {
      const xp128f Bbaba = xp128f::get(B.m_Pbaba);
      const xp128f Bbwbw = xp128f::get(B.m_Pbwbw);
      LA::AlignedMatrix3x3f::AddsATo(Bbwbw, A.m_Frbw, ABT->m_Prbw);
      LA::AlignedMatrix3x3f::AddsATo(Bbaba, A.m_Fvba, ABT->m_Pvba);
      LA::AlignedMatrix3x3f::AddsATo(Bbwbw, A.m_Fvbw, ABT->m_Pvbw);
      LA::AlignedMatrix3x3f::AddsATo(Bbaba, A.m_Fpba, ABT->m_Ppba);
      LA::AlignedMatrix3x3f::AddsATo(Bbwbw, A.m_Fpbw, ABT->m_Ppbw);
    }
    static inline void ABT(const Transition &A, const Covariance &B, DD *ABTdd, DB *ABTdb) {
      ABT(A.m_Fdd, B.m_Pdd, ABTdd);
      //ABT(A.m_Fdd, B.m_Pbd, ABTdb);
      AB(A.m_Fdd, B.m_Pdb, ABTdb);
      AddABTTo(A.m_Fdb, B.m_Pdb, ABTdd);
      AddABTTo(A.m_Fdb, B.m_Pbb, ABTdb);
    }
    static inline void ABTToUpper(const DD &A, const Transition::DD &B, DD *ABT) {
      ABT->m_Prr = A.m_Prr;
      //LA::AlignedMatrix3x3f::ABT(A.m_Prr, B.m_Fvr, ABT->m_Prv);
      SkewSymmetricMatrix::ABT(A.m_Prr, B.m_Fvr, ABT->m_Prv);
      ABT->m_Prv += A.m_Prv;
      //LA::AlignedMatrix3x3f::ABT(A.m_Prr, B.m_Fpr, ABT->m_Prp);
      SkewSymmetricMatrix::ABT(A.m_Prr, B.m_Fpr, ABT->m_Prp);
      //LA::AlignedMatrix3x3f::AddABTTo(A.m_Prv, B.m_Fpv, ABT->m_Prp);
      LA::AlignedMatrix3x3f::AddsATo(B.m_Fpv, A.m_Prv, ABT->m_Prp);
      ABT->m_Prp += A.m_Prp;
      ABT->m_Pvv = A.m_Pvv;
      //LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Pvr, B.m_Fvr, ABT->m_Pvv);
      SkewSymmetricMatrix::AddABTToUpper(A.m_Pvr, B.m_Fvr, ABT->m_Pvv);
      //LA::AlignedMatrix3x3f::ABT(A.m_Pvr, B.m_Fpr, ABT->m_Pvp);
      SkewSymmetricMatrix::ABT(A.m_Pvr, B.m_Fpr, ABT->m_Pvp);
      //LA::AlignedMatrix3x3f::AddABTTo(A.m_Pvv, B.m_Fpv, ABT->m_Pvp);
      LA::AlignedMatrix3x3f::AddsATo(B.m_Fpv, A.m_Pvv, ABT->m_Pvp);
      ABT->m_Pvp += A.m_Pvp;
      ABT->m_Ppp = A.m_Ppp;
      //LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Ppr, B.m_Fpr, ABT->m_Ppp);
      SkewSymmetricMatrix::AddABTToUpper(A.m_Ppr, B.m_Fpr, ABT->m_Ppp);
      //LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Ppv, B.m_Fpv, ABT->m_Ppp);
      LA::AlignedMatrix3x3f::AddsAToUpper(B.m_Fpv, A.m_Ppv, ABT->m_Ppp);
    }
    static inline void AddABTToUpper(const DB &A, const Transition::DB &B, DD *ABT) {
      LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Prbw, B.m_Frbw, ABT->m_Prr);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Prbw, B.m_Fvbw, ABT->m_Prv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Prbw, B.m_Fpbw, ABT->m_Prp);
      LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Pvba, B.m_Fvba, ABT->m_Pvv);
      LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Pvbw, B.m_Fvbw, ABT->m_Pvv);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Pvba, B.m_Fpba, ABT->m_Pvp);
      LA::AlignedMatrix3x3f::AddABTTo(A.m_Pvbw, B.m_Fpbw, ABT->m_Pvp);
      LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Ppba, B.m_Fpba, ABT->m_Ppp);
      LA::AlignedMatrix3x3f::AddABTToUpper(A.m_Ppbw, B.m_Fpbw, ABT->m_Ppp);
    }
    static inline void ABTToUpper(const DD &Add, const DB &Adb, const Transition &B, DD *ABT) {
      ABTToUpper(Add, B.m_Fdd, ABT);
      AddABTToUpper(Adb, B.m_Fdb, ABT);
    }
    static inline void FPFT(const Transition &F, const Covariance &P, DD *U, Covariance *FPFT) {
      ABT(F, P, U, &FPFT->m_Pdb);
      ABTToUpper(*U, FPFT->m_Pdb, F, &FPFT->m_Pdd);
      FPFT->m_Pbb = P.m_Pbb;
      FPFT->SetLowerFromUpper();
    }
   public:
    DD m_Pdd;
    DB m_Pdb;
    //BD m_Pbd;
    BB m_Pbb;
  };

  class Weight {
   public:
    inline const LA::AlignedMatrix3x3f* operator [] (const int i) const { return m_W[i]; }
    inline       LA::AlignedMatrix3x3f* operator [] (const int i)       { return m_W[i]; }
    inline void Set(const Covariance &P, AlignedVector<float> *work) {
      work->Resize(15 * 15);
      float *_P[15];
      _P[0] = work->Data();
      for (int i = 1; i < 15; ++i) {
        _P[i] = _P[i - 1] + 15;
      }
      P.m_Pdd.GetUpper(_P);
      P.m_Pdb.Get(_P, 9);
      P.m_Pbb.GetUpper(_P, 9, 9);
      if (LA::LS::InverseLDL<float>(15, _P)) {
        for (int i = 0, _i = 0; i < 5; ++i) {
          float *W0 = _P[_i++], *W1 = _P[_i++], *W2 = _P[_i++];
          for (int j = 0, _j = 0; j < 5; ++j, _j += 3) {
            m_W[i][j].Set(W0 + _j, W1 + _j, W2 + _j);
          }
        }
      } else {
        MakeZero();
      }
    }
    inline void GetScaled(const float w, Weight *W) const {
      const xp128f _w = xp128f::get(w);
      GetScaled(_w, W);
    }
    inline void GetScaled(const xp128f &w, Weight *W) const {
      for (int i = 0; i < 5; ++i) {
        m_W[i][i].GetScaledToUpper(w, W->m_W[i][i]);
        W->m_W[i][i].SetLowerFromUpper();
        for (int j = i + 1; j < 5; ++j) {
          m_W[i][j].GetScaled(w, W->m_W[i][j]);
          W->m_W[i][j].GetTranspose(W->m_W[j][i]);
        }
      }
    }
    inline void MakeZero() { memset(this, 0, sizeof(Weight)); }
    inline void SetLowerFromUpper() {
      for (int i = 0; i < 5; ++i) {
        m_W[i][i].SetLowerFromUpper();
        for (int j = i + 1; j < 5; ++j) {
          m_W[i][j].GetTranspose(m_W[j][i]);
        }
      }
    }
    inline bool AssertEqual(const Weight &W, const int verbose = 1,
                            const std::string str = "", const bool norm = true) const {
      bool scc = true;
      LA::AlignedMatrix3x3f W1, W2;
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
          W1 = m_W[i][j];
          W2 = W[i][j];
          if (norm) {
            const float n1 = sqrt(W1.SquaredFrobeniusNorm());
            const float n2 = sqrt(W2.SquaredFrobeniusNorm());
            const float n = std::max(n1, n2);
            if (n == 0.0f)
              continue;
            const float s = 1.0f / n;
            W1 *= s;
            W2 *= s;
          }
          scc = W1.AssertEqual(W2, verbose, str + UT::String(".W[%d][%d]", i, j)) && scc;
        }
      }
      return scc;
    }
   public:
    LA::AlignedMatrix3x3f m_W[5][5];
  };

  class Error {
   public:
    inline bool Valid() const { return m_er.Valid(); }
    inline bool Invalid() const { return m_er.Invalid(); }
    inline void Invalidate() { m_er.Invalidate(); }
    inline void Print(const bool e = false, const bool l = false) const {
      UT::PrintSeparator();
      m_er.Print(" er = ", e, l, true);
      m_ev.Print(" ev = ", e, l, true);
      m_ep.Print(" ep = ", e, l, true);
      m_eba.Print("eba = ", e, l, true);
      m_ebw.Print("ebw = ", e, l, true);
    }
    inline void Print(const std::string str, const bool e, const bool l) const {
      UT::PrintSeparator();
      const std::string _str(str.size(), ' ');
      m_er.Print( str + " er = ", e, l, true);
      m_ev.Print(_str + " ev = ", e, l, true);
      m_ep.Print(_str + " ep = ", e, l, true);
      m_eba.Print(_str + "eba = ", e, l, true);
      m_ebw.Print(_str + "ebw = ", e, l, true);
    }
   public:
    LA::AlignedVector3f m_er, m_ev, m_ep, m_eba, m_ebw;
  };
  class Jacobian {
   public:
    class Gravity {
    public:
      LA::AlignedMatrix2x3f m_JvgT, m_JpgT;
    };
    class FirstMotion {
     public:
      inline void GetTranspose(FirstMotion &JT) const {
        m_Jrbw1.GetTranspose(JT.m_Jrbw1);
        m_Jvv1.GetTranspose(JT.m_Jvv1);
        m_Jvba1.GetTranspose(JT.m_Jvba1);
        m_Jvbw1.GetTranspose(JT.m_Jvbw1);
        m_Jpv1.GetTranspose(JT.m_Jpv1);
        m_Jpba1.GetTranspose(JT.m_Jpba1);
        m_Jpbw1.GetTranspose(JT.m_Jpbw1);
      }
     public:
      LA::AlignedMatrix3x3f m_Jrbw1;
      LA::AlignedMatrix3x3f m_Jvv1, m_Jvba1, m_Jvbw1;
      LA::AlignedMatrix3x3f m_Jpv1, m_Jpba1, m_Jpbw1;
    };
    class Global : public FirstMotion {
     public:
      inline void GetTranspose(Global &JT) const {
        FirstMotion::GetTranspose(JT);
        m_Jrr1.GetTranspose(JT.m_Jrr1);
        m_Jvr1.GetTranspose(JT.m_Jvr1);
        m_Jpp1.GetTranspose(JT.m_Jpp1);
        m_Jpr1.GetTranspose(JT.m_Jpr1);
        m_Jpr2.GetTranspose(JT.m_Jpr2);
      }
     public:
      LA::AlignedMatrix3x3f m_Jrr1;
      LA::AlignedMatrix3x3f m_Jvr1;
      LA::AlignedMatrix3x3f m_Jpp1;
      LA::AlignedMatrix3x3f m_Jpr1;
      LA::AlignedMatrix3x3f m_Jpr2;
    };
    class RelativeLF : public Gravity, public Global {
     public:
      LA::AlignedMatrix3x3f m_Jvr2, m_Jvv2;
    };
    class RelativeKF : public Gravity, public FirstMotion {
     public:
      inline void GetTranspose(RelativeKF &JT) const {
        FirstMotion::GetTranspose(JT);
        m_Jrr2.GetTranspose(JT.m_Jrr2);
        m_Jvr2.GetTranspose(JT.m_Jvr2);
        m_Jvv2.GetTranspose(JT.m_Jvv2);
        m_Jpp2.GetTranspose(JT.m_Jpp2);
        m_Jpr2.GetTranspose(JT.m_Jpr2);
      }
     public:
      LA::AlignedMatrix3x3f m_Jrr2;
      LA::AlignedMatrix3x3f m_Jvr2, m_Jvv2;
      LA::AlignedMatrix3x3f m_Jpp2;
      LA::AlignedMatrix3x3f m_Jpr2;
    };
  };
  class ErrorJacobian {
   public:
    Error m_e;
    Jacobian::Global m_J;
  };
  class Factor {
   public:
    class Unitary {
     public:
      inline void Set(const LA::AlignedMatrix3x3f *Ap, const LA::AlignedMatrix3x3f *Ar,
                      const LA::AlignedMatrix3x3f *Av, const LA::AlignedMatrix3x3f *Aba,
                      const LA::AlignedMatrix3x3f *Abw, const LA::AlignedVector3f *b) {
        m_Acc.Set(Ap, Ar, b);
        m_Acm.Set(Ap + 2, Ar + 2);
        m_Amm.Set(Av + 2, Aba + 2, Abw + 2, b + 2);
      }
      inline void Set(const LA::AlignedMatrix3x3f *Av, const LA::AlignedMatrix3x3f *Aba,
                      const LA::AlignedMatrix3x3f *Abw, const LA::AlignedVector3f *b) {
        m_Acc.MakeZero();
        m_Acm.MakeZero();
        m_Amm.Set(Av, Aba, Abw, b);
      }
      inline void MakeMinus() {
        m_Acc.MakeMinus();
        m_Acm.MakeMinus();
        m_Amm.MakeMinus();
      }
      static inline void AmB(const Unitary &A, const Unitary &B, Unitary &AmB) {
        Camera::Factor::Unitary::CC::AmB(A.m_Acc, B.m_Acc, AmB.m_Acc);
        Camera::Factor::Unitary::CM::AmB(A.m_Acm, B.m_Acm, AmB.m_Acm);
        Camera::Factor::Unitary::MM::AmB(A.m_Amm, B.m_Amm, AmB.m_Amm);
      }
     public:
      Camera::Factor::Unitary::CC m_Acc;
      Camera::Factor::Unitary::CM m_Acm;
      Camera::Factor::Unitary::MM m_Amm;
    };
    class Auxiliary {
     public:
      class Global {
       public:
        void Set(const Jacobian::Global &J, const Error &e, const float w, const Weight &W,
                 const float Tpv);
        inline void Get(Unitary *A11, Unitary *A22, Camera::Factor::Binary *A12) const {

          //const LA::AlignedMatrix3x3f *A[10] = {&m_A[ 0], &m_A[ 9], &m_A[17], &m_A[24], &m_A[30],
          //                                      &m_A[35], &m_A[39], &m_A[42], &m_A[44], &m_A[45]};
          A11->Set(m_A, m_A + 9, m_A + 17, m_A + 24, m_A + 30, m_b);
          A22->Set(m_A + 40, m_A + 44, m_A + 47, m_A + 49, m_A + 50, m_b + 5);
          A12->Set(m_A + 5, m_A + 14, m_A + 22, m_A + 29, m_A + 35);
        }
       public:
        Jacobian::Global m_JT;
        Weight m_W;

        LA::AlignedMatrix3x3f m_JTW[10][5], m_A[55];
        LA::AlignedVector3f m_b[10];
      };
      class RelativeLF : public Global {
       public:
        void Set(const Jacobian::RelativeLF &J, const Error &e, const float w, const Weight &W,
                 const float Tpv);
       public:
        LA::AlignedMatrix3x3f m_Jvr2T, m_Jvv2T;
        LA::AlignedMatrix2x3f m_JTWg[5];
        LA::SymmetricMatrix2x2f m_Agg;
        LA::AlignedMatrix2x3f m_Agc[10];
        LA::Vector2f m_bg;
      };
      class RelativeKF {
       public:
        void Set(const Jacobian::RelativeKF &J, const Error &e, const float w, const Weight &W,
                 const float Tpv);

        inline void Get(Unitary *A11, Unitary *A22, Camera::Factor::Binary *A12) const {
          A11->Set(m_Ac, m_Ac + 7, m_Ac + 13, m_bc);
          A22->Set(m_Ac + 21, m_Ac + 25, m_Ac + 28, m_Ac + 30, m_Ac + 31, m_bc + 3);
          A12->Set(m_Ac + 3, m_Ac + 10, m_Ac + 16);
        }

       public:
        Jacobian::RelativeKF m_JT;
        Weight m_W;

        LA::AlignedMatrix3x3f m_JTWc[8][5], m_Ac[36];
        LA::AlignedVector3f m_bc[8];
        LA::AlignedMatrix2x3f m_JTWg[5];
        LA::SymmetricMatrix2x2f m_Agg;
        LA::AlignedMatrix2x3f m_Agc[8];
        LA::Vector2f m_bg;
      };
    };
   public:
    inline void MakeZero() { memset(this, 0, sizeof(Factor)); }
   public:
    ErrorJacobian m_Je;
    union {
      struct { float m_data[21], m_F; };
      struct { Unitary m_A11, m_A22; };
    };
  };
  class Reduction {
   public:
    Error m_e;
    float m_F, m_dF;
  };
  class ESError : public LA::AlignedVector3f {
   public:
    inline ESError() {}
    inline ESError(const LA::AlignedVector3f &e, const float s = 1.0f) {
      if (s == 1.0f) {
        *((LA::AlignedVector3f *) this) = e;
      } else {
        e.GetScaled(s, *this);
      }
    }
    inline void Print(const bool l = true) const {
      if (l) {
        UT::Print("%f %f %f", x(), y(), z());
      } else {
        UT::Print("%.2f %.2f %.2f", x(), y(), z());
      }
    }
  };
  class ES : public UT::ES<float, int> {
   public:
    inline void Initialize() {
      UT::ES<float, int>::Initialize();
      m_ESr.Initialize();
      m_ESp.Initialize();
      m_ESv.Initialize();
      m_ESba.Initialize();
      m_ESbw.Initialize();
    }
    inline void Accumulate(const Error &e, const float F, const int iFrm = -1) {
      UT::ES<float, int>::Accumulate(F, F, iFrm);
      m_ESr.Accumulate(ESError(e.m_er, UT_FACTOR_RAD_TO_DEG), -1.0f, iFrm);
      m_ESp.Accumulate(ESError(e.m_ep), -1.0f, iFrm);
      m_ESv.Accumulate(ESError(e.m_ev), -1.0f, iFrm);
      m_ESba.Accumulate(ESError(e.m_eba), -1.0f, iFrm);
      m_ESbw.Accumulate(ESError(e.m_ebw, UT_FACTOR_RAD_TO_DEG), -1.0f, iFrm);
    }
    inline void Print(const std::string str = "", const bool l = true) const {
      if (!Valid()) {
        return;
      }
      UT::ES<float, int>::Print(str + "ed = ", true, l);
      const std::string _str(str.size() + 17, ' ');
      if (m_ESr.Valid()) {
        m_ESr.Print(_str + "er  = ", false, l);
      }
      if (m_ESp.Valid()) {
        m_ESp.Print(_str + "ep  = ", false, l);
      }
      if (m_ESv.Valid()) {
        m_ESv.Print(_str + "ev  = ", false, l);
      }
      if (m_ESba.Valid()) {
        m_ESba.Print(_str + "eba = ", false, l);
      }
      if (m_ESbw.Valid()) {
        m_ESbw.Print(_str + "ebw = ", false, l);
      }
    }
  public:
    UT::ES<ESError, int> m_ESr, m_ESp, m_ESv, m_ESba, m_ESbw;
  };
 public:

  inline bool Valid() const { return m_RT.Valid(); }
  inline bool Invalid() const { return m_RT.Invalid(); }
  inline void Invalidate() { m_RT.Invalidate(); }

  inline Rotation3D GetRotationState(const Camera &C1, const Camera &C2) const {
    return Rotation3D(C1.m_T) / C2.m_T;
  }
  inline Rotation3D GetRotationMeasurement(const Camera &C1, const float eps) const {
    return m_RT / Rotation3D(m_Jrbw * (C1.m_bw - m_bw), eps);
  }
  inline Rotation3D GetRotationMeasurement(const LA::AlignedVector3f &dbw, const float eps) const {
    return m_RT / Rotation3D(m_Jrbw * dbw, eps);
  }
  inline LA::AlignedVector3f GetRotationError(const Camera &C1, const Camera &C2,
                                              const float eps) const {
    const Rotation3D eR = GetRotationMeasurement(C1, eps) / GetRotationState(C1, C2);
    return eR.GetRodrigues(eps);
  }

  inline LA::AlignedVector3f GetVelocityMeasurement(const Camera &C1) const {
    return m_v + m_Jvba * (C1.m_ba - m_ba) + m_Jvbw * (C1.m_bw - m_bw);
  }
  inline LA::AlignedVector3f GetVelocityState(const Camera &C1, const Camera &C2) const {
    LA::AlignedVector3f dv = C2.m_v - C1.m_v;
    if (!IMU_GRAVITY_EXCLUDED) {
      dv.z() += IMU_GRAVITY_MAGNITUDE * m_Tvg;
    }
    return C1.m_T.GetAppliedRotation(dv);
  }
  inline LA::AlignedVector3f GetVelocityError(const Camera &C1, const Camera &C2) const {
    return GetVelocityState(C1, C2) - GetVelocityMeasurement(C1);
  }

  inline LA::AlignedVector3f GetPositionState(const Camera &C1, const Camera &C2,
                                              const Point3D &pu) const {
    LA::AlignedVector3f dp = C2.m_p - C1.m_p - C1.m_v * m_Tpv;
    if (!IMU_GRAVITY_EXCLUDED) {
      dp.z() += IMU_GRAVITY_MAGNITUDE * m_Tpg;
    }
    dp += C2.m_T.GetAppliedRotationInversely(pu);
    dp = C1.m_T.GetAppliedRotation(dp);
    dp -= pu;
    return dp;
  }
  inline LA::AlignedVector3f GetPositionMeasurement(const Camera &C1) const {
    return m_p + m_Jpba * (C1.m_ba - m_ba) + m_Jpbw * (C1.m_bw - m_bw);
  }
  inline LA::AlignedVector3f GetPositionError(const Camera &C1, const Camera &C2,
                                              const Point3D &pu) const {
    return GetPositionState(C1, C2, pu) - GetPositionMeasurement(C1);
  }

  inline void GetError(const Camera &C1, const Camera &C2, const Point3D &pu, Error &e,
                       const float eps) const {
    e.m_er = GetRotationError(C1, C2, eps);
    e.m_ev = GetVelocityError(C1, C2);
    e.m_ep = GetPositionError(C1, C2, pu);
    e.m_eba = C1.m_ba - C2.m_ba;
    e.m_ebw = C1.m_bw - C2.m_bw;
  }
  inline Error GetError(const Camera &C1, const Camera &C2, const Point3D &pu,
                        const float eps) const {
    Error e;
    GetError(C1, C2, pu, e, eps);
    return e;
  }
  static inline void GetError(const ErrorJacobian &Je, const LA::AlignedVector3f *xp1,
                              const LA::AlignedVector3f *xr1, const LA::AlignedVector3f *xv1,
                              const LA::AlignedVector3f *xba1, const LA::AlignedVector3f *xbw1, 
                              const LA::AlignedVector3f *xp2, const LA::AlignedVector3f *xr2,
                              const LA::AlignedVector3f *xv2, const LA::AlignedVector3f *xba2,
                              const LA::AlignedVector3f *xbw2, Error &e) {
    e = Je.m_e;
    if (xp1) {
      if (xp2) {
        const LA::AlignedVector3f dxp = *xp1 - *xp2;
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpp1, dxp, (float *) &e.m_ep);
      } else {
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpp1, *xp1, (float *) &e.m_ep);
      }
    } else if (xp2) {
      LA::AlignedMatrix3x3f::SubtractAbFrom(Je.m_J.m_Jpp1, *xp2, (float *) &e.m_ep);
    }
    if (xr1) {
      if (xr2) {
        const LA::AlignedVector3f dxr = *xr1 - *xr2;
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jrr1, dxr, (float *) &e.m_er);
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpr2, *xr2, (float *) &e.m_ep);
      } else {
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jrr1, *xr1, (float *) &e.m_er);
      }
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jvr1, *xr1, (float *) &e.m_ev);
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpr1, *xr1, (float *) &e.m_ep);
    } else if (xr2) {
      LA::AlignedMatrix3x3f::SubtractAbFrom(Je.m_J.m_Jrr1, *xr2, (float *) &e.m_er);
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpr2, *xr2, (float *) &e.m_ep);
    }
    if (xv1) {
      if (xv2) {
        const LA::AlignedVector3f dxv = *xv1 - *xv2;
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jvv1, dxv, (float *) &e.m_ev);
      } else {
        LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jvv1, *xv1, (float *) &e.m_ev);
      }
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpv1, *xv1, (float *) &e.m_ep);
    } else if (xv2) {
      LA::AlignedMatrix3x3f::SubtractAbFrom(Je.m_J.m_Jvv1, *xv2, (float *) &e.m_ev);
    }
    if (xba1) {
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jvba1, *xba1, (float *) &e.m_ev);
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpba1, *xba1, (float *) &e.m_ep);
      e.m_eba += *xba1;
    }
    if (xbw1) {
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jrbw1, *xbw1, (float *) &e.m_er);
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jvbw1, *xbw1, (float *) &e.m_ev);
      LA::AlignedMatrix3x3f::AddAbTo(Je.m_J.m_Jpbw1, *xbw1, (float *) &e.m_ep);
      e.m_ebw += *xbw1;
    }
    if (xba2) {
      e.m_eba -= *xba2;
    }
    if (xbw2) {
      e.m_ebw -= *xbw2;
    }
  }
  inline void GetError(const Jacobian::RelativeLF &J, const LA::Vector2f &xg,
                       const LA::AlignedVector3f &xp1, const LA::AlignedVector3f &xr1,
                       const LA::AlignedVector3f &xv1, const LA::AlignedVector3f &xba1,
                       const LA::AlignedVector3f &xbw1, const LA::AlignedVector3f &xp2,
                       const LA::AlignedVector3f &xr2, const LA::AlignedVector3f &xv2,
                       const LA::AlignedVector3f &xba2, const LA::AlignedVector3f &xbw2,
                       Error &e) const {
    GetError(J, xg, xp1, xr1, xv1, xba1, xbw1, xp2, xr2, xv2, xba2, xbw2, m_Tpv, e);
  }
  static inline void GetError(const Jacobian::RelativeLF &J, const LA::Vector2f &xg,
                              const LA::AlignedVector3f &xp1, const LA::AlignedVector3f &xr1,
                              const LA::AlignedVector3f &xv1, const LA::AlignedVector3f &xba1,
                              const LA::AlignedVector3f &xbw1, const LA::AlignedVector3f &xp2,
                              const LA::AlignedVector3f &xr2, const LA::AlignedVector3f &xv2,
                              const LA::AlignedVector3f &xba2, const LA::AlignedVector3f &xbw2,
                              const float Tpv, Error &e) {
    const LA::AlignedVector3f dxr = xr1 - xr2;
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jrr1, dxr, (float *) &e.m_er);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jrbw1, xbw1, (float *) &e.m_er);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvr1, xr1, (float *) &e.m_ev);
    //LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvv1, xv1, (float *) &e.m_ev);
    e.m_ev -= xv1;
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvba1, xba1, (float *) &e.m_ev);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvbw1, xbw1, (float *) &e.m_ev);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvr2, xr2, (float *) &e.m_ev);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvv2, xv2, (float *) &e.m_ev);
    LA::AlignedMatrix2x3f::AddATbTo(J.m_JvgT, xg, e.m_ev);
    const LA::AlignedVector3f dxp = xp1 - xp2;
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpp1, dxp, (float *) &e.m_ep);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpr1, xr1, (float *) &e.m_ep);
    //LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpv1, xv1, e.m_ep);
    e.m_ep -= (xv1 * Tpv);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpba1, xba1, (float *) &e.m_ep);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpbw1, xbw1, (float *) &e.m_ep);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpr2, xr2, (float *) &e.m_ep);
    LA::AlignedMatrix2x3f::AddATbTo(J.m_JpgT, xg, e.m_ep);
    e.m_eba += xba1;
    e.m_eba -= xba2;
    e.m_ebw += xbw1;
    e.m_ebw -= xbw2;
  }
  inline void GetError(const Jacobian::RelativeKF &J, const LA::Vector2f &xg,
                       const LA::AlignedVector3f &xv1, const LA::AlignedVector3f &xba1,
                       const LA::AlignedVector3f &xbw1, const LA::AlignedVector3f &xp2,
                       const LA::AlignedVector3f &xr2, const LA::AlignedVector3f &xv2,
                       const LA::AlignedVector3f &xba2, const LA::AlignedVector3f &xbw2,
                       Error &e) const {
    GetError(J, xg, xv1, xba1, xbw1, xp2, xr2, xv2, xba2, xbw2, m_Tpv, e);
  }

  static inline void GetError(const Jacobian::RelativeKF &J, const LA::Vector2f &xg,
                              const LA::AlignedVector3f &xv1, const LA::AlignedVector3f &xba1,
                              const LA::AlignedVector3f &xbw1, const LA::AlignedVector3f &xp2,
                              const LA::AlignedVector3f &xr2, const LA::AlignedVector3f &xv2,
                              const LA::AlignedVector3f &xba2, const LA::AlignedVector3f &xbw2,
                              const float Tpv, Error &e) {
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jrbw1, xbw1, (float *) &e.m_er);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jrr2, xr2, (float *) &e.m_er);
    //LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvv1, xv1, (float *) &e.m_ev);
    e.m_ev -= xv1;
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvba1, xba1, (float *) &e.m_ev);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvbw1, xbw1, (float *) &e.m_ev);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvr2, xr2, (float *) &e.m_ev);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jvv2, xv2, (float *) &e.m_ev);
    LA::AlignedMatrix2x3f::AddATbTo(J.m_JvgT, xg, e.m_ev);
    //LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpv1, xv1, e.m_ep);
    e.m_ep -= (xv1 * Tpv);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpba1, xba1, (float *) &e.m_ep);
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpbw1, xbw1, (float *) &e.m_ep);
    //LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpp2, xp2, (float *) &e.m_ep);
    e.m_ep += xp2;
    LA::AlignedMatrix3x3f::AddAbTo(J.m_Jpr2, xr2, (float *) &e.m_ep);
    LA::AlignedMatrix2x3f::AddATbTo(J.m_JpgT, xg, e.m_ep);
    e.m_eba += xba1;
    e.m_eba -= xba2;
    e.m_ebw += xbw1;
    e.m_ebw -= xbw2;
  }
  inline void GetErrorJacobian(const Camera &C1, const Camera &C2, const Point3D &pu,
                               Error *e, Jacobian::Global *J, const float eps) const {
    const Rotation3D dR = GetRotationState(C1, C2);
    const LA::AlignedVector3f drbw = m_Jrbw * (C1.m_bw - m_bw);
    const Rotation3D eR = m_RT / Rotation3D(drbw, eps) / dR;
    eR.GetRodrigues(e->m_er, eps);
    Rotation3D::GetRodriguesJacobianInverse(e->m_er, J->m_Jrr1, eps);
    Rotation3D::GetRodriguesJacobian(drbw.GetMinus(), J->m_Jrbw1, eps);
    J->m_Jrr1.MakeMinus();
    J->m_Jrbw1 = J->m_Jrr1 * m_RT * J->m_Jrbw1 * m_Jrbw;
    J->m_Jrr1 = J->m_Jrr1 * eR * C1.m_T;

    e->m_ev = GetVelocityState(C1, C2);
    SkewSymmetricMatrix::AB(e->m_ev, C1.m_T, J->m_Jvr1);
    e->m_ev -= GetVelocityMeasurement(C1);
    C1.m_T.GetMinus(J->m_Jvv1);
    m_Jvba.GetMinus(J->m_Jvba1);
    m_Jvbw.GetMinus(J->m_Jvbw1);

    e->m_ep = C2.m_p - C1.m_p - C1.m_v * m_Tpv;
    if (!IMU_GRAVITY_EXCLUDED) {
      e->m_ep.z() += IMU_GRAVITY_MAGNITUDE * m_Tpg;
    }
    e->m_ep = C1.m_T.GetAppliedRotation(e->m_ep);
    const LA::AlignedVector3f R21pu = dR.GetApplied(pu);

    e->m_ep += R21pu;
    SkewSymmetricMatrix::AB(e->m_ep, C1.m_T, J->m_Jpr1);
    e->m_ep -= pu;
    e->m_ep -= GetPositionMeasurement(C1);
    C1.m_T.GetMinus(J->m_Jpp1);
    C1.m_T.GetScaled(-m_Tpv, J->m_Jpv1);
    m_Jpba.GetMinus(J->m_Jpba1);
    m_Jpbw.GetMinus(J->m_Jpbw1);
    SkewSymmetricMatrix::ATB(R21pu, C1.m_T, J->m_Jpr2);

    e->m_eba = C1.m_ba - C2.m_ba;
    e->m_ebw = C1.m_bw - C2.m_bw;

    //J->m_Jpp1.MakeZero();
    //J->m_Jpr1.MakeZero();
    //J->m_Jpba1.MakeZero();
    //J->m_Jpbw1.MakeZero();
    //J->m_Jpr2.MakeZero();
  }
  inline void GetErrorJacobian(const Camera &C1, const Camera &C2, const Point3D &pu,
                               const Rotation3D &Rg, Error *e, Jacobian::RelativeLF *J,
                               const float eps) const {
    const Rotation3D dR = GetRotationState(C1, C2);
    const LA::AlignedVector3f drbw = m_Jrbw * (C1.m_bw - m_bw);
    const Rotation3D eR = m_RT / Rotation3D(drbw, eps) / dR;
    eR.GetRodrigues(e->m_er, eps);
    Rotation3D::GetRodriguesJacobianInverse(e->m_er, J->m_Jrr1, eps);
    Rotation3D::GetRodriguesJacobian(drbw.GetMinus(), J->m_Jrbw1, eps);
    J->m_Jrr1.MakeMinus();
    J->m_Jrbw1 = J->m_Jrr1 * m_RT * J->m_Jrbw1 * m_Jrbw;
    J->m_Jrr1 = J->m_Jrr1 * eR;
    const Rotation3D R1T = Rg / C1.m_T;
    J->m_Jrr1 = LA::AlignedMatrix3x3f::GetABT(J->m_Jrr1, R1T);

    C1.m_T.ApplyRotation(C2.m_v, e->m_ev);
    J->m_Jvv2 = dR;
    LA::AlignedMatrix3x3f::Ab(J->m_Jvv2, pu, (float *) &e->m_ep);
    //const Rotation3D R1 = Rotation3D(C1.m_T) / Rg;
    const Rotation3D R1 = R1T.GetTranspose();
    SkewSymmetricMatrix::ATB(e->m_ev, R1, J->m_Jvr2);
    SkewSymmetricMatrix::ATB(e->m_ep, R1, J->m_Jpr2);
    if (IMU_GRAVITY_EXCLUDED) {
      J->m_JvgT.Invalidate();
      J->m_JpgT.Invalidate();
    } else {
      const LA::AlignedVector3f g1 = C1.m_T.GetColumn2();
      const LA::AlignedVector3f dv = g1 * (m_Tvg * IMU_GRAVITY_MAGNITUDE);
      e->m_ev += dv;
      SkewSymmetricMatrix::ATBT(C1.m_T, dv, J->m_JvgT);
      const LA::AlignedVector3f dp = g1 * (m_Tpg * IMU_GRAVITY_MAGNITUDE);
      e->m_ep += dp;
      SkewSymmetricMatrix::ATBT(C1.m_T, dp, J->m_JpgT);
    }
    SkewSymmetricMatrix::AB(e->m_ev, R1, J->m_Jvr1);
    const LA::AlignedVector3f v1 = C1.m_T.GetAppliedRotation(C1.m_v);
    e->m_ev -= v1;
    e->m_ev -= GetVelocityMeasurement(C1);
    //J->m_Jvv1.MakeDiagonal(-1.0);
    m_Jvba.GetMinus(J->m_Jvba1);
    m_Jvbw.GetMinus(J->m_Jvbw1);

    e->m_ep += C1.m_T.GetAppliedRotation(C2.m_p - C1.m_p);
    SkewSymmetricMatrix::AB(e->m_ep, R1, J->m_Jpr1);
    e->m_ep -= v1 * m_Tpv;
    e->m_ep -= pu;
    e->m_ep -= GetPositionMeasurement(C1);
    R1.GetMinus(J->m_Jpp1);
    //J->m_Jpv1.MakeDiagonal(-m_Tpv);
    m_Jpba.GetMinus(J->m_Jpba1);
    m_Jpbw.GetMinus(J->m_Jpbw1);

    e->m_eba = C1.m_ba - C2.m_ba;
    e->m_ebw = C1.m_bw - C2.m_bw;
  }
  inline void GetErrorJacobian(const Camera &C1, const Camera &C2, const Point3D &pu,
                               Error *e, Jacobian::RelativeKF *J, const float eps) const {
    const Rotation3D dR = GetRotationState(C1, C2);
    const LA::AlignedVector3f drbw = m_Jrbw * (C1.m_bw - m_bw);
    const Rotation3D eR = m_RT / Rotation3D(drbw, eps) / dR;
    eR.GetRodrigues(e->m_er, eps);
    Rotation3D::GetRodriguesJacobianInverse(e->m_er, J->m_Jrr2, eps);
    Rotation3D::GetRodriguesJacobian(drbw.GetMinus(), J->m_Jrbw1, eps);
    J->m_Jrbw1 = J->m_Jrr2 * m_RT * J->m_Jrbw1 * m_Jrbw;
    J->m_Jrbw1.MakeMinus();
    J->m_Jrr2 = J->m_Jrr2 * eR;

    C1.m_T.ApplyRotation(C2.m_v, e->m_ev);
    J->m_Jvv2 = dR;
    LA::AlignedMatrix3x3f::Ab(J->m_Jvv2, pu, (float *) &e->m_ep);
    SkewSymmetricMatrix::GetTranspose(e->m_ev, J->m_Jvr2);
    SkewSymmetricMatrix::GetTranspose(e->m_ep, J->m_Jpr2);
    if (IMU_GRAVITY_EXCLUDED) {
      J->m_JvgT.Invalidate();
      J->m_JpgT.Invalidate();
    } else {
      const LA::AlignedVector3f g1 = C1.m_T.GetColumn2();
      const LA::AlignedVector3f dv = g1 * (m_Tvg * IMU_GRAVITY_MAGNITUDE);
      e->m_ev += dv;
      SkewSymmetricMatrix::ATBT(C1.m_T, dv, J->m_JvgT);
      const LA::AlignedVector3f dp = g1 * (m_Tpg * IMU_GRAVITY_MAGNITUDE);
      e->m_ep += dp;
      SkewSymmetricMatrix::ATBT(C1.m_T, dp, J->m_JpgT);
    }
    const LA::AlignedVector3f v1 = C1.m_T.GetAppliedRotation(C1.m_v);
    e->m_ev -= v1;
    e->m_ev -= GetVelocityMeasurement(C1);
    //J->m_Jvv1.MakeDiagonal(-1.0);
    m_Jvba.GetMinus(J->m_Jvba1);
    m_Jvbw.GetMinus(J->m_Jvbw1);

    e->m_ep += C1.m_T.GetAppliedRotation(C2.m_p - C1.m_p);
    e->m_ep -= v1 * m_Tpv;
    e->m_ep -= pu;
    e->m_ep -= GetPositionMeasurement(C1);
    //J->m_Jpv1.MakeDiagonal(-m_Tpv);
    //J->m_Jpp2.MakeIdentity();
    m_Jpba.GetMinus(J->m_Jpba1);
    m_Jpbw.GetMinus(J->m_Jpbw1);

    e->m_eba = C1.m_ba - C2.m_ba;
    e->m_ebw = C1.m_bw - C2.m_bw;
  }
  inline void GetFactor(const float w, const Camera &C1, const Camera &C2, const Point3D &pu,
                        Factor *A, Camera::Factor::Binary *A12, Factor::Auxiliary::Global *U,
                        const float eps) const {
    GetErrorJacobian(C1, C2, pu, &A->m_Je.m_e, &A->m_Je.m_J, eps);
    A->m_F = GetCost(w, A->m_Je.m_e);
    U->Set(A->m_Je.m_J, A->m_Je.m_e, w, m_W, m_Tpv);
    U->Get(&A->m_A11, &A->m_A22, A12);
  }
  inline void GetFactor(const float w, const Camera &C1, const Camera &C2, const Point3D &pu,
                        const Rotation3D &Rg, Error *e, Jacobian::RelativeLF *J,
                        Factor::Auxiliary::RelativeLF *U, const float eps) const {
    GetErrorJacobian(C1, C2, pu, Rg, e, J, eps);
    U->Set(*J, *e, w, m_W, m_Tpv);
  }
  inline void GetFactor(const float w, const Camera &C1, const Camera &C2, const Point3D &pu,
                        Error *e, Jacobian::RelativeKF *J, Factor::Auxiliary::RelativeKF *U,
                        const float eps) const {
    GetErrorJacobian(C1, C2, pu, e, J, eps);
    U->Set(*J, *e, w, m_W, m_Tpv);
  }
  inline float GetCost(const float w, const Error &e) const {
    return GetCost(w, m_W, e);
  }
  static inline float GetCost(const float w, const Weight &W, const Error &e) {
    LA::AlignedVector3f We;
    float F = 0.0f;
    const LA::AlignedVector3f *_e = (LA::AlignedVector3f *) &e;
    for (int i = 0; i < 5; ++i) {
      We.MakeZero();
      const LA::AlignedMatrix3x3f *Wi = W[i];
      for (int j = 0; j < 5; ++j) {
        LA::AlignedMatrix3x3f::AddAbTo(Wi[j], _e[j], (float *) &We);
      }
      F += _e[i].Dot(We);
    }
    return w * F;
  }
  inline float GetCost(const float w, const ErrorJacobian &Je, const LA::AlignedVector3f *xp1,
                       const LA::AlignedVector3f *xr1, const LA::AlignedVector3f *xv1,
                       const LA::AlignedVector3f *xba1, const LA::AlignedVector3f *xbw1,
                       const LA::AlignedVector3f *xp2, const LA::AlignedVector3f *xr2,
                       const LA::AlignedVector3f *xv2, const LA::AlignedVector3f *xba2,
                       const LA::AlignedVector3f *xbw2, Error &e) const {
    GetError(Je, xp1, xr1, xv1, xba1, xbw1, xp2, xr2, xv2, xba2, xbw2, e);
    return GetCost(w, e);
  }
  inline void GetReduction(const float w, const Factor &A, const Camera &C1, const Camera &C2,
                           const Point3D &pu, const LA::AlignedVector3f *xp1,
                           const LA::AlignedVector3f *xr1, const LA::AlignedVector3f *xv1,
                           const LA::AlignedVector3f *xba1, const LA::AlignedVector3f *xbw1,
                           const LA::AlignedVector3f *xp2, const LA::AlignedVector3f *xr2, 
                           const LA::AlignedVector3f *xv2, const LA::AlignedVector3f *xba2,
                           const LA::AlignedVector3f *xbw2, Reduction &Ra, Reduction &Rp,
                           const float eps) const {
    GetError(C1, C2, pu, Ra.m_e, eps);
    GetError(A.m_Je, xp1, xr1, xv1, xba1, xbw1, xp2, xr2, xv2, xba2, xbw2, Rp.m_e);
    Ra.m_dF = A.m_F - (Ra.m_F = GetCost(w, Ra.m_e));
    Rp.m_dF = A.m_F - (Rp.m_F = GetCost(w, Rp.m_e));
  }

  inline bool AssertEqual(const Delta &D, const int verbose = 1, const std::string str = "",
    const bool normWeight = true) const {
    bool scc = true;
    scc = m_ba.AssertEqual(D.m_ba, verbose, str + ".m_ba") && scc;
    scc = m_bw.AssertEqual(D.m_bw, verbose, str + ".m_bw") && scc;
    scc = m_RT.AssertEqual(D.m_RT, verbose, str + ".m_R") && scc;
    scc = m_v.AssertEqual(D.m_v, verbose, str + ".m_v") && scc;
    scc = m_p.AssertEqual(D.m_p, verbose, str + ".m_p") && scc;
    scc = m_Jrbw.AssertEqual(D.m_Jrbw, verbose, str + ".m_Jrbw") && scc;
    scc = m_Jvba.AssertEqual(D.m_Jvba, verbose, str + ".m_Jvba") && scc;
    scc = m_Jvbw.AssertEqual(D.m_Jvbw, verbose, str + ".m_Jvbw") && scc;
    scc = m_Jpba.AssertEqual(D.m_Jpba, verbose, str + ".m_Jpba") && scc;
    scc = m_Jpbw.AssertEqual(D.m_Jpbw, verbose, str + ".m_Jpbw") && scc;
    scc = m_W.AssertEqual(D.m_W, verbose, str + ".m_W", normWeight) && scc;
    scc = UT::AssertEqual(m_Tvg, D.m_Tvg, verbose, str + ".m_Tvg") && scc;
    scc = UT::AssertEqual(m_Tpv, D.m_Tpv, verbose, str + ".m_Tpv") && scc;
    scc = UT::AssertEqual(m_Tpg, D.m_Tpg, verbose, str + ".m_Tpg") && scc;
    return scc;
  }
  inline void Print(const bool e = false) const {
    UT::PrintSeparator();
    m_ba.Print("  ba = ", e, false, true);
    m_bw.Print("  bw = ", e, false, true);
    m_RT.Print("  RT = ", e);
    m_v.Print("   v = ", e, false, true);
    m_p.Print("   p = ", e, false, true);
    m_Jrbw.Print("Jrbw = ", e);
    m_Jvba.Print("Jvba = ", e);
    m_Jvbw.Print("Jvbw = ", e);
    m_Jpba.Print("Jpba = ", e);
    m_Jpbw.Print("Jpbw = ", e);
    if (e) {
      UT::Print(" Tvg = %e\n", m_Tvg);
    } else {
      UT::Print(" Tvg = %f\n", m_Tvg);
    }
  }
 public:
  Measurement m_u1, m_u2;
  LA::AlignedVector3f m_ba, m_bw;
  Rotation3D m_RT;
  LA::AlignedVector3f m_v, m_p;
  LA::AlignedMatrix3x3f m_Jrbw, m_Jvba, m_Jvbw, m_Jpba, m_Jpbw;
  Weight m_W;
  float m_Tvg, m_Tpv, m_Tpg, m_r;
};

void InitializeCamera(const AlignedVector<Measurement> &us, Camera &C);
void PreIntegrate(const AlignedVector<Measurement> &us, const float t1, const float t2,
                  const Camera &C1, Delta *D, AlignedVector<float> *work, const bool jac/* = true*/,
                  const Measurement *u1/* = NULL*/, const Measurement *u2/* = NULL*/,
                  const float eps);
void PreIntegrate(const Measurement *us, const int N, const float t1, const float t2,
                  const Camera &C1, Delta *D, AlignedVector<float> *work, const bool jac/* = true*/,
                  const Measurement *u1/* = NULL*/, const Measurement *u2/* = NULL*/,
                  const float eps);
void Propagate(const Point3D &pu, const Delta &D, const Camera &C1, Camera &C2,
               const float eps);

};

#endif
