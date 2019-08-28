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
#ifndef _MATRIX_NxN_H_
#define _MATRIX_NxN_H_

#include "VectorN.h"
#include "Matrix4x4.h"
#include "Matrix2x8.h"
#include "Matrix3x6.h"
#include "Matrix6x6.h"
#include "Matrix8x8.h"
#include "Matrix9x9.h"
#include "AlignedMatrix.h"

namespace LA {
template<int M, int N>
class SIMD_ALIGN_DECLSPEC AlignedMatrixMxNf {
 public:
  class Row : public AlignedVectorNf<N> {
  };

 public:
  inline const Row& operator() (const int i) const { return m_rows[i]; }
  inline       Row& operator() (const int i)       { return m_rows[i]; }
  inline const float* operator[] (const int i) const { return m_rows[i].m_data; }
  inline       float* operator[] (const int i)       { return m_rows[i].m_data; }
  inline const float& operator() (const int i, const int j) const { return m_rows[i].m_data[j]; }
  inline       float& operator() (const int i, const int j)       { return m_rows[i].m_data[j]; }

  //inline AlignedMatrixMxNf() {}
  //inline AlignedMatrixMxNf(const AlignedMatrixMxNf<M, N> &B) { *this = B; }

  inline void operator += (const AlignedMatrixMxNf<M, N> &B) {
    for (int i = 0; i < M; ++i)
      m_rows[i] += B.m_rows[i];
  }
  inline void operator -= (const AlignedMatrixMxNf<M, N> &B) {
    for (int i = 0; i < M; ++i)
      m_rows[i] -= B.m_rows[i];
  }
  inline void operator *= (const float s) {
    const xp128f _s = xp128f::get(s);
    for (int i = 0; i < M; ++i)
      m_rows[i] *= _s;
  }
  inline void operator *= (const xp128f &s) {
    for (int i = 0; i < M; ++i)
      m_rows[i] *= s;
  }
  inline AlignedMatrixMxNf<M, N> operator + (const AlignedMatrixMxNf<M, N> &B) const {
    AlignedMatrixMxNf<M, N> _ApB;
    ApB(*this, B, _ApB);
    return _ApB;
  }
  inline AlignedMatrixMxNf<M, N> operator - (const AlignedMatrixMxNf<M, N> &B) const {
    AlignedMatrixMxNf<M, N> _AmB;
    AmB(*this, B, _AmB);
    return _AmB;
  }

  inline void MakeZero() { memset(this, 0, sizeof(AlignedMatrixMxNf<M, N>)); }
  inline void MakeMinus() {
    for (int i = 0; i < M; ++i)
      m_rows[i].MakeMinus();
  }

  inline void MakeIdentity() {
    MakeZero();
    const int min_MN = std::min(M, N);
    for (int i = 0; i < min_MN; ++i)
      m_rows[i].m_data[i] = 1.0f;
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix2x3f &B) {
    memcpy(m_rows[i].m_data + j, &B.m00(), 12);
    memcpy(m_rows[i + 1].m_data + j, &B.m10(), 12);
  }
  inline void SetBlock(const int i, const int j, const Matrix2x3f &B) {
    memcpy(m_rows[i].m_data + j, B[0], 12);
    memcpy(m_rows[i + 1].m_data + j, B[1], 12);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix3x3f &B) {
    memcpy(m_rows[i].m_data + j, &B.m00(), 12);
    memcpy(m_rows[i + 1].m_data + j, &B.m10(), 12);
    memcpy(m_rows[i + 2].m_data + j, &B.m20(), 12);
  }
  inline void SetBlock(const int i, const int j, const Matrix3x3f &B) {
    memcpy(m_rows[i].m_data + j, B[0], 12);
    memcpy(m_rows[i + 1].m_data + j, B[1], 12);
    memcpy(m_rows[i + 2].m_data + j, B[2], 12);
  }
  inline void SetBlock(const int i, const int j, const AlignedVector3f &v) {
    m_rows[i].m_data[j] = v.v0();
    m_rows[i + 1].m_data[j] = v.v1();
    m_rows[i + 2].m_data[j] = v.v2();
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix6x6f &B) {
    memcpy(m_rows[i] + j, B[0], 24);
    memcpy(m_rows[i + 1] + j, B[1], 24);
    memcpy(m_rows[i + 2] + j, B[2], 24);
    memcpy(m_rows[i + 3] + j, B[3], 24);
    memcpy(m_rows[i + 4] + j, B[4], 24);
    memcpy(m_rows[i + 5] + j, B[5], 24);
  }
  inline void SetBlockDiagonal(const int i, const SymmetricMatrix3x3f &B) {
    int k = i;
    memcpy(m_rows[k] + k, &B.m00(), 12);  ++k;
    memcpy(m_rows[k] + k, &B.m11(), 8);   ++k;
    m_rows[k][k] = B.m22();
  }
  inline void SetBlockDiagonal(const int i, const SymmetricMatrix6x6f &B) {
    int k = i;
    memcpy(m_rows[k] + k, &B.m00(), 24);  ++k;
    memcpy(m_rows[k] + k, &B.m11(), 20);  ++k;
    memcpy(m_rows[k] + k, &B.m22(), 16);  ++k;
    memcpy(m_rows[k] + k, &B.m33(), 12);  ++k;
    memcpy(m_rows[k] + k, &B.m44(), 8);   ++k;
    m_rows[k][k] = B.m55();
  }

  inline void GetBlock(const int i, const int j, AlignedMatrix3x3f &B) const {
    memcpy(&B.m00(), m_rows[i].m_data + j, 12);
    memcpy(&B.m10(), m_rows[i + 1].m_data + j, 12);
    memcpy(&B.m20(), m_rows[i + 2].m_data + j, 12);
  }
  inline void GetBlock(const int i, const int j, Matrix2x3f &B) const {
    memcpy(B[0], m_rows[i].m_data + j, 12);
    memcpy(B[1], m_rows[i + 1].m_data + j, 12);
  }
  inline void GetBlock(const int i, const int j, Matrix3x3f &B) const {
    memcpy(B[0], m_rows[i].m_data + j, 12);
    memcpy(B[1], m_rows[i + 1].m_data + j, 12);
    memcpy(B[2], m_rows[i + 2].m_data + j, 12);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix6x6f &B) const {
    memcpy(B[0], m_rows[i] + j, 24);
    memcpy(B[1], m_rows[i + 1] + j, 24);
    memcpy(B[2], m_rows[i + 2] + j, 24);
    memcpy(B[3], m_rows[i + 3] + j, 24);
    memcpy(B[4], m_rows[i + 4] + j, 24);
    memcpy(B[5], m_rows[i + 5] + j, 24);
  }
  inline void GetBlock(const int i, const int j, Vector6f &v) const {
    v.v0() = m_rows[i][j];
    v.v1() = m_rows[i + 1][j];
    v.v2() = m_rows[i + 2][j];
    v.v3() = m_rows[i + 3][j];
    v.v4() = m_rows[i + 4][j];
    v.v5() = m_rows[i + 5][j];
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix3x3f &B) const {
    int k = i;
    memcpy(&B.m00(), m_rows[k] + k, 12);  ++k;
    memcpy(&B.m11(), m_rows[k] + k, 8);   ++k;
    B.m22() = m_rows[k][k];
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix3x3d &B) const {
    B.m00() = double(m_rows[i][i]);
    B.m01() = double(m_rows[i][i + 1]);
    B.m02() = double(m_rows[i][i + 2]);
    B.m11() = double(m_rows[i + 1][i + 1]);
    B.m12() = double(m_rows[i + 1][i + 2]);
    B.m22() = double(m_rows[i + 2][i + 2]);
  }
  inline void GetBlockDiagonal(const int i, AlignedMatrix3x3f &B) const {
    memcpy(&B.m00(), m_rows[i] + i, 12);
    memcpy(&B.m10(), m_rows[i + 1] + i, 12);
    memcpy(&B.m20(), m_rows[i + 2] + i, 12);
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix6x6f &B) const {
    int k = i;
    memcpy(&B.m00(), m_rows[k] + k, 24);  ++k;
    memcpy(&B.m11(), m_rows[k] + k, 20);  ++k;
    memcpy(&B.m22(), m_rows[k] + k, 16);  ++k;
    memcpy(&B.m33(), m_rows[k] + k, 12);  ++k;
    memcpy(&B.m44(), m_rows[k] + k, 8);   ++k;
    B.m55() = m_rows[k][k];
  }

  inline void Increase(const int i, const int j, const AlignedMatrix3x3f &B) {
    float *A;
    A = m_rows[i].m_data + j;     A[0] += B.m00();  A[1] += B.m01();  A[2] += B.m02();
    A = m_rows[i + 1].m_data + j; A[0] += B.m10();  A[1] += B.m11();  A[2] += B.m12();
    A = m_rows[i + 2].m_data + j; A[0] += B.m20();  A[1] += B.m21();  A[2] += B.m22();
  }
  inline void IncreaseDiagonal(const int i, const SymmetricMatrix3x3f &B) {
    float *A;
    A = m_rows[i].m_data + i;     A[0] += B.m00();  A[1] += B.m01();  A[2] += B.m02();
    A = m_rows[i + 1].m_data + i;                   A[1] += B.m11();  A[2] += B.m12();
    A = m_rows[i + 2].m_data + i;                                     A[2] += B.m22();
  }

  inline void Transpose() {
    for (int i = 0; i < M; ++i) {
      for (int j = i + 1; j < M; ++j) {
        UT_SWAP(m_rows[i][j], m_rows[j][i]);
      }
    }
  }
  inline AlignedMatrixMxNf<N, M> GetTranspose() const {
    AlignedMatrixMxNf<N, M> MT;
    GetTranspose(MT);
    return MT;
  }
  inline void GetTranspose(AlignedMatrixMxNf<N, M> &_M) const {
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        _M[j][i] = m_rows[i][j];
      }
    }
  }

  inline void SetLowerFromUpper() {
    const int Nmin = std::min(M, N);
    for (int i = 0; i < Nmin; ++i) {
      for (int j = i; j < Nmin; ++j) {
        m_rows[j][i] = m_rows[i][j];
      }
    }
  }
  
  inline bool Valid() const { return m_rows[0][0] != FLT_MAX; }
  inline bool Invalid() const { return m_rows[0][0] == FLT_MAX; }
  inline void Invalidate() { m_rows[0][0] = FLT_MAX; }

  inline void GetScaled(const xp128f &w, AlignedMatrixMxNf<M, N> &_M) const {
    for (int i = 0; i < M; ++i) {
      m_rows[i].GetScaled(w, _M(i));
    }
  }

  inline void Print(const bool e = false) const {
    for (int i = 0; i < M; ++i) {
      const float *row = m_rows[i];
      for (int j = 0; j < N; ++j) {
        if (e) {
          UT::Print("%e ", row[j]);
        } else {
          UT::Print("%f ", row[j]);
        }
      }
      UT::Print("\n");
    }
  }
  inline void Print(const std::string str, const bool e) const {
    const std::string _str(str.size(), ' ');
    for (int i = 0; i < M; ++i) {
      UT::Print("%s", i == 0 ? str.c_str() : _str.c_str());
      const float *row = m_rows[i];
      for (int j = 0; j < N; ++j) {
        if (e) {
          UT::Print("%e ", row[j]);
        } else {
          UT::Print("%f ", row[j]);
        }
      }
      UT::Print("\n");
    }
  }
  inline void PrintDiagonal(const bool e = false) const {
    const int _N = std::min(M, N);
    for (int i = 0; i < _N; ++i) {
      if (e) {
        UT::Print("%e ", m_rows[i][i]);
      } else {
        UT::Print("%f ", m_rows[i][i]);
      }
    }
    UT::Print("\n");
  }
  inline void PrintDiagonal(const std::string str, const bool e) const {
    UT::Print("%s", str.c_str());
    PrintDiagonal(e);
  }

  inline bool AssertEqual(const AlignedMatrixMxNf<M, N> &_M,
                          const int verbose = 1, const std::string str = "",
                          const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    bool equal = true;
    for (int i = 0; i < M && equal; ++i) {
      equal = m_rows[i].AssertEqual(_M(i), verbose, UT::String("%s[%d]", str.c_str(), i),
                                    epsAbs, epsRel);
    }
    if (equal) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
      UT::PrintSeparator();
      _M.Print(verbose > 1);
      const AlignedMatrixMxNf<M, N> E = *this - _M;
      UT::PrintSeparator();
      E.Print(verbose > 1);
    }
    return false;
  }
  
  inline bool AssertZero(const int verbose = 1, const std::string str = "",
                         const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    bool zero = true;
    for (int i = 0; i < M && zero; ++i) {
      zero = m_rows[i].AssertZero(verbose, UT::String("%s[%d]", str.c_str(), i),
                                  epsAbs, epsRel);
    }
    if (zero) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
    }
    return false;
  }
  inline void AssertSymmetric() const {
    const int min_MN = std::min(M, N);
    for (int i = 0; i < min_MN; ++i) {
      for (int j = i; j < min_MN; ++j) {
        UT_ASSERT(m_rows[j][i] == m_rows[i][j]);
      }
    }
  }

  inline void Random(const float mMax) {
    for (int i = 0; i < M; ++i) {
      m_rows[i].Random(mMax);
    }
  }

  static inline void ApB(const AlignedMatrixMxNf<M, N> &A, const AlignedMatrixMxNf<M, N> &B,
                         AlignedMatrixMxNf<M, N> &ApB) {
    for (int i = 0; i < M; ++i)
      Row::apb(A.m_rows[i], B.m_rows[i], ApB.m_rows[i]);
  }
  static inline void AmB(const AlignedMatrixMxNf<M, N> &A, const AlignedMatrixMxNf<M, N> &B,
                         AlignedMatrixMxNf<M, N> &AmB) {
    for (int i = 0; i < M; ++i)
      Row::amb(A.m_rows[i], B.m_rows[i], AmB.m_rows[i]);
  }

 protected:
  Row m_rows[M];
};

class SIMD_ALIGN_DECLSPEC AlignedMatrix2x9f : public AlignedMatrixMxNf<2, 9> {
 public:
  inline void Set(const AlignedMatrix2x3f &M0, const AlignedMatrix2x3f &M1,
                  const AlignedMatrix2x3f &M2) {
    SetBlock(0, 0, M0);
    SetBlock(0, 3, M1);
    SetBlock(0, 6, M2);
  }
  static inline void AddABTTo(const AlignedMatrix2x9f &A, const AlignedMatrix2x9f &B,
                              SymmetricMatrix2x2f &ABT) {
    ABT.m00() += A(0).Dot(B(0));
    ABT.m01() += A(0).Dot(B(1));
    ABT.m11() += A(1).Dot(B(1));
  }
  static inline void Ab(const AlignedMatrix2x9f &A, const AlignedVector9f &b, float *Ab) {
    for (int i = 0; i < 2; ++i) {
      Ab[i] = A(i).Dot(b);
    }
  }
  static inline void AddAbTo(const AlignedMatrix2x9f &A, const AlignedVector9f &b, float *Ab) {
    for (int i = 0; i < 2; ++i) {
      Ab[i] += A(i).Dot(b);
    }
  }
  static inline void ATB(const SymmetricMatrix2x2f &A, const AlignedMatrix2x9f &B,
                         AlignedMatrix2x9f &ATB) {
    const float *B0 = B[0], *B1 = B[1];
    float *ATB0 = ATB[0], *ATB1 = ATB[1];
    for (int i = 0; i < 9; ++i) {
      ATB0[i] = A.m00() * B0[i] + A.m10() * B1[i];
      ATB1[i] = A.m01() * B0[i] + A.m11() * B1[i];
    }
  }
  static inline void AddATbTo(const AlignedMatrix2x9f &A, const Vector2f &b, AlignedVector9f &ATb) {
    const AlignedVectorNf<9> &A0 = A(0), &A1 = A(1);
    const xp128f b0 = xp128f::get(b.v0()), b1 = xp128f::get(b.v1());
    ATb.v0123() += A0.m_data4[0] * b0 + A1.m_data4[0] * b1;
    ATb.v4567() += A0.m_data4[1] * b0 + A1.m_data4[1] * b1;
    ATb.v8() += A0.m_data[8] * b.v0() + A1.m_data[8] * b.v1();
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix2x12f : public AlignedMatrixMxNf<2, 12> {
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix2x13f : public AlignedMatrixMxNf<2, 13> {
 public:
  inline AlignedMatrix2x13f() {}
  inline void Set(const Vector2f &M0, const AlignedMatrix2x6f &M1, const AlignedMatrix2x6f &M2) {
    float *r;
    r = m_rows[0];  r[0] = M0.v0();   memcpy(r + 1, M1[0], 24); memcpy(r + 7, M2[0], 24);
    r = m_rows[1];  r[0] = M0.v1();   memcpy(r + 1, M1[1], 24); memcpy(r + 7, M2[1], 24);
  }
  static inline void AB(const SymmetricMatrix2x2f &A, const AlignedMatrix2x13f &B,
                        AlignedMatrix2x13f &AB) {
    const AlignedMatrix2x13f::Row &B0 = B(0), &B1 = B(1);
    const float a00 = A.m00();
    const float a01 = A.m01();
    const float a11 = A.m11();
    AB(0).m_data4[0] = B0.m_data4[0] * a00 + B1.m_data4[0] * a01;
    AB(0).m_data4[1] = B0.m_data4[1] * a00 + B1.m_data4[1] * a01;
    AB(0).m_data4[2] = B0.m_data4[2] * a00 + B1.m_data4[2] * a01;
    AB(0).m_data[12] = B0.m_data[12] * a00 + B1.m_data[12] * a01;
    AB(1).m_data4[0] = B0.m_data4[0] * a01 + B1.m_data4[0] * a11;
    AB(1).m_data4[1] = B0.m_data4[1] * a01 + B1.m_data4[1] * a11;
    AB(1).m_data4[2] = B0.m_data4[2] * a01 + B1.m_data4[2] * a11;
    AB(1).m_data[12] = B0.m_data[12] * a01 + B1.m_data[12] * a11;
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix2x14f : public AlignedMatrixMxNf<2, 14> {
 public:
  inline AlignedMatrix2x14f() {}
  inline void Set(const AlignedMatrix2x13f &M0, const Vector2f &M1) {
    float *r;
    r = m_rows[0];  memcpy(r, M0[0], 52); r[13] = M1.v0();
    r = m_rows[1];  memcpy(r, M0[1], 52); r[13] = M1.v1();
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix3x9f : public AlignedMatrixMxNf<3, 9> {
 public:
  static inline void Ab(const AlignedMatrix3x9f &A, const AlignedVector9f &b, float *Ab) {
    for (int i = 0; i < 3; ++i) {
      Ab[i] = A(i).Dot(b);
    }
  }
  static inline void AddAbTo(const AlignedMatrix3x9f &A, const AlignedVector9f &b, float *Ab) {
    for (int i = 0; i < 3; ++i) {
      Ab[i] += A(i).Dot(b);
    }
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix3x12f : public AlignedMatrixMxNf<3, 12> {
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix3x13f : public AlignedMatrixMxNf<3, 13> {
 public:
  inline void Set(const AlignedVector3f &M0, const AlignedMatrix3x6f &M1,
                  const AlignedMatrix3x6f &M2) {
    float *r;
    r = m_rows[0];  r[0] = M0.v0();   memcpy(r + 1, &M1.m00(), 24); memcpy(r + 7, &M2.m00(), 24);
    r = m_rows[1];  r[0] = M0.v1();   memcpy(r + 1, &M1.m10(), 24); memcpy(r + 7, &M2.m10(), 24);
    r = m_rows[2];  r[0] = M0.v2();   memcpy(r + 1, &M1.m20(), 24); memcpy(r + 7, &M2.m20(), 24);
  }
  static inline void AB(const SymmetricMatrix2x2f &A0, const float A1, const AlignedMatrix3x13f &B,
                        AlignedMatrix3x13f &AB) {
    const float a00 = A0.m00(), a01 = A0.m01();
    const float a11 = A0.m11(), a22 = A1;
    const AlignedMatrix3x13f::Row &B0 = B(0), &B1 = B(1), &B2 = B(2);
    AB(0).m_data4[0] = B0.m_data4[0] * a00 + B1.m_data4[0] * a01;
    AB(0).m_data4[1] = B0.m_data4[1] * a00 + B1.m_data4[1] * a01;
    AB(0).m_data4[2] = B0.m_data4[2] * a00 + B1.m_data4[2] * a01;
    AB(0).m_data[12] = B0.m_data[12] * a00 + B1.m_data[12] * a01;
    AB(1).m_data4[0] = B0.m_data4[0] * a01 + B1.m_data4[0] * a11;
    AB(1).m_data4[1] = B0.m_data4[1] * a01 + B1.m_data4[1] * a11;
    AB(1).m_data4[2] = B0.m_data4[2] * a01 + B1.m_data4[2] * a11;
    AB(1).m_data[12] = B0.m_data[12] * a01 + B1.m_data[12] * a11;
    AB(2).m_data4[0] = B2.m_data4[0] * a22;
    AB(2).m_data4[1] = B2.m_data4[1] * a22;
    AB(2).m_data4[2] = B2.m_data4[2] * a22;
    AB(2).m_data[12] = B2.m_data[12] * a22;
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix3x14f : public AlignedMatrixMxNf<3, 14> {
 public:
  inline AlignedMatrix3x14f() {}
  inline void Set(const AlignedMatrix3x13f &M0, const Vector2f &M10, const float M11) {
    float *r;
    r = m_rows[0];  memcpy(r, M0[0], 52); r[13] = M10.v0();
    r = m_rows[1];  memcpy(r, M0[1], 52); r[13] = M10.v1();
    r = m_rows[2];  memcpy(r, M0[2], 52); r[13] = M11;
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix6x9f : public AlignedMatrixMxNf<6, 9> {
 public:
  inline void Set(const AlignedMatrix3x3f &M00, const AlignedMatrix3x3f &M01,
                  const AlignedMatrix3x3f &M02, const AlignedMatrix3x3f &M10,
                  const AlignedMatrix3x3f &M11, const AlignedMatrix3x3f &M12) {
    SetBlock(0, 0, M00);  SetBlock(0, 3, M01);  SetBlock(0, 6, M02);
    SetBlock(3, 0, M10);  SetBlock(3, 3, M11);  SetBlock(3, 6, M12);
  }
  inline void Set(const AlignedMatrix3x3f *M0, const AlignedMatrix3x3f *M1) {
    Set(M0[0], M0[1], M0[2], M1[0], M1[1], M1[2]);
  }
  inline void GetScaledColumn(const AlignedVector9f &s, AlignedMatrix6x9f &M) const {
    for (int i = 0; i < 6; ++i) {
      m_rows[i].GetScaled(s, M(i));
    }
  }
  inline void ScaleRow(const Vector6f &s) {
    for (int i = 0; i < 6; ++i) {
      m_rows[i].Scale(s[i]);
    }
  }
  inline void Increase3(const AlignedMatrix3x9f &A) {
    m_rows[3] += A(0);
    m_rows[4] += A(1);
    m_rows[5] += A(2);
  }
  template<typename TYPE>
  static inline void Ab(const AlignedMatrix6x9f &A, const AlignedVector9f &b, TYPE *Ab) {
    for (int i = 0; i < 6; ++i) {
      Ab[i] = A(i).Dot(b);
    }
  }
  template<typename TYPE>
  static inline void AddAbTo(const AlignedMatrix6x9f &A, const AlignedVector9f &b, TYPE *Ab) {
    for (int i = 0; i < 6; ++i) {
      Ab[i] += A(i).Dot(b);
    }
  }
  static inline void AddATBTo(const AlignedMatrix2x6f &A, const AlignedMatrix2x9f &B,
                              AlignedMatrix6x9f &ATB) {
    const float *A0 = A[0], *A1 = A[1], *B0 = B[0], *B1 = B[1];
    for (int i = 0; i < 6; ++i) {
      float *ATBi = ATB[i];
      for (int j = 0; j < 9; ++j) {
        ATBi[j] += A0[i] * B0[j] + A[1][i] * B1[j];
      }
    }
  }
  static inline void ATB(const AlignedMatrix6x6f &A, const AlignedMatrix6x9f &B,
                         AlignedMatrix6x9f &ATB) {
    const float *A0 = A[0], *A1 = A[1], *A2 = A[2], *A3 = A[3], *A4 = A[4], *A5 = A[5];
    const float *B0 = B[0], *B1 = B[1], *B2 = B[2], *B3 = B[3], *B4 = B[4], *B5 = B[5];
    for (int i = 0; i < 6; ++i) {
      float *ATBi = ATB[i];
      for (int j = 0; j < 9; ++j) {
        ATBi[j] = A0[i] * B0[j] + A1[i] * B1[j] + A2[i] * B2[j] +
                  A3[i] * B3[j] + A4[i] * B4[j] + A5[i] * B5[j];
      }
    }
  }
  static inline void AddATBTo(const AlignedMatrix6x6f &A, const AlignedMatrix6x9f &B,
                              AlignedMatrix6x9f &ATB) {
    const float *A0 = A[0], *A1 = A[1], *A2 = A[2], *A3 = A[3], *A4 = A[4], *A5 = A[5];
    const float *B0 = B[0], *B1 = B[1], *B2 = B[2], *B3 = B[3], *B4 = B[4], *B5 = B[5];
    for (int i = 0; i < 6; ++i) {
      float *ATBi = ATB[i];
      for (int j = 0; j < 9; ++j) {
        ATBi[j] += A0[i] * B0[j] + A1[i] * B1[j] + A2[i] * B2[j] +
                   A3[i] * B3[j] + A4[i] * B4[j] + A5[i] * B5[j];
      }
    }
  }
  static inline void AddABTToUpper(const AlignedMatrix6x9f &A, const AlignedMatrix6x9f &B,
                                   AlignedMatrix6x6f &ABT) {
    for (int i = 0; i < 6; ++i) {
      const AlignedMatrix6x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = i; j < 6; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTTo(const AlignedMatrix6x9f &A, const AlignedMatrix6x9f &B,
                              AlignedMatrix6x6f &ABT) {
    for (int i = 0; i < 6; ++i) {
      const AlignedMatrix6x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 6; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix7x8f : public AlignedMatrixMxNf<7, 8> {
 public:
  inline void Get(float &M00, Vector6f &M01, float &M02, SymmetricMatrix6x6f &M11,
                  Vector6f &M12) const {
    M00 = m_rows[0][0];
    memcpy(M01, m_rows[0] + 1, 24);
    M02 = m_rows[0][7];
    GetBlockDiagonal(1, M11);
    GetBlock(1, 7, M12);
  }
  static inline void ATBToUpper(const AlignedMatrix2x7f &A, const AlignedMatrix2x8f &B,
                                AlignedMatrix7x8f &ATB) {
    ATB(0).m_data4[0] = B.m_00_01_02_03() * A.m00() + B.m_10_11_12_13() * A.m10();
    ATB(0).m_data4[1] = B.m_04_05_06_07() * A.m00() + B.m_14_15_16_17() * A.m10();

    ATB(1).m_data4[0] = B.m_00_01_02_03() * A.m01() + B.m_10_11_12_13() * A.m11();
    ATB(1).m_data4[1] = B.m_04_05_06_07() * A.m01() + B.m_14_15_16_17() * A.m11();

    ATB(2).m_data[2] = A.m02() * B.m02() + A.m12() * B.m12();
    ATB(2).m_data[3] = A.m02() * B.m03() + A.m12() * B.m13();
    ATB(2).m_data4[1] = B.m_04_05_06_07() * A.m02() + B.m_14_15_16_17() * A.m12();

    ATB(3).m_data[3] = A.m03() * B.m03() + A.m13() * B.m13();
    ATB(3).m_data4[1] = B.m_04_05_06_07() * A.m03() + B.m_14_15_16_17() * A.m13();
    ATB(4).m_data4[1] = B.m_04_05_06_07() * A.m04() + B.m_14_15_16_17() * A.m14();

    ATB(5).m_data4[1] = B.m_04_05_06_07() * A.m05() + B.m_14_15_16_17() * A.m15();
    ATB(6).m_data[6] = A.m06() * B.m06() + A.m16() * B.m16();
    ATB(6).m_data[7] = A.m06() * B.m07() + A.m16() * B.m17();
  }
  static inline void AddATBToUpper(const AlignedMatrix2x7f &A, const AlignedMatrix2x8f &B,
                                   AlignedMatrix7x8f &ATB) {
    ATB(0).m_data4[0] += B.m_00_01_02_03() * A.m00() + B.m_10_11_12_13() * A.m10();
    ATB(0).m_data4[1] += B.m_04_05_06_07() * A.m00() + B.m_14_15_16_17() * A.m10();

    ATB(1).m_data4[0] += B.m_00_01_02_03() * A.m01() + B.m_10_11_12_13() * A.m11();
    ATB(1).m_data4[1] += B.m_04_05_06_07() * A.m01() + B.m_14_15_16_17() * A.m11();

    ATB(2).m_data[2] += A.m02() * B.m02() + A.m12() * B.m12();
    ATB(2).m_data[3] += A.m02() * B.m03() + A.m12() * B.m13();
    ATB(2).m_data4[1] += B.m_04_05_06_07() * A.m02() + B.m_14_15_16_17() * A.m12();

    ATB(3).m_data[3] += A.m03() * B.m03() + A.m13() * B.m13();
    ATB(3).m_data4[1] += B.m_04_05_06_07() * A.m03() + B.m_14_15_16_17() * A.m13();
    ATB(4).m_data4[1] += B.m_04_05_06_07() * A.m04() + B.m_14_15_16_17() * A.m14();

    ATB(5).m_data4[1] += B.m_04_05_06_07() * A.m05() + B.m_14_15_16_17() * A.m15();
    ATB(6).m_data[6] += A.m06() * B.m06() + A.m16() * B.m16();
    ATB(6).m_data[7] += A.m06() * B.m07() + A.m16() * B.m17();
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix9x6f : public AlignedMatrixMxNf<9, 6> {
 public:
  inline void Set(const AlignedMatrix3x3f &M00, const AlignedMatrix3x3f &M01,
                  const AlignedMatrix3x3f &M10, const AlignedMatrix3x3f &M11,
                  const AlignedMatrix3x3f &M20, const AlignedMatrix3x3f &M21) {
    SetBlock(0, 0, M00);  SetBlock(0, 3, M01);
    SetBlock(3, 0, M10);  SetBlock(3, 3, M11);
    SetBlock(6, 0, M20);  SetBlock(6, 3, M21);
  }
  inline void Set(const AlignedMatrix3x3f *M0, const AlignedMatrix3x3f *M1,
                  const AlignedMatrix3x3f *M2) {
    Set(M0[0], M0[1], M1[0], M1[1], M2[0], M2[1]);
  }
  inline void GetScaledColumn(const Vector6f &s, AlignedMatrix9x6f &M) const {
    AlignedVector6f _s;
    _s.Set(s);
    GetScaledColumn(_s, M);
  }
  inline void GetScaledColumn(const AlignedVector6f &s, AlignedMatrix9x6f &M) const {
    for (int i = 0; i < 9; ++i) {
      m_rows[i].GetScaled(s, M(i));
    }
  }
  inline void ScaleRow(const Vector9f &s) {
    for (int i = 0; i < 9; ++i) {
      m_rows[i].Scale(s[i]);
    }
  }
  static inline void ABT(const AlignedMatrix2x6f &A, const AlignedMatrix9x6f &B,
                         AlignedMatrix2x9f &ABT) {
    const xp128f a1 = xp128f::get(A[1]);
    for (int i = 0; i < 9; ++i) {
      const AlignedVectorNf<6> &Bi = B(i);
      ABT[0][i] = Bi.Dot(A.m_00_01_02_03(), A[0][4], A[0][5]);
      ABT[1][i] = Bi.Dot(a1, A[1][4], A[1][5]);
    }
  }
  static inline void AddABTTo(const AlignedMatrix2x6f &A, const AlignedMatrix9x6f &B,
                              AlignedMatrix2x9f &ABT) {
    const xp128f a1 = xp128f::get(A[1]);
    for (int i = 0; i < 9; ++i) {
      const AlignedVectorNf<6> &Bi = B(i);
      ABT[0][i] += Bi.Dot(A.m_00_01_02_03(), A[0][4], A[0][5]);
      ABT[1][i] += Bi.Dot(a1, A[1][4], A[1][5]);
    }
  }
  static inline void ABT(const AlignedMatrix6x6f &A, const AlignedMatrix9x6f &B,
                         AlignedMatrix6x9f &ABT) {
    const xp128f a1 = xp128f::get(A[1]), a3 = xp128f::get(A[3]), a5 = xp128f::get(A[5]);
    for (int i = 0; i < 9; ++i) { 
      const AlignedVectorNf<6> &Bi = B(i);
      ABT[0][i] = Bi.Dot(A.m_00_01_02_03(), A[0][4], A[0][5]);
      ABT[1][i] = Bi.Dot(a1, A[1][4], A[1][5]);
      ABT[2][i] = Bi.Dot(A.m_20_21_22_23(), A[2][4], A[2][5]);
      ABT[3][i] = Bi.Dot(a3, A[3][4], A[3][5]);
      ABT[4][i] = Bi.Dot(A.m_40_41_42_43(), A[4][4], A[4][5]);
      ABT[5][i] = Bi.Dot(a5, A[5][4], A[5][5]);
    }
  }
  static inline void AddABTTo(const AlignedMatrix6x6f &A, const AlignedMatrix9x6f &B,
                              AlignedMatrix6x9f &ABT) {
    const xp128f a1 = xp128f::get(A[1]), a3 = xp128f::get(A[3]), a5 = xp128f::get(A[5]);
    for (int i = 0; i < 9; ++i) { 
      const AlignedVectorNf<6> &Bi = B(i);
      ABT[0][i] += Bi.Dot(A.m_00_01_02_03(), A[0][4], A[0][5]);
      ABT[1][i] += Bi.Dot(a1, A[1][4], A[1][5]);
      ABT[2][i] += Bi.Dot(A.m_20_21_22_23(), A[2][4], A[2][5]);
      ABT[3][i] += Bi.Dot(a3, A[3][4], A[3][5]);
      ABT[4][i] += Bi.Dot(A.m_40_41_42_43(), A[4][4], A[4][5]);
      ABT[5][i] += Bi.Dot(a5, A[5][4], A[5][5]);
    }
  }
  static inline void ABT(const AlignedMatrix9x6f &A, const AlignedMatrix6x6f &B,
                         AlignedMatrix9x6f &ABT) {
    const xp128f b1 = xp128f::get(B[1]), b3 = xp128f::get(B[3]), b5 = xp128f::get(B[5]);
    for (int i = 0; i < 9; ++i) {
      const AlignedVectorNf<6> &Ai = A(i);
      ABT[i][0] = Ai.Dot(B.m_00_01_02_03(), B[0][4], B[0][5]);
      ABT[i][1] = Ai.Dot(b1, B[1][4], B[1][5]);
      ABT[i][2] = Ai.Dot(B.m_20_21_22_23(), B[2][4], B[2][5]);
      ABT[i][3] = Ai.Dot(b3, B[3][4], B[3][5]);
      ABT[i][4] = Ai.Dot(B.m_40_41_42_43(), B[4][4], B[4][5]);
      ABT[i][5] = Ai.Dot(b5, B[5][4], B[5][5]);
    }
  }
  static inline void AddABTTo(const AlignedMatrix9x6f &A, const AlignedMatrix6x6f &B,
                              AlignedMatrix9x6f &ABT) {
    const xp128f b1 = xp128f::get(B[1]), b3 = xp128f::get(B[3]), b5 = xp128f::get(B[5]);
    for (int i = 0; i < 9; ++i) {
      const AlignedVectorNf<6> &Ai = A(i);
      ABT[i][0] += Ai.Dot(B.m_00_01_02_03(), B[0][4], B[0][5]);
      ABT[i][1] += Ai.Dot(b1, B[1][4], B[1][5]);
      ABT[i][2] += Ai.Dot(B.m_20_21_22_23(), B[2][4], B[2][5]);
      ABT[i][3] += Ai.Dot(b3, B[3][4], B[3][5]);
      ABT[i][4] += Ai.Dot(B.m_40_41_42_43(), B[4][4], B[4][5]);
      ABT[i][5] += Ai.Dot(b5, B[5][4], B[5][5]);
    }
  }
  template<typename TYPE>
  static inline void Ab(const AlignedMatrix9x6f &A, const AlignedVector6f &b, TYPE *Ab) {
    for (int i = 0; i < 9; ++i) {
      Ab[i] = A(i).Dot(b);
    }
  }
  template<typename TYPE>
  static inline void AddAbTo(const AlignedMatrix9x6f &A, const AlignedVector6f &b, TYPE *Ab) {
    for (int i = 0; i < 9; ++i) {
      Ab[i] += A(i).Dot(b);
    }
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix9x9f : public AlignedMatrixMxNf<9, 9> {
 public:
  inline void Set(const Matrix9x9f &M) {
    AlignedMatrix9x9f &_M = *this;
    memcpy(_M[0], M[0], 36);
    memcpy(_M[1], M[1], 36);
    memcpy(_M[2], M[2], 36);
    memcpy(_M[3], M[3], 36);
    memcpy(_M[4], M[4], 36);
    memcpy(_M[5], M[5], 36);
    memcpy(_M[6], M[6], 36);
    memcpy(_M[7], M[7], 36);
    memcpy(_M[8], M[8], 36);
  }
  inline void Set(const AlignedMatrix3x3f &M00, const AlignedMatrix3x3f &M01,
                  const AlignedMatrix3x3f &M02, const AlignedMatrix3x3f &M10,
                  const AlignedMatrix3x3f &M11, const AlignedMatrix3x3f &M12,
                  const AlignedMatrix3x3f &M20, const AlignedMatrix3x3f &M21,
                  const AlignedMatrix3x3f &M22) {
    SetBlock(0, 0, M00);  SetBlock(0, 3, M01);  SetBlock(0, 6, M02);
    SetBlock(3, 0, M10);  SetBlock(3, 3, M11);  SetBlock(3, 6, M12);
    SetBlock(6, 0, M20);  SetBlock(6, 3, M21);  SetBlock(6, 6, M22);
  }
  inline void Set(const AlignedMatrix3x3f *M0, const AlignedMatrix3x3f *M1,
                  const AlignedMatrix3x3f *M2) {
    Set(M0[0], M0[1], M0[2], M1[0], M1[1], M1[2], M2[0], M2[1], M2[2]);
  }
  inline void Set(const Matrix9x9d &M) {
    AlignedMatrix9x9f &_M = *this;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        _M[i][j] = static_cast<float>(M[i][j]);
      }
    }
  }
  inline void Get(Matrix9x9d &M) const {
    const AlignedMatrix9x9f &_M = *this;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        M[i][j] = static_cast<float>(_M[i][j]);
      }
    }
  }
  inline void GetDiagonal(Vector9f &d) const {
    const AlignedMatrix9x9f &M = *this;
    d.v0() = M[0][0];   d.v1() = M[1][1];   d.v2() = M[2][2];
    d.v3() = M[3][3];   d.v4() = M[4][4];   d.v5() = M[5][5];
    d.v6() = M[6][6];   d.v7() = M[7][7];   d.v8() = M[8][8];
  }
  inline void GetDiagonal(Vector3f &d0, Vector3f &d1, Vector3f &d2) const {
    const AlignedMatrix9x9f &M = *this;
    d0.v0() = M[0][0];  d0.v1() = M[1][1];  d0.v2() = M[2][2];
    d1.v0() = M[3][3];  d1.v1() = M[4][4];  d1.v2() = M[5][5];
    d2.v0() = M[6][6];  d2.v1() = M[7][7];  d2.v2() = M[8][8];
  }
  inline void IncreaseDiagonal(const int i, const SymmetricMatrix3x3f &D) {
    AlignedMatrixMxNf<9, 9>::IncreaseDiagonal(i, D);
  }
  inline void IncreaseDiagonal(const float d012, const float d345, const float d678) {
    AlignedMatrix9x9f &M = *this;
    M[0][0] += d012;  M[1][1] += d012;  M[2][2] += d012;
    M[3][3] += d345;  M[4][4] += d345;  M[5][5] += d345;
    M[6][6] += d678;  M[7][7] += d678;  M[8][8] += d678;
  }
  inline void Set(const SymmetricMatrix9x9f &M) {
    AlignedMatrix9x9f &_M = *this;
    memcpy(_M[0], &M.m00(), 36);
    memcpy(_M[1] + 1, &M.m11(), 32);
    memcpy(_M[2] + 2, &M.m22(), 28);
    memcpy(_M[3] + 3, &M.m33(), 24);
    memcpy(_M[4] + 4, &M.m44(), 20);
    memcpy(_M[5] + 5, &M.m55(), 16);
    memcpy(_M[6] + 6, &M.m66(), 12);
    memcpy(_M[7] + 7, &M.m77(), 8);
    _M[8][8] = M.m88();
    _M.SetLowerFromUpper();
  }
  inline Matrix9x9f GetMatrix9x9f() const {
    Matrix9x9f M;
    const AlignedMatrix9x9f &_M = *this;
    memcpy(M[0], _M[0], 36);
    memcpy(M[1], _M[1], 36);
    memcpy(M[2], _M[2], 36);
    memcpy(M[3], _M[3], 36);
    memcpy(M[4], _M[4], 36);
    memcpy(M[5], _M[5], 36);
    memcpy(M[6], _M[6], 36);
    memcpy(M[7], _M[7], 36);
    memcpy(M[8], _M[8], 36);
    return M;
  }
  inline SymmetricMatrix9x9f GetSymmetricMatrix9x9f() const {
    SymmetricMatrix9x9f M;
    const AlignedMatrix9x9f &_M = *this;
    memcpy(&M.m00(), _M[0], 36);
    memcpy(&M.m11(), _M[1] + 1, 32);
    memcpy(&M.m22(), _M[2] + 2, 28);
    memcpy(&M.m33(), _M[3] + 3, 24);
    memcpy(&M.m44(), _M[4] + 4, 20);
    memcpy(&M.m55(), _M[5] + 5, 16);
    memcpy(&M.m66(), _M[6] + 6, 12);
    memcpy(&M.m77(), _M[7] + 7, 8);
    M.m88() = _M[8][8];
    return M;
  }
  inline void GetScaledColumn(const AlignedVector9f &s, AlignedMatrix9x9f &M) const {
    for (int i = 0; i < 9; ++i) {
      m_rows[i].GetScaled(s, M(i));
    }
  }
  inline void ScaleRow(const Vector9f &s) {
    for (int i = 0; i < 9; ++i) {
      m_rows[i].Scale(s[i]);
    }
  }
  inline bool SolveLDL(AlignedVector9f &b, const float *eps = NULL,
                       const bool decomposed = false) {
    AlignedMatrix9x9f &A = *this;
    float* _A[9] = {A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8]};
    return LS::SolveLDL(9, _A, &b.v0(), eps, decomposed);
  }
  inline int RankLDL(const float *eps = NULL) {
    AlignedMatrix9x9f &A = *this;
    float* _A[9] = {A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8]};
    return LS::RankLDL(9, _A, _A[8], eps);
  }
  inline bool InverseLDL(const float *eps = NULL, const bool decomposed = false) {
    return InverseLDL(*this, eps, decomposed);
  }
  static inline bool InverseLDL(AlignedMatrix9x9f &A, const float *eps = NULL,
                                const bool decomposed = false) {
    float* _A[9] = {A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8]};
    if (LS::InverseLDL<float>(9, _A, eps, decomposed)) {
      A.SetLowerFromUpper();
      return true;
    } else {
      A.Invalidate();
      return false;
    }
  }
  inline bool GetInverseLDL(AlignedMatrix9x9f &A, const float *eps = NULL) const {
    A = *this;
    return A.InverseLDL(eps);
  }
  inline AlignedMatrix9x9f GetInverseLDL(const float *eps = NULL) const {
    AlignedMatrix9x9f A = *this;
    A.InverseLDL(eps);
    return A;
  }
  template<typename TYPE>
  static inline void Ab(const AlignedMatrix9x9f &A, const AlignedVector9f &b, TYPE *Ab, const int i0 = 0) {
    for (int i = i0; i < 9; ++i) {
      Ab[i] = A(i).Dot(b);
    }
  }
  template<typename TYPE>
  static inline void AddAbTo(const AlignedMatrix9x9f &A, const AlignedVector9f &b, TYPE *Ab, const int i0 = 0) {
    for (int i = i0; i < 9; ++i) {
      Ab[i] += A(i).Dot(b);
    }
  }
  template<typename TYPE>
  static inline void SubtractAbFrom(const AlignedMatrix9x9f &A, const AlignedVector9f &b, TYPE *Ab) {
    for (int i = 0; i < 9; ++i) {
      Ab[i] -= A(i).Dot(b);
    }
  }
  static inline void ABT(const AlignedMatrix6x9f &A, const AlignedMatrix9x9f &B,
                         AlignedMatrix6x9f &ABT) {
    for (int i = 0; i < 6; ++i) {
      const AlignedMatrix6x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 9; ++j) {
        ABTi[j] = Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTTo(const AlignedMatrix6x9f &A, const AlignedMatrix9x9f &B,
                              AlignedMatrix6x9f &ABT) {
    for (int i = 0; i < 6; ++i) {
      const AlignedMatrix6x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 9; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void ABT(const AlignedMatrix9x6f &A, const AlignedMatrix9x6f &B,
                         AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x6f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 9; ++j) {
        ABTi[j] = Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTTo(const AlignedMatrix9x6f &A, const AlignedMatrix9x6f &B,
                              AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x6f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 9; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTToUpper(const AlignedMatrix9x6f &A, const AlignedMatrix9x6f &B,
                                   AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x6f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = i; j < 9; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void ABT(const AlignedMatrix9x9f &A, const AlignedMatrix6x9f &B,
                         AlignedMatrix9x6f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 6; ++j) {
        ABTi[j] = Ai.Dot(B(j));
      }
    }
  }
  static inline void ABT(const AlignedMatrix9x9f &A, const AlignedMatrix9x9f &B,
                         AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 9; ++j) {
        ABTi[j] = Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTTo(const AlignedMatrix9x9f &A, const AlignedMatrix9x9f &B,
                              AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 9; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTToUpper(const AlignedMatrix9x9f &A, const AlignedMatrix9x9f &B,
                                   AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = i; j < 9; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void SubtractABTFromUpper(const AlignedMatrix9x9f &A, const AlignedMatrix9x9f &B,
                                          AlignedMatrix9x9f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = i; j < 9; ++j) {
        ABTi[j] -= Ai.Dot(B(j));
      }
    }
  }
  static inline void AddABTTo(const AlignedMatrix9x9f &A, const AlignedMatrix6x9f &B,
                              AlignedMatrix9x6f &ABT) {
    for (int i = 0; i < 9; ++i) {
      const AlignedMatrix9x9f::Row &Ai = A(i);
      float *ABTi = ABT[i];
      for (int j = 0; j < 6; ++j) {
        ABTi[j] += Ai.Dot(B(j));
      }
    }
  }
  static inline void AddATBToUpper(const AlignedMatrix2x9f &A, const AlignedMatrix2x9f &B,
                                   AlignedMatrix9x9f &ATB) {
    const float *A0 = A[0], *A1 = A[1], *B0 = B[0], *B1 = B[1];
    for (int i = 0; i < 9; ++i) {
      float *ATBi = ATB[i];
      for (int j = i; j < 9; ++j) {
        ATBi[j] += A0[i] * B0[j] + A1[i] * B1[j];
      }
    }
  }
  static inline void AddATBToUpper(const AlignedMatrix6x9f &A, const AlignedMatrix6x9f &B,
                                   AlignedMatrix9x9f &ATB) {
    const float *A0 = A[0], *A1 = A[1], *A2 = A[2], *A3 = A[3], *A4 = A[4], *A5 = A[5];
    const float *B0 = B[0], *B1 = B[1], *B2 = B[2], *B3 = B[3], *B4 = B[4], *B5 = B[5];
    for (int i = 0; i < 9; ++i) {
      float *ATBi = ATB[i];
      for (int j = i; j < 9; ++j) {
        ATBi[j] += A0[i] * B0[j] + A1[i] * B1[j] + A2[i] * B2[j] +
                   A3[i] * B3[j] + A4[i] * B4[j] + A5[i] * B5[j];
      }
    }
  }
  static inline void AddATbTo(const AlignedMatrix6x9f &A, const float *b, float *ATb) {
    const float *A0 = A[0], *A1 = A[1], *A2 = A[2], *A3 = A[3], *A4 = A[4], *A5 = A[5];
    for (int i = 0; i < 9; ++i) {
      ATb[i] += A0[i] * b[0] + A1[i] * b[1] + A2[i] * b[2] +
                A3[i] * b[3] + A4[i] * b[4] + A5[i] * b[5];
    }
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix12x12f : public AlignedMatrixMxNf<12, 12> {
 public:
  inline void Set(const SymmetricMatrix6x6f &M00, const AlignedMatrix6x6f &M01,
                  const SymmetricMatrix6x6f &M11) {
    SetBlockDiagonal(0, M00);
    SetBlock(0, 6, M01);
    SetBlockDiagonal(6, M11);
    SetLowerFromUpper();
  }
  inline bool DecomposeLDL(const float *eps = NULL) {
    AlignedMatrix12x12f &A = *this;
    float *_A[12] = {A[0], A[1], A[2], A[3], A[4], A[5],
                     A[6], A[7], A[8], A[9], A[10], A[11]};
    return LS::DecomposeLDL(12, _A, eps);
  }
  inline bool SolveLDL(AlignedVector12f &b, const float *eps = NULL,
                       const bool decomposed = false) {
    AlignedMatrix12x12f &A = *this;
    float *_A[12] = {A[0], A[1], A[2], A[3], A[4], A[5],
                     A[6], A[7], A[8], A[9], A[10], A[11]};
    return LS::SolveLDL(12, _A, &b.v0(), eps, decomposed);
  }
  inline int RankLDL(const float *eps = NULL) {
    AlignedMatrix12x12f &A = *this;
    float *_A[12] = {A[0], A[1], A[2], A[3], A[4], A[5],
                     A[6], A[7], A[8], A[9], A[10], A[11]};
    return LS::RankLDL(12, _A, _A[11], eps);
  }
  static inline void Ab(const AlignedMatrix12x12f &A, AlignedVector12f &b, AlignedVector12f &Ab) {
    float *_Ab = Ab;
    for (int i = 0; i < 12; ++i) {
      _Ab[i] = A(i).Dot(b);
    }
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix13x14f : public AlignedMatrixMxNf<13, 14> {
 public:
  inline AlignedMatrix13x14f() {}
  inline void Get(float &M00, Vector6f &M01, Vector6f &M02, float &M03, SymmetricMatrix6x6f &M11,
                  AlignedMatrix6x6f &M12, Vector6f &M13,
                  SymmetricMatrix6x6f &M22, Vector6f &M23) const {
    M00 = m_rows[0][0];
    memcpy(M01, m_rows[0] + 1, 24);
    memcpy(M02, m_rows[0] + 7, 24);
    M03 = m_rows[0][13];
    GetBlockDiagonal(1, M11);
    GetBlock(1, 7, M12);
    GetBlock(1, 13, M13);
    GetBlockDiagonal(7, M22);
    GetBlock(7, 13, M23);
  }

  static inline void AddATBToUpper(const AlignedMatrix2x13f &A, const AlignedMatrix2x14f &B,
                                   AlignedMatrix13x14f &ATB) {
    float a0, a1;
    const float *A0 = A[0], *A1 = A[1];
    const AlignedMatrix2x14f::Row &B0 = B(0), &B1 = B(1);
    a0 = A0[0];
    a1 = A1[0];
    ATB(0).m_data4[0] += B0.m_data4[0] * a0 + B1.m_data4[0] * a1;
    ATB(0).m_data4[1] += B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(0).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(0).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(0).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[1];
    a1 = A1[1];
    ATB(1).m_data4[0] += B0.m_data4[0] * a0 + B1.m_data4[0] * a1;
    ATB(1).m_data4[1] += B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(1).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(1).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(1).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[2];
    a1 = A1[2];
    ATB(2).m_data[2]  += B0.m_data[2]  * a0 + B1.m_data[2]  * a1;
    ATB(2).m_data[3]  += B0.m_data[3]  * a0 + B1.m_data[3]  * a1;
    ATB(2).m_data4[1] += B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(2).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(2).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(2).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[3];
    a1 = A1[3];
    ATB(3).m_data[3]  += B0.m_data[3]  * a0 + B1.m_data[3]  * a1;
    ATB(3).m_data4[1] += B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(3).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(3).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(3).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;


    a0 = A0[4];
    a1 = A1[4];
    ATB(4).m_data4[1] += B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(4).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(4).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(4).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[5];
    a1 = A1[5];
    ATB(5).m_data4[1] += B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(5).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(5).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(5).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[6];
    a1 = A1[6];
    ATB(6).m_data[6]  += B0.m_data[6]  * a0 + B1.m_data[6]  * a1;
    ATB(6).m_data[7]  += B0.m_data[7]  * a0 + B1.m_data[7]  * a1;
    ATB(6).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(6).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(6).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[7];
    a1 = A1[7];
    ATB(7).m_data[7]  += B0.m_data[7]  * a0 + B1.m_data[7]  * a1;
    ATB(7).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(7).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(7).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[8];
    a1 = A1[8];
    ATB(8).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(8).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(8).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[9];
    a1 = A1[9];
    ATB(9).m_data4[2] += B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(9).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(9).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[10];
    a1 = A1[10];
    ATB(10).m_data[10] += B0.m_data[10] * a0 + B1.m_data[10] * a1;
    ATB(10).m_data[11] += B0.m_data[11] * a0 + B1.m_data[11] * a1;
    ATB(10).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(10).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[11];
    a1 = A1[11];
    ATB(11).m_data[11] += B0.m_data[11] * a0 + B1.m_data[11] * a1;
    ATB(11).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(11).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[12];
    a1 = A1[12];
    ATB(12).m_data[12] += B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(12).m_data[13] += B0.m_data[13] * a0 + B1.m_data[13] * a1;
  }
  static inline void ATBToUpper(const AlignedMatrix2x13f &A, const AlignedMatrix2x14f &B,
                                AlignedMatrix13x14f &ATB) {
    float a0, a1;
    const float *A0 = A[0], *A1 = A[1];
    const AlignedMatrix2x14f::Row &B0 = B(0), &B1 = B(1);
    a0 = A0[0];
    a1 = A1[0];
    ATB(0).m_data4[0] = B0.m_data4[0] * a0 + B1.m_data4[0] * a1;
    ATB(0).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(0).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(0).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(0).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[1];
    a1 = A1[1];
    ATB(1).m_data4[0] = B0.m_data4[0] * a0 + B1.m_data4[0] * a1;
    ATB(1).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(1).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(1).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(1).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[2];
    a1 = A1[2];
    ATB(2).m_data[2]  = B0.m_data[2]  * a0 + B1.m_data[2]  * a1;
    ATB(2).m_data[3]  = B0.m_data[3]  * a0 + B1.m_data[3]  * a1;
    ATB(2).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(2).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(2).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(2).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[3];
    a1 = A1[3];
    ATB(3).m_data[3]  = B0.m_data[3]  * a0 + B1.m_data[3]  * a1;
    ATB(3).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(3).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(3).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(3).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;


    a0 = A0[4];
    a1 = A1[4];
    ATB(4).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(4).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(4).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(4).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[5];
    a1 = A1[5];
    ATB(5).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1;
    ATB(5).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(5).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(5).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[6];
    a1 = A1[6];
    ATB(6).m_data[6]  = B0.m_data[6]  * a0 + B1.m_data[6]  * a1;
    ATB(6).m_data[7]  = B0.m_data[7]  * a0 + B1.m_data[7]  * a1;
    ATB(6).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(6).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(6).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[7];
    a1 = A1[7];
    ATB(7).m_data[7]  = B0.m_data[7]  * a0 + B1.m_data[7]  * a1;
    ATB(7).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(7).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(7).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[8];
    a1 = A1[8];
    ATB(8).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(8).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(8).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[9];
    a1 = A1[9];
    ATB(9).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1;
    ATB(9).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(9).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[10];
    a1 = A1[10];
    ATB(10).m_data[10] = B0.m_data[10] * a0 + B1.m_data[10] * a1;
    ATB(10).m_data[11] = B0.m_data[11] * a0 + B1.m_data[11] * a1;
    ATB(10).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(10).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[11];
    a1 = A1[11];
    ATB(11).m_data[11] = B0.m_data[11] * a0 + B1.m_data[11] * a1;
    ATB(11).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(11).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;

    a0 = A0[12];
    a1 = A1[12];
    ATB(12).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1;
    ATB(12).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1;
  }
  static inline void ATBToUpper(const AlignedMatrix3x13f &A, const AlignedMatrix3x14f &B,
                                AlignedMatrix13x14f &ATB) {
    float a0, a1, a2;
    const float *A0 = A[0], *A1 = A[1], *A2 = A[2];
    const AlignedMatrix3x14f::Row &B0 = B(0), &B1 = B(1), &B2 = B(2);
    a0 = A0[0];
    a1 = A1[0];
    a2 = A2[0];
    ATB(0).m_data4[0] = B0.m_data4[0] * a0 + B1.m_data4[0] * a1 + B2.m_data4[0] * a2;
    ATB(0).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1 + B2.m_data4[1] * a2;
    ATB(0).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(0).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(0).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[1];
    a1 = A1[1];
    a2 = A2[1];
    ATB(1).m_data4[0] = B0.m_data4[0] * a0 + B1.m_data4[0] * a1 + B2.m_data4[0] * a2;
    ATB(1).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1 + B2.m_data4[1] * a2;
    ATB(1).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(1).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(1).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[2];
    a1 = A1[2];
    a2 = A2[2];
    ATB(2).m_data[2] = B0.m_data[2] * a0 + B1.m_data[2] * a1 + B2.m_data[2] * a2;
    ATB(2).m_data[3] = B0.m_data[3] * a0 + B1.m_data[3] * a1 + B2.m_data[3] * a2;
    ATB(2).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1 + B2.m_data4[1] * a2;
    ATB(2).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(2).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(2).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[3];
    a1 = A1[3];
    a2 = A2[3];
    ATB(3).m_data[3] = B0.m_data[3] * a0 + B1.m_data[3] * a1 + B2.m_data[3] * a2;
    ATB(3).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1 + B2.m_data4[1] * a2;
    ATB(3).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(3).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(3).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[4];
    a1 = A1[4];
    a2 = A2[4];
    ATB(4).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1 + B2.m_data4[1] * a2;
    ATB(4).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(4).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(4).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[5];
    a1 = A1[5];
    a2 = A2[5];
    ATB(5).m_data4[1] = B0.m_data4[1] * a0 + B1.m_data4[1] * a1 + B2.m_data4[1] * a2;
    ATB(5).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(5).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(5).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[6];
    a1 = A1[6];
    a2 = A2[6];
    ATB(6).m_data[6] = B0.m_data[6] * a0 + B1.m_data[6] * a1 + B2.m_data[6] * a2;
    ATB(6).m_data[7] = B0.m_data[7] * a0 + B1.m_data[7] * a1 + B2.m_data[7] * a2;
    ATB(6).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(6).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(6).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[7];
    a1 = A1[7];
    a2 = A2[7];
    ATB(7).m_data[7] = B0.m_data[7] * a0 + B1.m_data[7] * a1 + B2.m_data[7] * a2;
    ATB(7).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(7).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(7).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[8];
    a1 = A1[8];
    a2 = A2[8];
    ATB(8).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(8).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(8).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[9];
    a1 = A1[9];
    a2 = A2[9];
    ATB(9).m_data4[2] = B0.m_data4[2] * a0 + B1.m_data4[2] * a1 + B2.m_data4[2] * a2;
    ATB(9).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(9).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[10];
    a1 = A1[10];
    a2 = A2[10];
    ATB(10).m_data[10] = B0.m_data[10] * a0 + B1.m_data[10] * a1 + B2.m_data[10] * a2;
    ATB(10).m_data[11] = B0.m_data[11] * a0 + B1.m_data[11] * a1 + B2.m_data[11] * a2;
    ATB(10).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(10).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[11];
    a1 = A1[11];
    a2 = A2[11];
    ATB(11).m_data[11] = B0.m_data[11] * a0 + B1.m_data[11] * a1 + B2.m_data[11] * a2;
    ATB(11).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(11).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;

    a0 = A0[12];
    a1 = A1[12];
    a2 = A2[12];
    ATB(12).m_data[12] = B0.m_data[12] * a0 + B1.m_data[12] * a1 + B2.m_data[12] * a2;
    ATB(12).m_data[13] = B0.m_data[13] * a0 + B1.m_data[13] * a1 + B2.m_data[13] * a2;
  }

  static inline void AddabTToUpper(const AlignedVector13f &a, const AlignedVector14f &b,
                                   AlignedMatrix13x14f &abT) {
    AlignedVector14f::AddsATo(a[0], b, abT(0));
    AlignedVector14f::AddsATo(a[1], b, abT(1));
    AlignedVector14f::AddsATo(a[2], b, abT(2));
    AlignedVector14f::AddsATo(a[3], b, abT(3));
    AlignedVector14f::AddsATo(a[4], b, abT(4));
    AlignedVector14f::AddsATo(a[5], b, abT(5));
    AlignedVector14f::AddsATo(a[6], b, abT(6));
    AlignedVector14f::AddsATo(a[7], b, abT(7));
    AlignedVector14f::AddsATo(a[8], b, abT(8));
    AlignedVector14f::AddsATo(a[9], b, abT(9));
    AlignedVector14f::AddsATo(a[10], b, abT(10));
    AlignedVector14f::AddsATo(a[11], b, abT(11));
    AlignedVector14f::AddsATo(a[12], b, abT(12));
  }
};

template<int M, int N>
class MatrixMxNf {

 public:

  inline const float* operator[] (const int i) const { return m_data[i]; }
  inline       float* operator[] (const int i)       { return m_data[i]; }

  inline MatrixMxNf<M, N> operator - (const MatrixMxNf<M, N> &B) const {
    MatrixMxNf<M, N> AmB;
    const MatrixMxNf<M, N> &A = *this;
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        AmB[i][j] = A[i][j] - B[i][j];
      }
    }
    return AmB;
  }

  inline void SetBlockDiagonal(const int i, const SymmetricMatrix3x3f &B) {
    int k = i;
    memcpy(m_data[k] + k, &B.m00(), 12);  ++k;
    memcpy(m_data[k] + k, &B.m11(), 8); ++k;
    m_data[k][k] = B.m22();
  }
  inline void SetBlockDiagonal(const int i, const SymmetricMatrix8x8f &B) {
    int k = i;
    memcpy(m_data[k] + k, &B.m00(), 32);  ++k;
    memcpy(m_data[k] + k, &B.m11(), 28);  ++k;
    memcpy(m_data[k] + k, &B.m22(), 24);  ++k;
    memcpy(m_data[k] + k, &B.m33(), 20);  ++k;
    memcpy(m_data[k] + k, &B.m44(), 16);  ++k;
    memcpy(m_data[k] + k, &B.m55(), 12);  ++k;
    memcpy(m_data[k] + k, &B.m66(), 8);   ++k;
    m_data[k][k] = B.m77();
  }

  inline void GetBlockDiagonal(const int i, SymmetricMatrix6x6f &B) const {
    int k = i;
    memcpy(&B.m00(), m_data[k] + k, 24);  ++k;
    memcpy(&B.m11(), m_data[k] + k, 20);  ++k;
    memcpy(&B.m22(), m_data[k] + k, 16);  ++k;
    memcpy(&B.m33(), m_data[k] + k, 12);  ++k;
    memcpy(&B.m44(), m_data[k] + k, 8);   ++k;
    B.m55() = m_data[k][k];
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix9x9f &B) const {
    int k = i;
    memcpy(&B.m00(), m_data[k] + k, 36);  ++k;
    memcpy(&B.m11(), m_data[k] + k, 32);  ++k;
    memcpy(&B.m22(), m_data[k] + k, 28);  ++k;
    memcpy(&B.m33(), m_data[k] + k, 24);  ++k;
    memcpy(&B.m44(), m_data[k] + k, 20);  ++k;
    memcpy(&B.m55(), m_data[k] + k, 16);  ++k;
    memcpy(&B.m66(), m_data[k] + k, 12);  ++k;
    memcpy(&B.m77(), m_data[k] + k, 8);   ++k;
    B.m88() = m_data[k][k];
  }

  inline void SetBlock(const int i, const int j, const Matrix3x3f &B) {
    memcpy(m_data[i] + j, &B.m00(), 12);
    memcpy(m_data[i + 1] + j, &B.m10(), 12);
    memcpy(m_data[i + 2] + j, &B.m20(), 12);
  }
  inline void SetBlock(const int i, const int j, const Vector3f &v) {
    m_data[i][j] = v.v0();
    m_data[i + 1][j] = v.v1();
    m_data[i + 2][j] = v.v2();
  }
  inline void SetBlock(const int i, const int j, const AlignedVector8f &v) {
    m_data[i][j] = v.v0();
    m_data[i + 1][j] = v.v1();
    m_data[i + 2][j] = v.v2();
    m_data[i + 3][j] = v.v3();
    m_data[i + 4][j] = v.v4();
    m_data[i + 5][j] = v.v5();
    m_data[i + 6][j] = v.v6();
    m_data[i + 7][j] = v.v7();
  }

  inline void GetBlock(const int i, const int j, AlignedMatrix3x3f &B) const {
    memcpy(&B.m00(), m_data[i] + j, 12);
    memcpy(&B.m10(), m_data[i + 1] + j, 12);
    memcpy(&B.m20(), m_data[i + 2] + j, 12);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix4x4f &B) const {
    memcpy(&B.m00(), m_data[i] + j, 16);
    memcpy(&B.m10(), m_data[i + 1] + j, 16);
    memcpy(&B.m20(), m_data[i + 2] + j, 16);
    memcpy(&B.m30(), m_data[i + 3] + j, 16);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix6x6f &B) const {
    memcpy(B[0], m_data[i] + j, 24);
    memcpy(B[1], m_data[i + 1] + j, 24);
    memcpy(B[2], m_data[i + 2] + j, 24);
    memcpy(B[3], m_data[i + 3] + j, 24);
    memcpy(B[4], m_data[i + 4] + j, 24);
    memcpy(B[5], m_data[i + 5] + j, 24);
  }
  inline void GetBlock(const int i, const int j, AlignedVector3f &v) const {
    v.v0() = m_data[i][j];
    v.v1() = m_data[i + 1][j];
    v.v2() = m_data[i + 2][j];
  }
  inline void GetBlock(const int i, const int j, Vector6f &v) const {
    v.v0() = m_data[i][j];
    v.v1() = m_data[i + 1][j];
    v.v2() = m_data[i + 2][j];
    v.v3() = m_data[i + 3][j];
    v.v4() = m_data[i + 4][j];
    v.v5() = m_data[i + 5][j];
  }
  inline void GetBlock(const int i, const int j, Vector9f &v) const {
    v.v0() = m_data[i][j];
    v.v1() = m_data[i + 1][j];
    v.v2() = m_data[i + 2][j];
    v.v3() = m_data[i + 3][j];
    v.v4() = m_data[i + 4][j];
    v.v5() = m_data[i + 5][j];
    v.v6() = m_data[i + 6][j];
    v.v7() = m_data[i + 7][j];
    v.v8() = m_data[i + 8][j];
  }

  inline void MakeZero() { memset(this, 0, sizeof(MatrixMxNf<M, N>)); }

  inline void Print(const bool e = false) const {
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        if (e) {
          UT::Print("%e ", m_data[i][j]);
        } else {
          UT::Print("%.4f ", m_data[i][j]);
        }
      }
      UT::Print("\n");
    }
  }

  inline bool AssertEqual(const MatrixMxNf<M, N> &_M, const int verbose = 1, const std::string str = "",
                          const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    if (UT::VectorAssertEqual(&m_data[0][0], &_M.m_data[0][0], M * N, verbose, str,
                              epsAbs, epsRel)) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
      UT::PrintSeparator();
      _M.Print(verbose > 1);
      const MatrixMxNf<M, N> E = *this - _M;
      UT::PrintSeparator();
      E.Print(verbose > 1);
    }
    return false;
  }

 protected:
  float m_data[M][N];
};

class Matrix9x10f : public MatrixMxNf<9, 10> {
};
class Matrix12x13f : public MatrixMxNf<12, 13> {
};
class Matrix18x19f : public MatrixMxNf<18, 19> {
};

template<int N>
class SIMD_ALIGN_DECLSPEC AlignedMatrixNx3f : public AlignedMatrixMxNf<N, 3> {
 public:
  inline void SetBlock(const int i, const AlignedMatrix3x3f &M) {
    memcpy(this->m_rows[i].m_data, &M.m00(), 12);
    memcpy(this->m_rows[i + 1].m_data, &M.m10(), 12);
    memcpy(this->m_rows[i + 2].m_data, &M.m20(), 12);
  }

  static inline void ABT(const AlignedMatrixNx3f<N> &A, const AlignedMatrix3x3f &B,
                         AlignedMatrixNx3f<N> &ABT) {
    for (int i = 0; i < N; ++i) {
      const xp128f &a = A(i).m_data4[0];
      ABT[i][0] = (B.m_00_01_02_r0() * a).vsum_012();
      ABT[i][1] = (B.m_10_11_12_r1() * a).vsum_012();
      ABT[i][2] = (B.m_20_21_22_r2() * a).vsum_012();
    }
  }
  static inline void ABTToUpper(const AlignedMatrixNx3f<N> &A, const AlignedMatrixNx3f<N> &B0,
                                const AlignedVector3f &B1, MatrixMxNf < N, N + 1 > &ABT) {
    for (int i = 0; i < N; ++i) {
      const xp128f &a = A(i).m_data4[0];
      for (int j = i; j < N; ++j)
        ABT[i][j] = (B0(j).m_data4[0] * a).vsum_012();
      ABT[i][N] = (B1.v012r() * a).vsum_012();
    }
  }

};

class SIMD_ALIGN_DECLSPEC AlignedMatrix9x3f : public AlignedMatrixNx3f<9> {
 public:
  inline void Set(const AlignedMatrix3x3f &M0, const AlignedMatrix3x3f &M1,
                  const AlignedMatrix3x3f &M2) {
    SetBlock(0, M0);
    SetBlock(3, M1);
    SetBlock(6, M2);
  }
  static inline void Ab(const AlignedMatrix9x3f &A, const AlignedVector3f &b, float *Ab) {
    for (int i = 0; i < 9; ++i) {
      Ab[i] = A(i).Dot(b);
    }
  }
  static inline void AddAbTo(const AlignedMatrix9x3f &A, const AlignedVector3f &b, float *Ab) {
    for (int i = 0; i < 9; ++i) {
      Ab[i] += A(i).Dot(b);
    }
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix12x3f : public AlignedMatrixNx3f<12> {
 public:
  inline void Set(const AlignedMatrix3x3f &M0, const AlignedMatrix3x3f &M1,
                  const AlignedMatrix3x3f &M2, const AlignedMatrix3x3f &M3) {
    SetBlock(0, M0);
    SetBlock(3, M1);
    SetBlock(6, M2);
    SetBlock(9, M3);
  }
};
class SIMD_ALIGN_DECLSPEC AlignedMatrix18x3f : public AlignedMatrixNx3f<18> {
 public:
  inline void Set(const AlignedMatrix3x3f &M0, const AlignedMatrix3x3f &M1,
                  const AlignedMatrix3x3f &M2,
                  const AlignedMatrix3x3f &M3, const AlignedMatrix3x3f &M4, const AlignedMatrix3x3f &M5) {
    SetBlock(0, M0);
    SetBlock(3, M1);
    SetBlock(6, M2);
    SetBlock(9, M3);
    SetBlock(12, M4);
    SetBlock(15, M5);
  }
};

class AlignedMatrixXf : public AlignedMatrixX<float> {

 public:

  inline AlignedMatrixXf() {}
  inline AlignedMatrixXf(const AlignedMatrixXf &M) { *this = M; }
  inline void operator = (const AlignedMatrixXf &M) { Set(M); }

  inline void operator *= (const float s) {
    const xp128f _s = xp128f::get(s);
    Scale(_s);
  }
  inline void operator *= (const xp128f &s) { Scale(s); }
  inline AlignedMatrixXf operator - (const AlignedMatrixXf &B) const {
    const AlignedMatrixXf &A = *this;
    AlignedMatrixXf AmB;
    AmB.Resize(m_Nr, m_Nc, m_symmetric);
    if (m_symmetric) {
      SIMD::Subtract((m_Nc * (m_Nc + 1)) >> 1, A.Data(), B.Data(), AmB.Data());
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        SIMD::Subtract(m_Nc, A[i], B[i], AmB[i]);
      }
    }
    return AmB;
  }

  inline bool Valid() const { return m_rows[0][0] != FLT_MAX; }
  inline bool Invalid() const { return m_rows[0][0] == FLT_MAX; }
  inline void Invalidate() { m_rows[0][0] = FLT_MAX; }

  inline void SetLowerFromUpper() {
    const int N = std::min(m_Nr, m_Nc);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < i; ++j) {
        m_rows[i][j] = m_rows[j][i];
      }
    }
  }

  inline void SetBlockDiagonal(const int i, const int N, const float m) {
    const int i2 = i + N;
    for (int _i = i; _i < i2; ++_i) {
      m_rows[_i][_i] = m;
    }
  }
  inline void SetBlockDiagonal(const int i, const SymmetricMatrix2x2f &M) {
    int k = i;
    memcpy(m_rows[k] + k, &M.m00(), 8);   ++k;
    m_rows[k][k] = M.m11();
  }
  inline void SetBlockDiagonal(const int i, const SymmetricMatrix6x6f &M) {
    int k = i;
    memcpy(m_rows[k] + k, &M.m00(), 24);  ++k;
    memcpy(m_rows[k] + k, &M.m11(), 20);  ++k;
    memcpy(m_rows[k] + k, &M.m22(), 16);  ++k;
    memcpy(m_rows[k] + k, &M.m33(), 12);  ++k;
    memcpy(m_rows[k] + k, &M.m44(), 8);   ++k;
    m_rows[k][k] = M.m55();
  }

  inline void SetBlock(const int i, const int j, const Matrix2x3f &M) {
    memcpy(m_rows[i] + j, M[0], 12);
    memcpy(m_rows[i + 1] + j, M[1], 12);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix2x6f &M) {
    memcpy(m_rows[i] + j, M[0], 24);
    memcpy(m_rows[i + 1] + j, M[1], 24);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix2x9f &M) {
    memcpy(m_rows[i] + j, M[0], 36);
    memcpy(m_rows[i + 1] + j, M[1], 36);
  }
  inline void SetBlock(const int i, const int j, const Matrix3x2f &M) {
    memcpy(m_rows[i] + j, M[0], 8);
    memcpy(m_rows[i + 1] + j, M[1], 8);
    memcpy(m_rows[i + 2] + j, M[2], 8);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix3x3f &M) {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(m_rows[i] + k++, &M.m00(), 12);
      memcpy(m_rows[i + 1] + k++, &M.m11(), 8);
      m_rows[k][k] = M.m22();
    } else {
      memcpy(m_rows[i] + j, &M.m00(), 12);
      memcpy(m_rows[i + 1] + j, &M.m10(), 12);
      memcpy(m_rows[i + 2] + j, &M.m20(), 12);
    }
  }
  inline void SetBlock(const int i, const int j, const Matrix3x3f &M) {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(m_rows[i] + k++, M[0], 12);
      memcpy(m_rows[i + 1] + k++, M[1] + 1, 8);
      m_rows[k][k] = M[2][2];
    } else {
      memcpy(m_rows[i] + j, M[0], 12);
      memcpy(m_rows[i + 1] + j, M[1], 12);
      memcpy(m_rows[i + 2] + j, M[2], 12);
    }
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix6x6f &M) {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(m_rows[i] + k++, &M[0][0], 24);
      memcpy(m_rows[i + 1] + k++, &M[1][1], 20);
      memcpy(m_rows[i + 2] + k++, &M[2][2], 16);
      memcpy(m_rows[i + 3] + k++, &M[3][3], 12);
      memcpy(m_rows[i + 4] + k++, &M[4][4], 8);
      m_rows[k][k] = M[5][5];
    } else {
      memcpy(m_rows[i] + j, M[0], 24);
      memcpy(m_rows[i + 1] + j, M[1], 24);
      memcpy(m_rows[i + 2] + j, M[2], 24);
      memcpy(m_rows[i + 3] + j, M[3], 24);
      memcpy(m_rows[i + 4] + j, M[4], 24);
      memcpy(m_rows[i + 5] + j, M[5], 24);
    }
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix6x9f &M) {
    memcpy(m_rows[i] + j, M[0], 36);
    memcpy(m_rows[i + 1] + j, M[1], 36);
    memcpy(m_rows[i + 2] + j, M[2], 36);
    memcpy(m_rows[i + 3] + j, M[3], 36);
    memcpy(m_rows[i + 4] + j, M[4], 36);
    memcpy(m_rows[i + 5] + j, M[5], 36);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix9x3f &M) {
    memcpy(m_rows[i] + j, M[0], 12);
    memcpy(m_rows[i + 1] + j, M[1], 12);
    memcpy(m_rows[i + 2] + j, M[2], 12);
    memcpy(m_rows[i + 3] + j, M[3], 12);
    memcpy(m_rows[i + 4] + j, M[4], 12);
    memcpy(m_rows[i + 5] + j, M[5], 12);
    memcpy(m_rows[i + 6] + j, M[6], 12);
    memcpy(m_rows[i + 7] + j, M[7], 12);
    memcpy(m_rows[i + 8] + j, M[8], 12);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix9x6f &M) {
    memcpy(m_rows[i] + j, M[0], 24);
    memcpy(m_rows[i + 1] + j, M[1], 24);
    memcpy(m_rows[i + 2] + j, M[2], 24);
    memcpy(m_rows[i + 3] + j, M[3], 24);
    memcpy(m_rows[i + 4] + j, M[4], 24);
    memcpy(m_rows[i + 5] + j, M[5], 24);
    memcpy(m_rows[i + 6] + j, M[6], 24);
    memcpy(m_rows[i + 7] + j, M[7], 24);
    memcpy(m_rows[i + 8] + j, M[8], 24);
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrix9x9f &M) {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(m_rows[i] + k++, &M[0][0], 36);
      memcpy(m_rows[i + 1] + k++, &M[1][1], 32);
      memcpy(m_rows[i + 2] + k++, &M[2][2], 28);
      memcpy(m_rows[i + 3] + k++, &M[3][2], 24);
      memcpy(m_rows[i + 4] + k++, &M[4][3], 20);
      memcpy(m_rows[i + 5] + k++, &M[5][4], 16);
      memcpy(m_rows[i + 6] + k++, &M[6][5], 12);
      memcpy(m_rows[i + 7] + k++, &M[7][6], 8);
      m_rows[k][k] = M[8][8];
    } else {
      memcpy(m_rows[i] + j, M[0], 36);
      memcpy(m_rows[i + 1] + j, M[1], 36);
      memcpy(m_rows[i + 2] + j, M[2], 36);
      memcpy(m_rows[i + 3] + j, M[3], 36);
      memcpy(m_rows[i + 4] + j, M[4], 36);
      memcpy(m_rows[i + 5] + j, M[5], 36);
      memcpy(m_rows[i + 6] + j, M[6], 36);
      memcpy(m_rows[i + 7] + j, M[7], 36);
      memcpy(m_rows[i + 8] + j, M[8], 36);
    }
  }
  inline void SetBlock(const int i, const int j, const AlignedMatrixXf &M) {
    if (m_symmetric && i == j) {
      size_t size = sizeof(float) * M.m_Nc;
      for (int k = 0; k < M.m_Nr; ++k, size -= sizeof(float)) {
        memcpy(m_rows[i + k] + j + k, M[k] + k, size);
      }
    } else {
      const size_t size = sizeof(float) * M.m_Nc;
      for (int k = 0; k < M.m_Nr; ++k) {
        memcpy(m_rows[i + k] + j, M[k], size);
      }
    }
  }
  inline void SetBlock(const int i, const int j, const Vector6f &v) {
    m_rows[i][j] = v.v0();
    m_rows[i + 1][j] = v.v1();
    m_rows[i + 2][j] = v.v2();
    m_rows[i + 3][j] = v.v3();
    m_rows[i + 4][j] = v.v4();
    m_rows[i + 5][j] = v.v5();
  }
  inline void SetBlock(const int i, const int j, const Vector9f &v) {
    m_rows[i][j] = v.v0();
    m_rows[i + 1][j] = v.v1();
    m_rows[i + 2][j] = v.v2();
    m_rows[i + 3][j] = v.v3();
    m_rows[i + 4][j] = v.v4();
    m_rows[i + 5][j] = v.v5();
    m_rows[i + 6][j] = v.v6();
    m_rows[i + 7][j] = v.v7();
    m_rows[i + 8][j] = v.v8();
  }
  inline void SetBlock(const int i, const int j, const AlignedVectorXf &V) {
    const int N = V.Size();
    for (int k = 0; k < N; ++k) {
      m_rows[i + k][j] = V[k];
    }
  }
  
  inline void GetBlockDiagonal(const int i, SymmetricMatrix2x2f &M) const {
    int k = i;
    memcpy(&M.m00(), m_rows[k] + k, 8);   ++k;
    M.m11() = m_rows[k][k];
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix3x3f &M) const {
    int k = i;
    memcpy(&M.m00(), m_rows[k] + k, 12);  ++k;
    memcpy(&M.m11(), m_rows[k] + k, 8);   ++k;
    M.m22() = m_rows[k][k];
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix6x6f &M) const {
    int k = i;
    memcpy(&M.m00(), m_rows[k] + k, 24);  ++k;
    memcpy(&M.m11(), m_rows[k] + k, 20);  ++k;
    memcpy(&M.m22(), m_rows[k] + k, 16);  ++k;
    memcpy(&M.m33(), m_rows[k] + k, 12);  ++k;
    memcpy(&M.m44(), m_rows[k] + k, 8);   ++k;
    M.m55() = m_rows[k][k];
  }
  inline void GetBlock(const int i, const int j, Matrix2x3f &M) const {
    memcpy(M[0], m_rows[i] + j, 12);
    memcpy(M[1], m_rows[i + 1] + j, 12);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix2x6f &M) const {
    memcpy(M[0], m_rows[i] + j, 24);
    memcpy(M[1], m_rows[i + 1] + j, 24);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix2x9f &M) const {
    memcpy(M[0], m_rows[i] + j, 36);
    memcpy(M[1], m_rows[i + 1] + j, 36);
  }
  inline void GetBlock(const int i, const int j, Matrix3x2f &M) const {
    memcpy(M[0], m_rows[i] + j, 8);
    memcpy(M[1], m_rows[i + 1] + j, 8);
    memcpy(M[2], m_rows[i + 2] + j, 8);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix3x3f &M) const {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(&M.m00(), m_rows[i] + k++, 12);
      memcpy(&M.m11(), m_rows[i + 1] + k++, 8);
      M.m22() = m_rows[k][k];
      M.SetLowerFromUpper();
    } else {
      memcpy(&M.m00(), m_rows[i] + j, 12);
      memcpy(&M.m10(), m_rows[i + 1] + j, 12);
      memcpy(&M.m20(), m_rows[i + 2] + j, 12);
    }
  }
  inline void GetBlock(const int i, const int j, Matrix3x3f &M) const {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(M[0], m_rows[i] + k++, 12);
      memcpy(M[1] + 1, m_rows[i + 1] + k++, 8);
      M[2][2] = m_rows[k][k];
      M.SetLowerFromUpper();
    } else {
      memcpy(M[0], m_rows[i] + j, 12);
      memcpy(M[1], m_rows[i + 1] + j, 12);
      memcpy(M[2], m_rows[i + 2] + j, 12);
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix6x6f &M) const {
    if (m_symmetric && i == j) {
      int k = j;
      memcpy(&M[0][0], m_rows[i] + k++, 24);
      memcpy(&M[1][1], m_rows[i + 1] + k++, 20);
      memcpy(&M[2][2], m_rows[i + 2] + k++, 16);
      memcpy(&M[3][3], m_rows[i + 3] + k++, 12);
      memcpy(&M[4][4], m_rows[i + 4] + k++, 8);
      M[5][5] = m_rows[k][k];
      M.SetLowerFromUpper();
    } else {
      memcpy(M[0], m_rows[i] + j, 24);
      memcpy(M[1], m_rows[i + 1] + j, 24);
      memcpy(M[2], m_rows[i + 2] + j, 24);
      memcpy(M[3], m_rows[i + 3] + j, 24);
      memcpy(M[4], m_rows[i + 4] + j, 24);
      memcpy(M[5], m_rows[i + 5] + j, 24);
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix6x9f &M) const {
    memcpy(M[0], m_rows[i] + j, 36);
    memcpy(M[1], m_rows[i + 1] + j, 36);
    memcpy(M[2], m_rows[i + 2] + j, 36);
    memcpy(M[3], m_rows[i + 3] + j, 36);
    memcpy(M[4], m_rows[i + 4] + j, 36);
    memcpy(M[5], m_rows[i + 5] + j, 36);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix9x9f &M) const {
    memcpy(M[0], m_rows[i] + j, 36);
    memcpy(M[1], m_rows[i + 1] + j, 36);
    memcpy(M[2], m_rows[i + 2] + j, 36);
    memcpy(M[3], m_rows[i + 3] + j, 36);
    memcpy(M[4], m_rows[i + 4] + j, 36);
    memcpy(M[5], m_rows[i + 5] + j, 36);
    memcpy(M[6], m_rows[i + 6] + j, 36);
    memcpy(M[7], m_rows[i + 7] + j, 36);
    memcpy(M[8], m_rows[i + 8] + j, 36);
  }
  inline void GetBlock(const int i, const int j, AlignedMatrixXf &M) const {
    if (m_symmetric && i == j) {
      size_t size = sizeof(float) * M.m_Nc;
      for (int k = 0; k < M.m_Nr; ++k, size -= sizeof(float)) {
        memcpy(M[k] + k, m_rows[i + k] + j + k, size);
      }
    } else {
      const size_t size = sizeof(float) * M.m_Nc;
      for (int k = 0; k < M.m_Nr; ++k) {
        memcpy(M[k], m_rows[i + k] + j, size);
      }
    }
  }
  inline void GetBlock(const int i, const int j, Vector6f &v) const {
    v.v0() = m_rows[i][j];
    v.v1() = m_rows[i + 1][j];
    v.v2() = m_rows[i + 2][j];
    v.v3() = m_rows[i + 3][j];
    v.v4() = m_rows[i + 4][j];
    v.v5() = m_rows[i + 5][j];
  }
  inline void GetBlock(const int i, const int j, AlignedVectorXf &V) const {
    const int N = V.Size();
    for (int k = 0; k < N; ++k) {
      V[k] = m_rows[i + k][j];
    }
  }
  
  inline void GetDiagonal(const int i, Vector2f &v) const {
    v.v0() = m_rows[i][i];
    v.v1() = m_rows[i + 1][i + 1];
  }
  inline void GetDiagonal(const int i, AlignedVector3f &v) const {
    v.v0() = m_rows[i][i];
    v.v1() = m_rows[i + 1][i + 1];
    v.v2() = m_rows[i + 2][i + 2];
  }
  inline void GetDiagonal(const int i, Vector3f &v) const {
    v.v0() = m_rows[i][i];
    v.v1() = m_rows[i + 1][i + 1];
    v.v2() = m_rows[i + 2][i + 2];
  }

  inline void IncreaseBlockDiagonal(const int i, const SymmetricMatrix2x2f &M) {
    int k = i;
    m_rows[k][k] += M.m00();
    m_rows[k][k + 1] += M.m01();
    ++k;
    m_rows[k][k] += M.m11();
  }
  inline void IncreaseBlockDiagonal(const int i, const SymmetricMatrix6x6f &M) {
    const float *_M = M;
    for (int _i = 0, k = 0; _i < 6; ++_i) {
      for (int _j = _i; _j < 6; ++_j, ++k) {
        m_rows[i + _i][i + _j] += _M[k];
      }
    }
  }

  inline void IncreaseBlock(const int i, const int j, const AlignedMatrix2x3f &M) {
    const float *_M[2] = {&M.m00(), &M.m10()};
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] += _M[_i][_j];
      }
    }
  }
  inline void IncreaseBlock(const int i, const int j, const Matrix3x2f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = 0; _j < 2; ++_j) {
        m_rows[i + _i][j + _j] += M[_i][_j];
      }
    }
  }
  inline void IncreaseBlock(const int i, const int j, const AlignedMatrix3x3f &M) {
    const float *_M[3] = {&M.m00(), &M.m10(), &M.m20()};
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] += _M[_i][_j];
      }
    }
  }
  inline void IncreaseBlock(const int i, const int j, const AlignedMatrix6x6f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 6; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 6; ++_j) {
        m_rows[i + _i][j + _j] += M[_i][_j];
      }
    }
  }

  inline void Scale(const xp128f &s) {
    if (m_symmetric) {
      SIMD::Multiply<0>((m_Nc * (m_Nc + 1)) >> 1, s, Data());
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        SIMD::Multiply<0>(m_Nc, s, m_rows[i]);
      }
    }
  }

  inline bool SolveLDL(AlignedVectorXf &b, const float *eps = NULL,
                       const bool decomposed = false) {
    return LS::SolveLDL(m_Nr, m_rows.data(), b.Data(), eps, decomposed);
  }
  inline int RankLDL(AlignedVectorXf *work, const float *eps = NULL) {
    work->Resize(m_Nc);
    return LS::RankLDL(m_Nr, m_rows.data(), work->Data(), eps);
  }
  inline bool InverseLDL(const float *eps = NULL, const bool decomposed = false) {
    if (LS::InverseLDL(m_Nr, m_rows.data(), eps, decomposed)) {
      SetLowerFromUpper();
      return true;
    } else {
      Invalidate();
      return false;
    }
  }
  inline bool MarginalizeLDL(const int i, const int N, AlignedVectorXf &b,
                             AlignedVectorXf *work, const float *eps = NULL,
                             const bool erase = true) {
    work->Resize(m_Nc);
    const bool scc = LS::MarginalizeLDL(m_Nr, m_rows.data(), b.Data(), i, N, work->Data(), eps);
    if (erase) {
      Erase(i, N);
      b.Erase(i, N);
    } else {
      MakeZero(i, N);
      b.MakeZero(i, N);
    }
    return scc;
  }

  inline void Print(const bool e = false) const {
    if (m_symmetric) {
      for (int i = 0; i < m_Nr; ++i) {
        for (int j = 0; j < i; ++j) {
          if (e) {
            UT::Print("%e ", m_rows[j][i]);
          } else {
            UT::Print("%f ", m_rows[j][i]);
          }
        }
        for (int j = i; j < m_Nc; ++j) {
          if (e) {
            UT::Print("%e ", m_rows[i][j]);
          } else {
            UT::Print("%f ", m_rows[i][j]);
          }
        }
        UT::Print("\n");
      }
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        for (int j = 0; j < m_Nc; ++j) {
          if (e) {
            UT::Print("%e ", m_rows[i][j]);
          } else {
            UT::Print("%f ", m_rows[i][j]);
          }
        }
        UT::Print("\n");
      }
    }
  }
  inline void PrintRow(const int i, const bool e = false) const {
    for (int j = 0; j < m_Nc; ++j) {
      if (e) {
        UT::Print("%e ", m_rows[i][j]);
      } else {
        UT::Print("%f ", m_rows[i][j]);
      }
    }
    UT::Print("\n");
  }
  inline void PrintDiagonal(const bool e = false) const {
    const int N = std::min(m_Nr, m_Nc);
    for (int i = 0; i < N; ++i) {
      if (e) {
        UT::Print("%e ", m_rows[i][i]);
      } else {
        UT::Print("%f ", m_rows[i][i]);
      }
    }
    UT::Print("\n");
  }
  inline void PrintBlock(const int i, const int j, const int Nr, const int Nc,
                         const bool e = false) const {
    for (int _i = i; _i < Nr; ++_i) {
      for (int _j = j; _j < Nc; ++_j) {
        if (e) {
          UT::Print("%e ", m_rows[_i][_j]);
        } else {
          UT::Print("%f ", m_rows[_i][_j]);
        }
      }
      UT::Print("\n");
    }
  }

  inline void Random(const float mMax) {
    if (m_symmetric) {
      for (int i = 0; i < m_Nr; ++i) {
        UT::Random(m_rows[i] + i, m_Nc - i, mMax);
      }
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        UT::Random(m_rows[i], m_Nc, mMax);
      }
    }
  }

  inline bool AssertEqual(const AlignedMatrixXf &M,
                          const int verbose = 1, const std::string str = "",
                          const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    AlignedVectorXf V1, V2;
    bool equal = M.m_Nr == m_Nr && M.m_Nc == m_Nc && M.m_symmetric == m_symmetric;
    //const int verbose1 = verbose & 15, verbose2 = verbose >> 4;
    const int verbose1 = verbose, verbose2 = verbose;
    for (int i = 0; i < m_Nr && equal; ++i) {
      if (m_symmetric) {
        const int N = m_Nc - i;
        if (N == 0) {
          break;
        }
        V1.Bind(m_rows[i] + i, N);
        V2.Bind((float *) (M[i] + i), N);
      } else {
        V1.Bind(m_rows[i], m_Nc);
        V2.Bind((float *) M[i], m_Nc);
      }
      equal = V1.AssertEqual(V2, verbose1, UT::String("%s %d", str.c_str(), i),
                             epsAbs, epsRel);
    }
    if (equal) {
      return true;
    } else if (verbose2) {
      UT::PrintSeparator();
      Print(verbose2 > 1);
      UT::PrintSeparator();
      M.Print(verbose2 > 1);
      const AlignedMatrixXf E = *this - M;
      UT::PrintSeparator();
      E.Print(verbose2 > 1);
    }
    return false;
  }
  inline bool AssertEqual(const std::string fileName,
                          const int verbose = 1, const std::string str = "",
                          const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    AlignedMatrixXf M;
    M.LoadB(fileName);
    return AssertEqual(M, verbose, str, epsAbs, epsRel);
  }
  inline bool AssertZero(const int verbose = 1, const std::string str = "",
                         const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    bool zero = true;
    //const int verbose1 = verbose & 15, verbose2 = verbose >> 4;
    const int verbose1 = verbose, verbose2 = verbose;
    for (int i = 0; i < m_Nr && zero; ++i) {
      const AlignedVectorXf V(m_rows[i], m_Nc, false);
      zero = V.AssertZero(verbose1, UT::String("%s %d", str.c_str(), i),
                          epsAbs, epsRel);
    }
    if (zero) {
      return true;
    } else if (verbose2) {
      UT::PrintSeparator();
      Print(verbose2 > 1);
    }
    return false;
  }

  static inline void AddAbTo(const AlignedMatrixXf &A, const AlignedVectorXf &b, float *Ab) {
    const int Nr = A.GetRows();
    if (A.m_symmetric) {
      const int Nc = A.GetColumns();
      for (int i = 0; i < Nr; ++i) {
        for (int j = 0; j < Nc; ++j) {
          Ab[i] += A[i][j] * b[j];
        }
      }
    } else {
      for (int i = 0; i < Nr; ++i) {
        Ab[i] += SIMD::Dot<0>(Nr, A[i], b.Data());
      }
    }
  }
};

template<int M, int N>
inline void Convert(const AlignedMatrixXf &M1, AlignedMatrixMxNf<M, N> &M2) {
  const int _N = sizeof(float) * N;
  for (int i = 0; i < M; ++i) {
    memcpy(M2[i], M1[i], _N);
  }
}
template<int N>
inline void Convert(const AlignedVectorXf &V1, float *V2) {
  memcpy(V2, V1.Data(), sizeof(float) * N);
}

class AlignedMatrixXd : public AlignedMatrixX<double> {

 public:

  inline AlignedMatrixXd() {}
  inline AlignedMatrixXd(const AlignedMatrixXd &M) { *this = M; }
  inline void operator = (const AlignedMatrixXd &M) { Set(M); }

  inline void operator *= (const double s) { Scale(s); }
  inline AlignedMatrixXd operator - (const AlignedMatrixXd &B) const {
    const AlignedMatrixXd &A = *this;
    AlignedMatrixXd AmB;
    AmB.Resize(m_Nr, m_Nc, m_symmetric);
    if (m_symmetric) {
      SIMD::Subtract((m_Nc * (m_Nc + 1)) >> 1, A.Data(), B.Data(), AmB.Data());
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        SIMD::Subtract(m_Nc, A[i], B[i], AmB[i]);
      }
    }
    return AmB;
  }

  inline bool Valid() const { return m_rows[0][0] != DBL_MAX; }
  inline bool Invalid() const { return m_rows[0][0] == DBL_MAX; }
  inline void Invalidate() { m_rows[0][0] = DBL_MAX; }

  inline void Set(const AlignedMatrixXf &A) {
    Resize(A.GetRows(), A.GetColumns(), A.Symmetric());
    for (int i = 0; i < m_Nr; ++i) {
      for (int j = m_symmetric ? i : 0; j < m_Nc; ++j) {
        m_rows[i][j] = static_cast<double>(A[i][j]);
      }
    }
  }
  inline void Set(const AlignedMatrixXd &V) { AlignedMatrixX<double>::Set(V); }
  inline void Get(AlignedMatrixXf &A) const {
    A.Resize(m_Nr, m_Nc, m_symmetric);
    for (int i = 0; i < m_Nr; ++i) {
      for (int j = m_symmetric ? i : 0; j < m_Nc; ++j) {
        A[i][j] = static_cast<float>(m_rows[i][j]);
      }
    }
  }

  inline void SetBlockDiagonal(const int i, const int N, const float m) {
    const double _m = static_cast<double>(m);
    const int i2 = i + N;
    for (int _i = i; _i < i2; ++_i) {
      m_rows[_i][_i] = _m;
    }
  }

  inline void SetBlockDiagonal(const int i, const SymmetricMatrix2x2f &M) {
    int k = i;
    m_rows[k][k] = static_cast<double>(M.m00());
    m_rows[k][k + 1] = static_cast<double>(M.m01());
    ++k;
    m_rows[k][k] = static_cast<double>(M.m11());
  }
  inline void GetBlockDiagonal(const int i, SymmetricMatrix2x2f &M) const {
    int k = i;
    M.m00() = static_cast<float>(m_rows[k][k]);
    M.m01() = static_cast<float>(m_rows[k][k + 1]);
    ++k;
    M.m11() = static_cast<float>(m_rows[k][k]);
  }
  inline void IncreaseBlockDiagonal(const int i, const SymmetricMatrix2x2f &M) {
    int k = i;
    m_rows[k][k] += static_cast<double>(M.m00());
    m_rows[k][k + 1] += static_cast<double>(M.m01());
    ++k;
    m_rows[k][k] += static_cast<double>(M.m11());
  }

  inline void SetBlock(const int i, const int j, const Matrix2x3f &M) {
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, Matrix2x3f &M) const {
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 3; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix2x6f &M) {
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 6; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix2x6f &M) const {
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 6; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix2x9f &M) {
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 9; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix2x9f &M) const {
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 9; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
  }

  inline void SetBlock(const int i, const int j, const Matrix3x2f &M) {
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = 0; _j < 2; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, Matrix3x2f &M) const {
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = 0; _j < 2; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix3x3f &M) {
    const float *_M[3] = {&M.m00(), &M.m10(), &M.m20()};
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(_M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix3x3f &M) const {
    float *_M[3] = {&M.m00(), &M.m10(), &M.m20()};
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 3; ++_j) {
        _M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
    if (symmetric) {
      M.SetLowerFromUpper();
    }
  }

  inline void SetBlock(const int i, const int j, const Matrix3x3f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, Matrix3x3f &M) const {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 3; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
    if (symmetric) {
      M.SetLowerFromUpper();
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix6x6f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 6; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 6; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix6x6f &M) const {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 6; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 6; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
    if (symmetric) {
      M.SetLowerFromUpper();
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix6x9f &M) {
    for (int _i = 0; _i < 6; ++_i) {
      for (int _j = 0; _j < 9; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix6x9f &M) const {
    for (int _i = 0; _i < 6; ++_i) {
      for (int _j = 0; _j < 9; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix9x6f &M) {
    for (int _i = 0; _i < 9; ++_i) {
      for (int _j = 0; _j < 6; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix9x6f &M) const {
    for (int _i = 0; _i < 9; ++_i) {
      for (int _j = 0; _j < 6; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
  }

  inline void SetBlock(const int i, const int j, const AlignedMatrix9x9f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 9; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 9; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrix9x9f &M) const {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 9; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 9; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
    if (symmetric) {
      M.SetLowerFromUpper();
    }
  }
  
  inline void SetBlock(const int i, const int j, const AlignedMatrixXf &M) {
    const int Nr = M.GetRows(), Nc = M.GetColumns();
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < Nr; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < Nc; ++_j) {
        m_rows[i + _i][j + _j] = static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void GetBlock(const int i, const int j, AlignedMatrixXf &M) const {
    const int Nr = M.GetRows(), Nc = M.GetColumns();
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < Nr; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < Nc; ++_j) {
        M[_i][_j] = static_cast<float>(m_rows[i + _i][j + _j]);
      }
    }
    if (symmetric) {
      M.SetLowerFromUpper();
    }
  }
  
  inline void IncreaseBlock(const int i, const int j, const AlignedMatrix2x3f &M) {
    const float *_M[2] = {&M.m00(), &M.m10()};
    for (int _i = 0; _i < 2; ++_i) {
      for (int _j = 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] += static_cast<double>(_M[_i][_j]);
      }
    }
  }
  inline void IncreaseBlock(const int i, const int j, const Matrix3x2f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = 0; _j < 2; ++_j) {
        m_rows[i + _i][j + _j] += static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void IncreaseBlock(const int i, const int j, const AlignedMatrix3x3f &M) {
    const float *_M[3] = {&M.m00(), &M.m10(), &M.m20()};
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 3; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 3; ++_j) {
        m_rows[i + _i][j + _j] += static_cast<double>(_M[_i][_j]);
      }
    }
  }

  inline void IncreaseBlock(const int i, const int j, const AlignedMatrix6x6f &M) {
    const bool symmetric = m_symmetric && i == j;
    for (int _i = 0; _i < 6; ++_i) {
      for (int _j = symmetric ? _i : 0; _j < 6; ++_j) {
        m_rows[i + _i][j + _j] += static_cast<double>(M[_i][_j]);
      }
    }
  }
  inline void IncreaseBlockDiagonal(const int i, const SymmetricMatrix6x6f &M) {
    const float *_M = M;
    for (int _i = 0, k = 0; _i < 6; ++_i) {
      for (int _j = _i; _j < 6; ++_j, ++k) {
        m_rows[i + _i][i + _j] += static_cast<double>(_M[k]);
      }
    }
  }

  inline void SetLowerFromUpper() {
    const int N = std::min(m_Nr, m_Nc);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < i; ++j) {
        m_rows[i][j] = m_rows[j][i];
      }
    }
  }

  inline void Scale(const double s) {
    if (m_symmetric) {
      SIMD::Multiply((m_Nc * (m_Nc + 1)) >> 1, s, Data());
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        SIMD::Multiply(m_Nc, s, m_rows[i]);
      }
    }
  }

  inline bool SolveLDL(AlignedVectorXd &b, const double *eps = NULL,
                       const bool decomposed = false) {
    return LS::SolveLDL(m_Nr, m_rows.data(), b.Data(), eps, decomposed);
  }
  inline int RankLDL(AlignedVectorXd *work, const double *eps = NULL) {
    work->Resize(m_Nc);
    return LS::RankLDL(m_Nr, m_rows.data(), work->Data(), eps);
  }
  inline bool InverseLDL(const double *eps = NULL, const bool decomposed = false) {
    if (LS::InverseLDL(m_Nr, m_rows.data(), eps, decomposed)) {
      SetLowerFromUpper();
      return true;
    } else {
      Invalidate();
      return false;
    }
  }
  inline bool MarginalizeLDL(const int i, const int N, AlignedVectorXd &b,
                             AlignedVectorXd *work, const double *eps = NULL) {
    work->Resize(m_Nc);
    const bool scc = LS::MarginalizeLDL(m_Nr, m_rows.data(), b.Data(), i, N, work->Data(), eps);
    Erase(i, N);
    b.Erase(i, N);
    return scc;
  }
  inline void Print(const bool e = false) const {
    if (m_symmetric) {
      for (int i = 0; i < m_Nr; ++i) {
        for (int j = 0; j < i; ++j) {
          if (e) {
            UT::Print("%e ", m_rows[j][i]);
          } else {
            UT::Print("%f ", m_rows[j][i]);
          }
        }
        for (int j = i; j < m_Nc; ++j) {
          if (e) {
            UT::Print("%e ", m_rows[i][j]);
          } else {
            UT::Print("%f ", m_rows[i][j]);
          }
        }
        UT::Print("\n");
      }
    } else {
      for (int i = 0; i < m_Nr; ++i) {
        for (int j = 0; j < m_Nc; ++j) {
          if (e) {
            UT::Print("%e ", m_rows[i][j]);
          } else {
            UT::Print("%f ", m_rows[i][j]);
          }
        }
        UT::Print("\n");
      }
    }
  }
  inline void PrintRow(const int i, const bool e = false) const {
    for (int j = 0; j < m_Nc; ++j) {
      if (e) {
        UT::Print("%e ", m_rows[i][j]);
      } else {
        UT::Print("%f ", m_rows[i][j]);
      }
    }
    UT::Print("\n");
  }
  inline bool AssertEqual(const AlignedMatrixXd &M,
                          const int verbose = 1, const std::string str = "",
                          const double epsAbs = 0.0, const double epsRel = 0.0) const {
    AlignedVectorXd V1, V2;
    bool equal = M.m_Nr == m_Nr && M.m_Nc == m_Nc && M.m_symmetric == m_symmetric;
    //const int verbose1 = verbose & 15, verbose2 = verbose >> 4;
    const int verbose1 = verbose, verbose2 = verbose;
    for (int i = 0; i < m_Nr && equal; ++i) {
      if (m_symmetric) {
        const int N = m_Nc - i;
        if (N == 0) {
          break;
        }
        V1.Bind(m_rows[i] + i, N);
        V2.Bind((float *) (M[i] + i), N);
      } else {
        V1.Bind(m_rows[i], m_Nc);
        V2.Bind((float *) M[i], m_Nc);
      }
      equal = V1.AssertEqual(V2, verbose1, UT::String("%s %d", str.c_str(), i),
                             epsAbs, epsRel);
    }
    if (equal) {
      return true;
    } else if (verbose2) {
      UT::PrintSeparator();
      Print(verbose2 > 1);
      UT::PrintSeparator();
      M.Print(verbose2 > 1);
      const AlignedMatrixXd E = *this - M;
      UT::PrintSeparator();
      E.Print(verbose2 > 1);
    }
    return false;
  }
  inline bool AssertZero(const int verbose = 1, const std::string str = "",
                         const double epsAbs = 0.0, const double epsRel = 0.0) const {
    bool zero = true;
    //const int verbose1 = verbose & 15, verbose2 = verbose >> 4;
    const int verbose1 = verbose, verbose2 = verbose;
    for (int i = 0; i < m_Nr && zero; ++i) {
      const AlignedVectorXd V(m_rows[i], m_Nc, false);
      zero = V.AssertZero(verbose1, UT::String("%s %d", str.c_str(), i),
                          epsAbs, epsRel);
    }
    if (zero) {
      return true;
    } else if (verbose2) {
      UT::PrintSeparator();
      Print(verbose2 > 1);
    }
    return false;
  }

  static inline void AddAbTo(const AlignedMatrixXd &A, const AlignedVectorXd &b, double *Ab) {
    const int Nr = A.GetRows();
    if (A.m_symmetric) {
      const int Nc = A.GetColumns();
      for (int i = 0; i < Nr; ++i) {
        for (int j = 0; j < Nc; ++j) {
          Ab[i] += A[i][j] * b[j];
        }
      }
    } else {
      for (int i = 0; i < Nr; ++i) {
        Ab[i] += SIMD::Dot<0>(Nr, A[i], b.Data());
      }
    }
  }
};
}

#endif
