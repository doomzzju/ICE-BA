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
#ifndef _MATRIX_4x4_H_
#define _MATRIX_4x4_H_

#include "Matrix3x3.h"
#include "Vector4.h"
#include "LinearSystem.h"

namespace LA {

class AlignedMatrix4x4f {

 public:

  inline const xp128f& m_00_01_02_03() const { return m_data[0]; }  inline xp128f& m_00_01_02_03() { return m_data[0]; }
  inline const xp128f& m_10_11_12_13() const { return m_data[1]; }  inline xp128f& m_10_11_12_13() { return m_data[1]; }
  inline const xp128f& m_20_21_22_23() const { return m_data[2]; }  inline xp128f& m_20_21_22_23() { return m_data[2]; }
  inline const xp128f& m_30_31_32_33() const { return m_data[3]; }  inline xp128f& m_30_31_32_33() { return m_data[3]; }

  inline const float& m00() const { return m_data[0][0]; }        inline float& m00() { return m_data[0][0]; }
  inline const float& m01() const { return m_data[0][1]; }        inline float& m01() { return m_data[0][1]; }
  inline const float& m02() const { return m_data[0][2]; }        inline float& m02() { return m_data[0][2]; }
  inline const float& m03() const { return m_data[0][3]; }        inline float& m03() { return m_data[0][3]; }
  inline const float& m10() const { return m_data[1][0]; }        inline float& m10() { return m_data[1][0]; }
  inline const float& m11() const { return m_data[1][1]; }        inline float& m11() { return m_data[1][1]; }
  inline const float& m12() const { return m_data[1][2]; }        inline float& m12() { return m_data[1][2]; }
  inline const float& m13() const { return m_data[1][3]; }        inline float& m13() { return m_data[1][3]; }
  inline const float& m20() const { return m_data[2][0]; }        inline float& m20() { return m_data[2][0]; }
  inline const float& m21() const { return m_data[2][1]; }        inline float& m21() { return m_data[2][1]; }
  inline const float& m22() const { return m_data[2][2]; }        inline float& m22() { return m_data[2][2]; }
  inline const float& m23() const { return m_data[2][3]; }        inline float& m23() { return m_data[2][3]; }
  inline const float& m30() const { return m_data[3][0]; }        inline float& m30() { return m_data[3][0]; }
  inline const float& m31() const { return m_data[3][1]; }        inline float& m31() { return m_data[3][1]; }
  inline const float& m32() const { return m_data[3][2]; }        inline float& m32() { return m_data[3][2]; }
  inline const float& m33() const { return m_data[3][3]; }        inline float& m33() { return m_data[3][3]; }

  inline operator const float* () const { return (const float *) this; }
  inline operator       float* ()       { return (      float *) this; }
  inline const float& operator() (int row, int col) const {
    return m_data[row][col];
  }

  inline float& operator() (int row, int col) {
    return m_data[row][col];
  }
  inline void operator = (const AlignedMatrix3x3f &M) {
    memcpy(&m00(), &M.m00(), 12); m03() = 0.0f;
    memcpy(&m10(), &M.m10(), 12); m13() = 0.0f;
    memcpy(&m20(), &M.m20(), 12); m23() = 0.0f;
    m_30_31_32_33().vset_all_lane(0.0f, 0.0f, 0.0f, 1.0f);
  }
  inline AlignedMatrix4x4f operator - (const AlignedMatrix4x4f &B) const {
    AlignedMatrix4x4f AmB;
    AmB.m_00_01_02_03() = m_00_01_02_03() - B.m_00_01_02_03();
    AmB.m_10_11_12_13() = m_10_11_12_13() - B.m_10_11_12_13();
    AmB.m_20_21_22_23() = m_20_21_22_23() - B.m_20_21_22_23();
    AmB.m_30_31_32_33() = m_30_31_32_33() - B.m_30_31_32_33();
    return AmB;
  }

  inline void Get(AlignedMatrix3x3f &M) const {
    memcpy(&M, this, sizeof(AlignedMatrix3x3f));
  }

  inline void MakeZero() {
    memset(this, 0, sizeof(AlignedMatrix4x4f));
  }
  inline void MakeIdentity() {
    MakeZero();
    m00() = m11() = m22() = m33() = 1.0f;
  }
  inline void IncreaseDiagonal(const float d) {
    m00() = d + m00();
    m11() = d + m11();
    m22() = d + m22();
    m33() = d + m33();
  }

  inline void SetLowerFromUpper() {
    m10() = m01();
    m20() = m02();
    m30() = m03();
    m21() = m12();
    m31() = m13();
    m32() = m23();
  }

  inline void Random(const float mMax) { UT::Random(&m00(), 16, mMax); }

  inline void Print(const bool e = false) const {
    if (e) {
      UT::Print("%e %e %e %e\n", m00(), m01(), m02(), m03());
      UT::Print("%e %e %e %e\n", m10(), m11(), m12(), m13());
      UT::Print("%e %e %e %e\n", m20(), m21(), m22(), m23());
      UT::Print("%e %e %e %e\n", m30(), m31(), m32(), m33());
    } else {
      UT::Print("%f %f %f %f\n", m00(), m01(), m02(), m03());
      UT::Print("%f %f %f %f\n", m10(), m11(), m12(), m13());
      UT::Print("%f %f %f %f\n", m20(), m21(), m22(), m23());
      UT::Print("%f %f %f %f\n", m30(), m31(), m32(), m33());
    }
  }
  inline bool AssertEqual(const AlignedMatrix4x4f &M,
                          const int verbose = 1, const std::string str = "",
                          const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    if (UT::VectorAssertEqual(&m00(), &M.m00(), 16, verbose, str, epsAbs, epsRel)) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
      UT::PrintSeparator();
      M.Print(verbose > 1);
      UT::PrintSeparator();
      const AlignedMatrix4x4f E = *this - M;
      E.Print(verbose > 1);
    }
    return false;
  }

 protected:
   xp128f m_data[4];
};

template<typename TYPE> class SymmetricMatrix4x4 {
 public:
  inline const TYPE& m00() const { return m_data[0]; }  inline TYPE& m00() { return m_data[0]; }
  inline const TYPE& m01() const { return m_data[1]; }  inline TYPE& m01() { return m_data[1]; }
  inline const TYPE& m02() const { return m_data[2]; }  inline TYPE& m02() { return m_data[2]; }
  inline const TYPE& m03() const { return m_data[3]; }  inline TYPE& m03() { return m_data[3]; }
  inline const TYPE& m11() const { return m_data[4]; }  inline TYPE& m11() { return m_data[4]; }
  inline const TYPE& m12() const { return m_data[5]; }  inline TYPE& m12() { return m_data[5]; }
  inline const TYPE& m13() const { return m_data[6]; }  inline TYPE& m13() { return m_data[6]; }
  inline const TYPE& m22() const { return m_data[7]; }  inline TYPE& m22() { return m_data[7]; }
  inline const TYPE& m23() const { return m_data[8]; }  inline TYPE& m23() { return m_data[8]; }
  inline const TYPE& m33() const { return m_data[9]; }  inline TYPE& m33() { return m_data[9]; }

  inline operator const TYPE* () const { return m_data; }
  inline operator     TYPE* ()     { return m_data; }

  inline void MakeZero() { memset(this, 0, sizeof(SymmetricMatrix4x4<TYPE>)); }

  inline void Get(AlignedMatrix4x4f &M) const;

  static inline void AAT(const AlignedMatrix4x4f &A, SymmetricMatrix4x4<TYPE> &AAT);

 protected:

  TYPE m_data[10];
};

typedef SymmetricMatrix4x4<float> SymmetricMatrix4x4f;
typedef SymmetricMatrix4x4<double> SymmetricMatrix4x4d;

template<> inline void SymmetricMatrix4x4f::Get(AlignedMatrix4x4f &M) const {
  M.m_00_01_02_03().vset_all_lane(m00(), m01(), m02(), m03());
  M.m_10_11_12_13().vset_all_lane(m01(), m11(), m12(), m13());
  M.m_20_21_22_23().vset_all_lane(m02(), m12(), m22(), m23());
  M.m_30_31_32_33().vset_all_lane(m03(), m13(), m23(), m33());
}

}

#endif
