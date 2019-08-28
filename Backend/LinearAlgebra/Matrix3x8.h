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
#ifndef _MATRIX_3x8_H_
#define _MATRIX_3x8_H_

#include "Matrix3x7.h"

namespace LA {

class AlignedMatrix3x8f {

 public:
  inline const xp128f& m_00_01_02_03() const { return m_data[0]; }  inline xp128f& m_00_01_02_03() { return m_data[0]; }
  inline const xp128f& m_04_05_06_07() const { return m_data[1]; }  inline xp128f& m_04_05_06_07() { return m_data[1]; }
  inline const xp128f& m_10_11_12_13() const { return m_data[2]; }  inline xp128f& m_10_11_12_13() { return m_data[2]; }
  inline const xp128f& m_14_15_16_17() const { return m_data[3]; }  inline xp128f& m_14_15_16_17() { return m_data[3]; }
  inline const xp128f& m_20_21_22_23() const { return m_data[4]; }  inline xp128f& m_20_21_22_23() { return m_data[4]; }
  inline const xp128f& m_24_25_26_27() const { return m_data[5]; }  inline xp128f& m_24_25_26_27() { return m_data[5]; }
  inline const float& m00() const { return m_data[0][0]; } inline float& m00() { return m_data[0][0]; }
  inline const float& m01() const { return m_data[0][1]; } inline float& m01() { return m_data[0][1]; }
  inline const float& m02() const { return m_data[0][2]; } inline float& m02() { return m_data[0][2]; }
  inline const float& m03() const { return m_data[0][3]; } inline float& m03() { return m_data[0][3]; }
  inline const float& m04() const { return m_data[1][0]; } inline float& m04() { return m_data[1][0]; }
  inline const float& m05() const { return m_data[1][1]; } inline float& m05() { return m_data[1][1]; }
  inline const float& m06() const { return m_data[1][2]; } inline float& m06() { return m_data[1][2]; }
  inline const float& m07() const { return m_data[1][3]; } inline float& m07() { return m_data[1][3]; }
  inline const float& m10() const { return m_data[2][0]; } inline float& m10() { return m_data[2][0]; }
  inline const float& m11() const { return m_data[2][1]; } inline float& m11() { return m_data[2][1]; }
  inline const float& m12() const { return m_data[2][2]; } inline float& m12() { return m_data[2][2]; }
  inline const float& m13() const { return m_data[2][3]; } inline float& m13() { return m_data[2][3]; }
  inline const float& m14() const { return m_data[3][0]; } inline float& m14() { return m_data[3][0]; }
  inline const float& m15() const { return m_data[3][1]; } inline float& m15() { return m_data[3][1]; }
  inline const float& m16() const { return m_data[3][2]; } inline float& m16() { return m_data[3][2]; }
  inline const float& m17() const { return m_data[3][3]; } inline float& m17() { return m_data[3][3]; }
  inline const float& m20() const { return m_data[4][0]; } inline float& m20() { return m_data[4][0]; }
  inline const float& m21() const { return m_data[4][1]; } inline float& m21() { return m_data[4][1]; }
  inline const float& m22() const { return m_data[4][2]; } inline float& m22() { return m_data[4][2]; }
  inline const float& m23() const { return m_data[4][3]; } inline float& m23() { return m_data[4][3]; }
  inline const float& m24() const { return m_data[5][0]; } inline float& m24() { return m_data[5][0]; }
  inline const float& m25() const { return m_data[5][1]; } inline float& m25() { return m_data[5][1]; }
  inline const float& m26() const { return m_data[5][2]; } inline float& m26() { return m_data[5][2]; }
  inline const float& m27() const { return m_data[5][3]; } inline float& m27() { return m_data[5][3]; }

  inline operator const float* () const { return (const float *) this; }
  inline operator       float* ()       { return (      float *) this; }

  inline const float& operator() (int row, int col) const {
    return data[row * 8 + col];
  }

  inline float& operator() (int row, int col) {
    return data[row * 8 + col];
  }

  inline void MakeZero() { memset(this, 0, sizeof(AlignedMatrix3x8f)); }

  inline void Print(const bool e = false) const {
    if (e) {
      UT::Print("%e %e %e %e %e %e %e %e\n", m00(), m01(), m02(), m03(), m04(), m05(), m06(), m07());
      UT::Print("%e %e %e %e %e %e %e %e\n", m10(), m11(), m12(), m13(), m14(), m15(), m16(), m17());
      UT::Print("%e %e %e %e %e %e %e %e\n", m20(), m21(), m22(), m23(), m24(), m25(), m26(), m27());
    } else {
      UT::Print("%f %f %f %f %f %f %f %f\n", m00(), m01(), m02(), m03(), m04(), m05(), m06(), m07());
      UT::Print("%f %f %f %f %f %f %f %f\n", m10(), m11(), m12(), m13(), m14(), m15(), m16(), m17());
      UT::Print("%f %f %f %f %f %f %f %f\n", m20(), m21(), m22(), m23(), m24(), m25(), m26(), m27());
    }
  }
  inline bool AssertEqual(const AlignedMatrix3x8f &M, const int verbose = 1) const {
    if (UT::VectorAssertEqual(&m00(), &M.m00(), 24, 0)) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
      UT::PrintSeparator();
      M.Print(verbose > 1);
    }
    return false;
  }

 protected:
  union {
    xp128f m_data[6];
    float data[24];
  };
};
}

#endif
