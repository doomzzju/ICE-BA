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
#ifndef _VECTOR_2_H_
#define _VECTOR_2_H_

#include "Utility.h"

namespace LA {
template<typename TYPE> class Vector2 {
 public:
  inline Vector2() {}
  inline Vector2(const TYPE *v) { Set(v[0], v[1]); }
  inline Vector2(const TYPE v0, const TYPE v1) { Set(v0, v1); }

  inline const TYPE& v0() const { return m_data[0]; }   inline TYPE& v0() { return m_data[0]; }
  inline const TYPE& v1() const { return m_data[1]; }   inline TYPE& v1() { return m_data[1]; }
  inline const TYPE& x () const { return m_data[0]; }   inline TYPE& x () { return m_data[0]; }
  inline const TYPE& y () const { return m_data[1]; }   inline TYPE& y () { return m_data[1]; }

  inline operator const TYPE* () const { return m_data; }
  inline operator       TYPE* ()       { return m_data; }
  inline const float& operator() (const int row, const int col = 0) const {
    return m_data[row];
  }
  inline float& operator() (const int row, const int col = 0) {
    return m_data[row];
  }
  inline bool operator == (const Vector2<TYPE> &v) const {
    return v0() == v.v0() && v1() == v.v1();
  }
  inline void operator += (const Vector2<TYPE> &v) {
    v0() = v.v0() + v0(); v1() = v.v1() + v1();
  }
  inline void operator -= (const Vector2<TYPE> &v) {
    v0() = -v.v0() + v0(); v1() = -v.v1() + v1();
  }
  inline void operator *= (const TYPE s) { v0() *= s; v1() *= s; }
  inline void operator *= (const Vector2<TYPE> &s) { v0() *= s.v0(); v1() *= s.v1(); }
  inline Vector2<TYPE> operator + (const Vector2<TYPE> &v) const {
    return Vector2<TYPE>(v0() + v.v0(), v1() + v.v1());
  }
  inline Vector2<TYPE> operator - (const Vector2<TYPE> &v) const {
    return Vector2<TYPE>(v0() - v.v0(), v1() - v.v1());
  }
  inline Vector2<TYPE> operator * (const TYPE s) const {
    Vector2<TYPE> v;
    GetScaled(s, v);
    return v;
  }
  inline Vector2<TYPE> operator * (const Vector2<TYPE> &s) const {
    Vector2<TYPE> v;
    GetScaled(s, v);
    return v;
  }

  inline void Set(const TYPE v0, const TYPE v1) { this->v0() = v0; this->v1() = v1; }
  inline void Set(const float *v);
  inline void Set(const double *v);

  inline void MakeZero() { memset(this, 0, sizeof(Vector2<TYPE>)); }
  inline void MakeMinus() { v0() = -v0(); v1() = -v1(); }
  inline void MakeSquareRoot() { v0() = sqrt(v0()); v1() = sqrt(v1()); }

  inline bool Valid() const { return v0() != UT::Invalid<TYPE>(); }
  inline bool Invalid() const { return v0() == UT::Invalid<TYPE>(); }
  inline void Invalidate() { v0() = UT::Invalid<TYPE>(); }

  inline void Get(float *v) const;
  inline void Get(double *v) const;
  inline void GetScaled(const TYPE s, Vector2<TYPE> &v) const { v.v0() = s * v0(); v.v1() = s * v1(); }
  inline void GetScaled(const Vector2<TYPE> &s, Vector2<TYPE> &v) const { v.v0() = s.v0() * v0(); v.v1() = s.v1() * v1(); }
  inline TYPE SquaredLength() const { return v0() * v0() + v1() * v1(); }
  inline TYPE Dot(const Vector2<TYPE> &v) const { return v0() * v.v0() + v1() * v.v1(); }
  inline TYPE Dot(const TYPE v0, const TYPE v1) const { return this->v0() * v0 + this->v1() * v1; }
  inline TYPE Sum() const { return v0() + v1(); }
  inline void SaveB(FILE *fp) const { UT::SaveB(*this, fp); }
  inline void LoadB(FILE *fp) { UT::LoadB(*this, fp); }
  inline void Print(const bool e = false, const bool n = true) const {
    UT::PrintValue<TYPE>(v0(), e);
    UT::Print(" ");
    UT::PrintValue<TYPE>(v1(), e);
    if (n)
      UT::Print("\n");
  }
  inline void Print(const std::string str, const bool e, const bool n) const {
    UT::Print("%s", str.c_str());
    Print(e, n);
  }
  inline bool AssertEqual(const Vector2<TYPE> &v,
                          const int verbose = 1, const std::string str = "",
                          const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    if (UT::VectorAssertEqual<TYPE>(&v0(), &v.v0(), 2, verbose, str, epsAbs, epsRel)) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
      v.Print(verbose > 1);
      const Vector2<TYPE> e = *this - v;
      e.Print(verbose > 1);
    }
    return false;
  }
  inline bool AssertZero(const int verbose = 1, const std::string str = "",
                         const float epsAbs = 0.0f, const float epsRel = 0.0f) const {
    if (UT::VectorAssertZero(&v0(), 2, verbose, str, epsAbs, epsRel)) {
      return true;
    } else if (verbose) {
      UT::PrintSeparator();
      Print(verbose > 1);
    }
    return false;
  }

  inline void Random(const TYPE vMax) { UT::Random<TYPE>(m_data, 2, -vMax, vMax); }
  inline void Random(const TYPE vMin, const TYPE vMax) { UT::Random<TYPE>(m_data, 2, vMin, vMax); }
  static inline Vector2<TYPE> GetRandom(const TYPE vMax) { Vector2<TYPE> v; v.Random(vMax); return v; }
  static inline Vector2<TYPE> GetRandom(const TYPE vMin, const TYPE vMax) { Vector2<TYPE> v; v.Random(vMin, vMax); return v; }

  static inline void apb(const Vector2<TYPE> &a, const Vector2<TYPE> &b, Vector2<TYPE> &apb) {
    apb.v0() = a.v0() + b.v0();
    apb.v1() = a.v1() + b.v1();
  }
  static inline void amb(const Vector2<TYPE> &a, const Vector2<TYPE> &b, Vector2<TYPE> &amb) {
    amb.v0() = a.v0() - b.v0();
    amb.v1() = a.v1() - b.v1();
  }

 protected:

  TYPE m_data[2];
};

typedef Vector2<short>  Vector2s;
typedef Vector2<ushort> Vector2us;
typedef Vector2<int   > Vector2i;
typedef Vector2<float > Vector2f;
typedef Vector2<double> Vector2d;

template<> inline void Vector2f::Set(const float *v) { memcpy(this, v, sizeof(Vector2f)); }
template<> inline void Vector2d::Set(const double *v) { memcpy(this, v, sizeof(Vector2d)); }
template<> inline void Vector2f::Set(const double *v) {
  v0() = static_cast<float>(v[0]);
  v1() = static_cast<float>(v[1]);
}
template<> inline void Vector2d::Set(const float *v) {
  v0() = static_cast<double>(v[0]);
  v1() = static_cast<double>(v[1]);
}
template<> inline void Vector2f::Get(float *v) const { memcpy(v, this, sizeof(Vector2f)); }
template<> inline void Vector2d::Get(double *v) const { memcpy(v, this, sizeof(Vector2d)); }
template<> inline void Vector2f::Get(double *v) const {
  v[0] = static_cast<double>(v0());
  v[1] = static_cast<double>(v1());
}
template<> inline void Vector2d::Get(float *v) const {
  v[0] = static_cast<float>(v0());
  v[1] = static_cast<float>(v1());
}

}

#endif
