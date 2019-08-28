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
//#define CAMERA_PRIOR_DEBUG_EIGEN
#include "CameraPrior.h"

namespace CameraPrior {

bool Pose::Marginalize(const int i, AlignedVector<float> *work, const float *eps) {
  const int N = static_cast<int>(m_iKFs.size());
  Matrix::X A;
  Vector::X b, _eps, _work;
  const int Npg = m_Zps.Size() == N ? 0 : 2, Np = Npg + N * 6;
  work->Resize((A.BindSize(Np, Np, true) + b.BindSize(Np) + _work.BindSize(Np) +
               (eps ? _eps.BindSize(6) : 0)) / sizeof(float));
  A.Bind(work->Data(), Np, Np, true);
  b.Bind(A.BindNext(), Np);
  GetPriorEquation(&A, &b);
  _work.Bind(b.BindNext(), Np);
  if (eps) {
    _eps.Bind(_work.BindNext(), 6);
    _eps.Set(eps);
  }
  const bool scc = A.MarginalizeLDL(Npg + i * 6, 6, b, &_work, eps ? _eps.Data() : NULL);
  Erase(i);
  SetPriorEquation(A, b);
  return scc;
}

bool Pose::MarginalizeUninformative(const float w, const float s2p, const float s2r,
                                    std::vector<int> *iks, AlignedVector<float> *work,
                                    const float *eps) {
  LA::AlignedMatrixXf S;
  if (!GetPriorMeasurement(1.0f, &S, NULL, NULL, work, eps)) {
    return false;
  }
  LA::Vector3f s2pi, s2ri;
  const float _s2p = s2p / w, _s2r = s2r / w;
  const int N = static_cast<int>(m_iKFs.size());
  iks->resize(N);
  for (int i = 0, ip = m_Zps.Size() == N ? 0 : 2; i < N; ++i) {
    S.GetDiagonal(ip, s2pi);  ip += 3;
    S.GetDiagonal(ip, s2ri);  ip += 3;
    if (s2pi.Maximal() > _s2p || s2ri.Maximal() > _s2r) {
      iks->at(i) = i;
    } else {
      iks->at(i) = -1;
    }
  }
  for (int i = 0; i < N; ++i) {
    const int _i = iks->at(i);
    if (_i == -1) {
      continue;
    }
    Marginalize(_i, work, eps);
    iks->at(i) = -1;
    for (int j = i + 1; j < N; ++j) {
      if (iks->at(j) != -1) {
        --iks->at(j);
      }
    }
  }
  //if (m_iKFs.empty()) {
  //  return false;
  //}
  LA::AlignedVectorXf x;
  return GetPriorMeasurement(1.0f, &S, &x, &m_xTb, work, eps);
}

void Pose::SetPriorEquation(const Matrix::X &A, const Vector::X &b, const bool g) {
  const int N = static_cast<int>(m_iKFs.size());
  const bool _g = m_Zps.Size() != N;
  const int ipg = 0, ipc = _g ? 2 : 0;
  if (g && _g) {
    A.GetBlockDiagonal(ipg, m_Arr);
    m_Arc.Resize(N);
    for (int i = 0, ip = ipc; i < N; ++i, ip += 6) {
      A.GetBlock(ipg, ip, m_Arc[i]);
    }
    b.GetBlock(ipg, m_br);
  } else {
    m_Arc.Resize(0);
  }
  m_Acc.Resize(N, N, true);
  m_bc.Resize(N);
  for (int i = 0, ip = ipc; i < N; ++i, ip += 6) {
    for (int j = i, jp = ip; j < N; ++j, jp += 6) {
      A.GetBlock(ip, jp, m_Acc[i][j]);
    }
    b.GetBlock(ip, m_bc[i]);
  }
}

void Pose::GetPriorEquation(Matrix::X *A, Vector::X *b, const bool symmetric, const bool g) const {
  const int N = static_cast<int>(m_iKFs.size());
  const bool _g = m_Zps.Size() != N;
  const int Npg = _g ? 2 : 0;
  const int Np = Npg + N * 6;
  const int ipg = 0, ipc = Npg;
  A->Resize(Np, Np, symmetric);
  if (b) {
    b->Resize(Np);
  }

  if (g && _g) {
    A->SetBlockDiagonal(ipg, m_Arr);
    for (int i = 0, ip = ipc; i < N; ++i, ip += 6) {
      A->SetBlock(ipg, ip, m_Arc[i]);
    }
    if (b) {
      b->SetBlock(ipg, m_br);
    }
  }
  for (int i = 0, ip = ipc; i < N; ++i, ip += 6) {
    for (int j = i, jp = ip; j < N; ++j, jp += 6) {
      A->SetBlock(ip, jp, m_Acc[i][j]);
    }
    if (b) {
      b->SetBlock(ip, m_bc[i]);
    }
  }
}

bool Pose::GetPriorMeasurement(const Element::T w, Matrix::X *S, Vector::X *x,
                               Element::T *xTb, const Element::T *eps) const {
  GetPriorEquation(S, x, false);
  if (x) {
    Vector::X b;
    if (xTb) {
      const int Np = x->Size(), NpC = SIMD::Ceil<Element::T>(Np);
      b.Bind(x->Data() + NpC, Np);
      b.Copy(*x);
    }
    if (!S->SolveLDL(*x, eps)) {
      return false;
    }
    x->MakeMinus();
    S->InverseLDL(eps, true);
    if (xTb) {
      *xTb = x->Dot(b);
    }
  } else {
    if (!S->InverseLDL(eps)) {
      return false;
    }
  }
  if (w != 1) {
    *S *= w;
    if (xTb) {
      *xTb /= w;
    }
  }
  return true;
}

bool Pose::GetPriorMeasurement(const float w, LA::AlignedMatrixXf *S, LA::AlignedVectorXf *x,
                               float *xTb, AlignedVector<float> *work, const float *eps) const {
  if (Invalid()) {
    return false;
  }
  const int N = static_cast<int>(m_iKFs.size());
  const bool g = m_Zps.Size() != N;
  const int Npg = g ? 2 : 0;
  const int Np = Npg + N * 6;
  const int ipg = 0;

  Matrix::X _S;
  Vector::X _x, _eps;
  Element::T _xTb;
  const int _Nx = xTb ? (SIMD::Ceil<Element::T>(Np) + Np) : Np;
  work->Resize((S->BindSize(Np, Np, true) + (x ? x->BindSize(Np) : 0) +
                _S.BindSize(Np, Np) + (x ? _x.BindSize(_Nx) : 0) +
               (eps ? _eps.BindSize(Np) : 0)) / sizeof(float));
  S->Bind(work->Data(), Np, Np, true);
  if (x) {
    x->Bind(S->BindNext(), Np);
  }
  _S.Bind(x ? x->BindNext() : S->BindNext(), Np, Np);
  if (x) {
    _x.Bind(_S.BindNext(), _Nx);
  }
  if (eps) {
    _eps.Bind(x ? _x.BindNext() : _S.BindNext(), Np);
    if (g) {
      LA::Vector2<Element::T> *eg = (LA::Vector2<Element::T> *) (_eps.Data() + ipg);
      //eg->MakeZero();
      eg->Set(eps + 3);
    }
    LA::Vector6<Element::T> *ecs = (LA::Vector6<Element::T> *) (_eps.Data() + Npg);
    for (int i = 0; i < N; ++i) {
      ecs[i].Set(eps);
    }
  }
  if (!GetPriorMeasurement(static_cast<Element::T>(w), &_S, x ? &_x : NULL,
                           xTb ? &_xTb : NULL, eps ? _eps.Data() : NULL)) {
    S->Resize(0, 0, true);
    if (x) {
      x->Resize(0);
      if (xTb) {
        *xTb = 0.0f;
      }
    }
    return false;
  }
  if (x) {
    _x.Get(*x);
    if (xTb) {
      *xTb = static_cast<float>(_xTb);
    }
  }
  _S.Get(*S);

  return true;
}

bool Motion::GetPriorMeasurement(const float w, LA::AlignedMatrixXf *S, LA::AlignedVectorXf *x,
                                 float *xTb, AlignedVector<float> *work, const float *eps) const {
  if (Invalid()) {
    return false;
  }

  Matrix::X _S;
  Vector::X _x, _eps;
  Element::T _xTb;
  const int _Nx = xTb ? (SIMD::Ceil<Element::T>(9) + 9) : 9;
  work->Resize((S->BindSize(9, 9) + (x ? x->BindSize(9) : 0) +
                _S.BindSize(9, 9) + (x ? _x.BindSize(_Nx) : 0) +
               (eps ? _eps.BindSize(9) : 0)) / sizeof(float));
  S->Bind(work->Data(), 9, 9);
  if (x) {
    x->Bind(S->BindNext(), 9);
  }
  _S.Bind(x ? x->BindNext() : S->BindNext(), 9, 9);
  if (x) {
    _x.Bind(_S.BindNext(), _Nx);
  }
  if (eps) {
    _eps.Bind(x ? _x.BindNext() : _S.BindNext(), 9);
    _eps.Set(eps);
  }
  _S.SetBlock(0, 0, m_Amm);
  if (x) {
    _x.SetBlock(0, m_bm);
  }

  if (x) {
    Vector::X _b;
    if (xTb) {
      _x.Resize(9);
      _b.Bind(_x.Data() + 9, 9);
      _b.Copy(_x);
    }
    if (!_S.SolveLDL(_x, eps ? _eps.Data() : NULL)) {
      S->Resize(0, 0);
      if (x) {
        x->Resize(0);
        if (xTb) {
          *xTb = 0.0f;
        }
      }
      return false;
    }
    _x.MakeMinus();
    _S.InverseLDL(eps ? _eps.Data() : NULL, true);
    if (xTb) {
      _xTb = _x.Dot(_b);
    }
  } else {
    if (!_S.InverseLDL(eps ? _eps.Data() : NULL)) {
      return false;
    }
  }
  if (w != 1.0) {
    const Element::T _w = static_cast<Element::T>(w);
    _S *= _w;
    if (xTb) {
      _xTb /= _w;
    }
  }
  _S.GetBlock(0, 0, *S);
  if (x) {
    _x.GetBlock(0, *x);
    if (xTb) {
      *xTb = static_cast<float>(_xTb);
    }
  }
  return true;
}

void Joint::SetPriorEquation(const Matrix::X &A, const Vector::X &b) {
  const int N = static_cast<int>(m_iKFs.size());
  const bool g = m_Zps.Size() != N;
  const int Npg = g ? 2 : 0;

  const int Npc = N * 6;
  const int ipg = 0, ipc = Npg, ipm = ipc + Npc;
  Pose::SetPriorEquation(A, b);
  if (g) {
    A.GetBlock(ipg, ipm, m_Arm);
  }
  m_Acm.Resize(N);
  for (int i = 0, ip = ipc; i < N; ++i, ip += 6) {
    A.GetBlock(ip, ipm, m_Acm[i]);
  }
  A.GetBlock(ipm, ipm, m_Amm);
  b.GetBlock(ipm, m_bm);
}

void Joint::GetPriorEquation(Matrix::X *A, Vector::X *b, const bool symmetric) const {
  const int N = static_cast<int>(m_iKFs.size());
  const bool g = m_Zps.Size() != N;
  const int Npg = g ? 2 : 0;

  const int Npc = N * 6, Np = Npg + Npc + 9;
  const int ipg = 0, ipc = Npg, ipm = ipc + Npc;
  Pose::GetPriorEquation(A, b, symmetric);
  A->Resize(Np, Np, symmetric, true);
  if (b) {
    b->Resize(Np, true);
  }
  if (g) {
    A->SetBlock(ipg, ipm, m_Arm);
  }
  for (int i = 0, ip = ipc; i < N; ++i, ip += 6) {
    A->SetBlock(ip, ipm, m_Acm[i]);
  }
  A->SetBlock(ipm, ipm, m_Amm);
  if (b) {
    b->SetBlock(ipm, m_bm);
  }
}

bool Joint::GetPriorMeasurement(const Element::T w, Matrix::X *S, Vector::X *x,
                                Element::T *xTb, const Element::T *eps) const {
  GetPriorEquation(S, x, false);
  if (x) {
    Vector::X b;
    if (xTb) {
      const int Np = x->Size(), NpC = SIMD::Ceil<Element::T>(Np);
      b.Bind(x->Data() + NpC, Np);
      b.Copy(*x);
    }
    if (!S->SolveLDL(*x, eps)) {
      return false;
    }
    x->MakeMinus();
    S->InverseLDL(eps, true);
    if (xTb) {
      *xTb = x->Dot(b);
    }
  } else {
    if (!S->InverseLDL(eps)) {
      return false;
    }
  }

  if (w != 1) {
    *S *= w;
    if (xTb) {
      *xTb /= w;
    }
  }
  return true;
}

bool Joint::GetPriorMeasurement(const float w, LA::AlignedMatrixXf *S, LA::AlignedVectorXf *x,
                                float *xTb, AlignedVector<float> *work, const float *eps) const {
  if (Pose::Invalid() || Motion::Invalid()) {
    return false;
  }
  const int N = static_cast<int>(m_iKFs.size());
  const bool g = m_Zps.Size() != N;
  const int Npg = g ? 2 : 0;

  const int Npgc = Npg + N * 6, Np = Npgc + 9;
  const int ipg1 = 0, ipc = Npg, ipm = Npgc;

  Matrix::X _S;
  Vector::X _x, _eps;
  Element::T _xTb;
  const int _Nx = xTb ? (SIMD::Ceil<Element::T>(Np) + Np) : Np;
  work->Resize((S->BindSize(Np, Np, true) + (x ? x->BindSize(Np) : 0) +
                _S.BindSize(Np, Np) + (x ? _x.BindSize(_Nx) : 0) +
               (eps ? _eps.BindSize(Np) : 0)) / sizeof(float));
  S->Bind(work->Data(), Np, Np, true);
  if (x) {
    x->Bind(S->BindNext(), Np);
  }
  _S.Bind(x ? x->BindNext() : S->BindNext(), Np, Np);
  if (x) {
    _x.Bind(_S.BindNext(), _Nx);
  }
  if (eps) {
    _eps.Bind(x ? _x.BindNext() : _S.BindNext(), Np);
    if (g) {
      LA::Vector2<Element::T> *eg = (LA::Vector2<Element::T> *) (_eps.Data() + ipg1);
      //eg->MakeZero();
      eg->Set(eps + 3);
    }
    LA::Vector6<Element::T> *ecs = (LA::Vector6<Element::T> *) (_eps.Data() + ipc);
    for (int i = 0; i < N; ++i) {
      ecs[i].Set(eps);
    }
    LA::Vector9<Element::T> *em = (LA::Vector9<Element::T> *) (_eps.Data() + ipm);
    em->Set(eps + 6);
  }
  if (!GetPriorMeasurement(static_cast<Element::T>(w), &_S, x ? &_x : NULL,
                           xTb ? &_xTb : NULL, eps ? _eps.Data() : NULL)) {
    S->Resize(0, 0, true);
    if (x) {
      x->Resize(0);
      if (xTb) {
        *xTb = 0.0f;
      }
    }
    return false;
  }
  if (x) {
    _x.Get(*x);
    if (xTb) {
      *xTb = static_cast<float>(_xTb);
    }
  }
  _S.Get(*S);
  return true;
}

bool Joint::Invertible(AlignedVector<float> *work, const float *eps) const {
  const int N = static_cast<int>(m_iKFs.size());
  const bool g = m_Zps.Size() != N;
  const int Npg = g ? 2 : 0;

  const int Npgc = Npg + N * 6, Np = Npgc + 9;
  const int ipg = 0, ipc = Npg, ipm = Npgc;

  Matrix::X S;
  Vector::X _work, _eps;
  work->Resize((S.BindSize(Np, Np, true) + _work.BindSize(Np) +
               (eps ? _eps.BindSize(Np) : 0)) / sizeof(float));
  S.Bind(work->Data(), Np, Np, true);
  _work.Bind(S.BindNext(), Np);
  if (eps) {
    _eps.Bind(_work.BindNext(), Np);
    if (g) {
      LA::Vector2<Element::T> *eg = (LA::Vector2<Element::T> *) (_eps.Data() + ipg);
      //eg->MakeZero();
      eg->Set(eps + 3);
    }
    LA::Vector6<Element::T> *ecs = (LA::Vector6<Element::T> *) (_eps.Data() + ipc);
    for (int i = 0; i < N; ++i) {
      ecs[i].Set(eps);
    }
    LA::Vector9<Element::T> *em = (LA::Vector9<Element::T> *) (_eps.Data() + ipm);
    em->Set(eps + 6);
  }
  GetPriorEquation(&S);
  const int rank = S.RankLDL(&_work, eps ? _eps.Data() : NULL);
  return rank == Np;
}

bool Joint::PropagateLF(const Rigid3D &Tr, const Camera &C,
                        const IMU::Delta::Factor::Auxiliary::RelativeLF &A,
                        AlignedVector<float> *work, const float *eps) {
  const int N = static_cast<int>(m_iKFs.size()), Nk = N - 1;
  SetPose(Tr, Nk, C.m_T);
  SetMotion(C.m_T, C.m_v, C.m_ba, C.m_bw);

  Matrix::X _A;
  Vector::X _b, _work, _eps;
  const int Npgk = 2 + Nk * 6, Npcm = 15, Np1 = Npgk + Npcm, Np2 = Np1 + Npcm;
  work->Resize((_work.BindSize(Np2) + _A.BindSize(Np2, Np2, true) + _b.BindSize(Np2) +
               (eps ? _eps.BindSize(Npcm) : 0)) / sizeof(float));
  _work.Bind(work->Data(), Np2);
  _A.Bind(_work.BindNext(), Np2, Np2, true);
  _b.Bind(_A.BindNext(), Np2);
  if (eps) {
    _eps.Bind(_b.BindNext(), Npcm);
    _eps.Set(eps);
  }
  GetPriorEquation(&_A, &_b);
  const int ipg = 0, ipc1 = Npgk, ipc2 = ipc1 + Npcm;
  _A.InsertZero(ipc2, Npcm, NULL);
  _b.InsertZero(ipc2, Npcm, NULL);
  //_A.Resize(Np2, Np2, true, true);
  //_b.Resize(Np2, true);
  _A.IncreaseBlockDiagonal(ipg, A.m_Agg);
  _b.IncreaseBlock(ipg, A.m_bg);
  for (int i = 0, ip = ipc1; i < 10; ++i, ip += 3) {
    _A.IncreaseBlock(0, ip, A.m_Agc[i]);
    _b.IncreaseBlock(ip, A.m_b[i]);
  }
  for (int i = 0, ip = ipc1, k = 0; i < 10; ++i, ip += 3) {
    for (int j = i, jp = ip; j < 10; ++j, jp += 3, ++k) {
      _A.IncreaseBlock(ip, jp, A.m_A[k]);
    }
  }
  //UT::PrintSeparator();
  //_b.Print(true);
  const bool scc = _A.MarginalizeLDL(Npgk, Npcm, _b, &_work, eps ? _eps.Data() : NULL);
  //UT::PrintSeparator();
  //_b.Print(true);
  SetPriorEquation(_A, _b);
  return scc;
}

bool Joint::PropagateLF(const IMU::Delta::Factor::Auxiliary::RelativeLF &A, LA::AlignedVectorXf *x,
                        AlignedVector<float> *work, const float *eps) const {
  const int N1 = static_cast<int>(m_iKFs.size()), N2 = N1 + 1, Nk = N1 - 1;
  Matrix::X _A;
  Vector::X _b, _eps;
  const int Npgk = 2 + Nk * 6, Npcm = 15, Np1 = Npgk + Npcm, Np2 = Np1 + Npcm;
  work->Resize((x->BindSize(Np2) + _A.BindSize(Np2, Np2) + _b.BindSize(Np2) +
               (eps ? _eps.BindSize(Np2) : 0)) / sizeof(float));
  x->Bind(work->Data(), Np2);
  _A.Bind(x->BindNext(), Np2, Np2);
  _b.Bind(_A.BindNext(), Np2);
  if (eps) {
    _eps.Bind(_b.BindNext(), Np2);
  }
  GetPriorEquation(&_A, &_b, false);
  const int ipg = 0, ipc1 = Npgk, ipc2 = ipc1 + Npcm;
  _A.InsertZero(ipc2, Npcm, NULL);
  _b.InsertZero(ipc2, Npcm, NULL);
  //_A.Resize(Np2, Np2, true, true);
  //_b.Resize(Np2, true);
  _A.IncreaseBlockDiagonal(ipg, A.m_Agg);
  _b.IncreaseBlock(ipg, A.m_bg);
  for (int i = 0, ip = ipc1; i < 10; ++i, ip += 3) {
    _A.IncreaseBlock(0, ip, A.m_Agc[i]);
    _b.IncreaseBlock(ip, A.m_b[i]);
  }
  for (int i = 0, ip = ipc1, k = 0; i < 10; ++i, ip += 3) {
    for (int j = i, jp = ip; j < 10; ++j, jp += 3, ++k) {
      _A.IncreaseBlock(ip, jp, A.m_A[k]);
    }
  }
  if (eps) {
    LA::Vector2<Element::T> *eg = (LA::Vector2<Element::T> *) (_eps.Data() + ipg);
    //eg->MakeZero();
    eg->Set(eps + 3);
    LA::Vector6<Element::T> *eks = (LA::Vector6<Element::T> *) (eg + 1);
    for (int i = 0; i < Nk; ++i) {
      eks[i].Set(eps);
    }
    LA::Vector6<Element::T> *ec1 = eks + Nk;
    LA::Vector9<Element::T> *em1 = (LA::Vector9<Element::T> *) (ec1 + 1);
    LA::Vector6<Element::T> *ec2 = (LA::Vector6<Element::T> *) (em1 + 1);
    LA::Vector9<Element::T> *em2 = (LA::Vector9<Element::T> *) (ec2 + 1);
    ec1->Set(eps);
    ec2->Set(eps);
    em1->Set(eps + 6);
    em2->Set(eps + 6);
  }
  if (!_A.SolveLDL(_b, eps ? _eps.Data() : NULL)) {
    x->Resize(0);
    return false;
  }
  _b.MakeMinus();
  //x->Copy(_b);
  _b.GetBlock(0, *x);

  return true;
}

bool Joint::PropagateKF(const Rigid3D &Tr, const Camera &C,
                        const IMU::Delta::Factor::Auxiliary::RelativeKF &A,
                        AlignedVector<float> *work, const float *eps) {
  Matrix::X _A;
  Vector::X _b, _work, _eps;
  const int Npg = 2, Npc = 6, Npm = 9, Np1 = Npg + Npm, Npcm = Npc + Npm, Np2 = Np1 + Npcm;
  work->Resize((_work.BindSize(Np2) + _A.BindSize(Np2, Np2, true) + _b.BindSize(Np2) +
               (eps ? _eps.BindSize(Npm) : 0)) / sizeof(float));
  _work.Bind(work->Data(), Np2);
  _A.Bind(_work.BindNext(), Np2, Np2, true);
  _b.Bind(_A.BindNext(), Np2);
  if (eps) {
    _eps.Bind(_b.BindNext(), Npm);
    _eps.Set(eps + Npc);
  }
  GetPriorEquation(&_A, &_b);

  const int ipg = 0, ipc1 = Npg, ipc2 = ipc1 + Npm;
  _A.InsertZero(ipc2, Npcm, NULL);
  _b.InsertZero(ipc2, Npcm, NULL);
  //_A.Resize(Np2, Np2, true, true);
  //_b.Resize(Np2, true);
  _A.IncreaseBlockDiagonal(ipg, A.m_Agg);
  _b.IncreaseBlock(ipg, A.m_bg);
  for (int i = 0, ip = ipc1; i < 8; ++i, ip += 3) {
    _A.IncreaseBlock(0, ip, A.m_Agc[i]);
    _b.IncreaseBlock(ip, A.m_bc[i]);
  }
  for (int i = 0, ip = ipc1, k = 0; i < 8; ++i, ip += 3) {
    for (int j = i, jp = ip; j < 8; ++j, jp += 3, ++k) {
      _A.IncreaseBlock(ip, jp, A.m_Ac[k]);
    }
  }
  const bool scc = _A.MarginalizeLDL(Npg, Npm, _b, &_work, eps ? _eps.Data() : NULL);

  m_iKFs.resize(1, INT_MAX);
  //m_Zps.Resize(1);
  const Rotation3D RrT = m_Zps.Back();
  m_Zps.Resize(2);
  m_Zps.Back() = RrT;
  SetPose(Tr, 0, C.m_T);
  SetMotion(C.m_T, C.m_v, C.m_ba, C.m_bw);
  SetPriorEquation(_A, _b);

  return scc;
}

bool Joint::PropagateKF(const IMU::Delta::Factor::Auxiliary::RelativeKF &A, LA::AlignedVectorXf *x,
                        AlignedVector<float> *work, const float *eps) const {
  Matrix::X _A;
  Vector::X _b, _eps;
  const int Npg = 2, Npc = 6, Npm = 9, Np1 = Npg + Npm, Npcm = Npc + Npm, Np2 = Np1 + Npcm;
  work->Resize((x->BindSize(Np2) + _A.BindSize(Np2, Np2) + _b.BindSize(Np2) +
               (eps ? _eps.BindSize(Np2) : 0)) / sizeof(float));
  x->Bind(work->Data(), Np2);
  _A.Bind(x->BindNext(), Np2, Np2);
  _b.Bind(_A.BindNext(), Np2);
  if (eps) {
    _eps.Bind(_b.BindNext(), Np2);
  }
  GetPriorEquation(&_A, &_b, false);

  const int ipg = 0, ipc1 = Npg, ipc2 = ipc1 + Npm;
  _A.InsertZero(ipc2, Npcm, NULL);
  _b.InsertZero(ipc2, Npcm, NULL);
  //_A.Resize(Np2, Np2, true, true);
  //_b.Resize(Np2, true);
  _A.IncreaseBlockDiagonal(ipg, A.m_Agg);
  _b.IncreaseBlock(ipg, A.m_bg);
  for (int i = 0, ip = ipc1; i < 8; ++i, ip += 3) {
    _A.IncreaseBlock(0, ip, A.m_Agc[i]);
    _b.IncreaseBlock(ip, A.m_bc[i]);
  }
  for (int i = 0, ip = ipc1, k = 0; i < 8; ++i, ip += 3) {
    for (int j = i, jp = ip; j < 8; ++j, jp += 3, ++k) {
      _A.IncreaseBlock(ip, jp, A.m_Ac[k]);
    }
  }
  if (eps) {
    LA::Vector2<Element::T> *eg = (LA::Vector2<Element::T> *) (_eps.Data() + ipg);
    //eg->MakeZero();
    eg->Set(eps + 3);
    LA::Vector9<Element::T> *em1 = (LA::Vector9<Element::T> *) (_eps.Data() + ipc1);
    LA::Vector6<Element::T> *ec2 = (LA::Vector6<Element::T> *) (em1 + 1);
    LA::Vector9<Element::T> *em2 = (LA::Vector9<Element::T> *) (ec2 + 1);
    ec2->Set(eps);
    em1->Set(eps + 6);
    em2->Set(eps + 6);
  }
  if (!_A.SolveLDL(_b, eps ? _eps.Data() : NULL)) {
    x->Resize(0);
    return false;
  }
  _b.MakeMinus();
  //x->Copy(_b);
  _b.GetBlock(0, *x);

  return true;
}

bool Joint::GetPriorPose(const int iKF, Pose *Zp, AlignedVector<float> *work, const float *eps) const {
  Zp->m_iKFr = m_iKFr;
  Zp->m_iKFs = m_iKFs;
  Zp->m_iKFs.back() = iKF;
  Zp->m_Zps.Set(m_Zps);

  Matrix::X A;
  Vector::X b, _work, _eps;
  const int N = static_cast<int>(m_iKFs.size());
  const int Npg = m_Zps.Size() == N ? 0 : 2, Npc = N * 6, Np = Npg + Npc + 9;
  work->Resize((A.BindSize(Np, Np, true) + b.BindSize(Np) + _work.BindSize(Np) +
               (eps ? _eps.BindSize(15) : 0)) / sizeof(float));
  A.Bind(work->Data(), Np, Np, true);
  b.Bind(A.BindNext(), Np);
  _work.Bind(b.BindNext(), Np);
  if (eps) {
    _eps.Bind(_work.BindNext(), 15);
    _eps.Set(eps);
  }
  bool scc = true;
  GetPriorEquation(&A, &b);
  if (iKF == INT_MAX) {
    const int Nk = N - 1;
    Zp->m_iKFs.resize(Nk);
    //Zp->m_Zps.Resize(Nk);
    const Rotation3D RrT = Zp->m_Zps.Back();
    Zp->m_Zps.Resize(Nk + 1);
    Zp->m_Zps.Back() = RrT;
    scc = A.MarginalizeLDL(Np - 15, 15, b, &_work, eps ? _eps.Data() : NULL);
  } else {
    scc = A.MarginalizeLDL(Np - 9, 9, b, &_work, eps ? _eps.Data() + 6 : NULL);
  }
  Zp->SetPriorEquation(A, b);

  return scc;
}

bool Joint::GetPriorMotion(Motion *Zp, AlignedVector<float> *work, const float *eps) const {
  Zp->m_v = m_v;
  Zp->m_ba = m_ba;
  Zp->m_bw = m_bw;
  Matrix::X A;
  Vector::X b, _work, _eps;
  const int N = static_cast<int>(m_iKFs.size());
  const bool g = m_Zps.Size() != N;
  const int Npg = g ? 2 : 0;

  const int Npgc = Npg + N * 6, Np = Npgc + 9;
  const int ipg = 0;
  work->Resize((A.BindSize(Np, Np, true) + b.BindSize(Np) + _work.BindSize(Np) +
               (eps ? _eps.BindSize(Np) : 0)) / sizeof(float));
  A.Bind(work->Data(), Np, Np, true);
  b.Bind(A.BindNext(), Np);
  _work.Bind(b.BindNext(), Np);
  if (eps) {
    _eps.Bind(_work.BindNext(), Np);
    if (g) {
      LA::Vector2<Element::T> *eg = (LA::Vector2<Element::T> *) (_eps.Data() + ipg);
      //eg->MakeZero();
      eg->Set(eps + 3);
    }
    LA::Vector6<Element::T> *ecs = (LA::Vector6<Element::T> *) (_eps.Data() + Npg);
    for (int i = 0; i < N; ++i) {
      ecs[i].Set(eps);
    }
  }
  GetPriorEquation(&A, &b);
  const bool scc = A.MarginalizeLDL(0, Npgc, b, &_work, eps ? _eps.Data() : NULL);
  A.GetBlock(0, 0, Zp->m_Amm);
  b.GetBlock(0, Zp->m_bm);
  if (!scc) {
    return false;
  }
  LA::AlignedMatrix9x9f Amm = Zp->m_Amm;
  const int rank = Amm.RankLDL(eps ? eps + 6 : NULL);
  return rank == 9;
}

}  // namespace CameraPrior
