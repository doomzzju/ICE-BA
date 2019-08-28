/******************************************************************************
 * Copyright 2017 Baidu Robotic Vision Authors. All Rights Reserved.
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
#include "IBA/IBA.h"
#include "feature_utils.h"
#include "image_utils.h"
#include "xp_quaternion.h"
#include "param.h"  // calib
#include "basic_datatype.h"
#include "iba_helper.h"
#include "pose_viewer.h"
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <opencv2/core.hpp>

#include <algorithm>
#include <string>
#include <fstream>
#include <vector>

namespace fs = boost::filesystem;
using std::string;
using std::vector;
DEFINE_string(imgs_folder, "", "The folder containing l and r folders, and the calib.yaml");
DEFINE_int32(grid_row_num, 1, "Number of rows of detection grids");
DEFINE_int32(grid_col_num, 1, "Number of cols of detection grids");
DEFINE_int32(max_num_per_grid, 70, "Max number of points per grid");
DEFINE_double(feat_quality, 0.07, "Tomasi-Shi feature quality level");
DEFINE_double(feat_min_dis, 10, "Tomasi-Shi feature minimal distance");
DEFINE_bool(not_use_fast, false, "Whether or not use FAST");
DEFINE_int32(pyra_level, 2, "Total pyramid levels");
DEFINE_int32(start_idx, 0, "The image index of the first detection (from 0)");
DEFINE_int32(end_idx, -1, "The image index of the last detection");
DEFINE_double(uniform_radius, 40, "< 5 disables uniformaty enforcement");
DEFINE_int32(ft_len, 125, "The feature track length threshold when dropout kicks in");
DEFINE_double(ft_droprate, 0.05, "The drop out rate when a feature track exceeds ft_len");
DEFINE_bool(show_feat_only, false, "wether or not show detection results only");
DEFINE_int32(fast_thresh, 10, "FAST feature threshold (only meaningful if use_fast=true)");
DEFINE_double(min_feature_distance_over_baseline_ratio,
              4, "Used for slave image feature detection");
DEFINE_double(max_feature_distance_over_baseline_ratio,
              3000, "Used for slave image feature detection");
DEFINE_string(iba_param_path, "", "iba parameters path");
DEFINE_string(gba_camera_save_path, "", "Save the camera states to when finished");
DEFINE_bool(stereo, false, "monocular or stereo mode");
DEFINE_bool(save_feature, false, "Save features to .dat file");

uint64_t t_last_imu = 0;
double timedelay = 0.0;

size_t load_image_data(const string& image_folder,
                       std::vector<string> &limg_name,
                       std::vector<string> &rimg_name) {
  LOG(INFO) << "Loading " << image_folder;
  std::string l_path = image_folder + "/mav0/cam0/data.csv";
  std::string r_path = image_folder + "/mav0/cam1/data.csv";
  std::string l_img_prefix = image_folder + "/mav0/cam0/data/";
  std::string r_img_prefix = image_folder + "/mav0/cam1/data/";
  std::ifstream limg_file(l_path);
  std::ifstream rimg_file(r_path);
  if (!limg_file.is_open() /*|| !rimg_file.is_open()*/) {
    LOG(WARNING) << image_folder << " cannot be opened";
    return 0;
  }
  std::string line;
  std::string time;
  while (getline(limg_file,line)) {
    if (line[0] == '#')
      continue;
    std::istringstream is(line);
    int i = 0;
    while (getline(is, time, ',')){
      bool is_exist = boost::filesystem::exists(l_img_prefix + time + ".png");
      uint64_t t;
      std::stringstream ss;
      ss << time;
      ss >> t;
      bool is_valid = (t <= t_last_imu);
      if (i == 0 && is_exist && is_valid){
        limg_name.push_back(time + ".png");
//         rimg_name.push_back(time + ".png");
      }
      i++;
    }
  }
  limg_file.close();
  rimg_file.close();
  LOG(INFO)<< "loaded " << limg_name.size() << " images";
  return limg_name.size();
}

size_t load_imu_data(const string& imu_file_str,
                     std::list<XP::ImuData>* imu_samples_ptr,
                     uint64_t &offset_ts_ns) {
  CHECK(imu_samples_ptr != NULL);
  LOG(INFO) << "Loading " << imu_file_str;
  std::ifstream imu_file(imu_file_str.c_str());
  if (!imu_file.is_open()) {
    LOG(WARNING) << imu_file_str << " cannot be opened";
    return 0;
  }
  std::list<XP::ImuData>& imu_samples = *imu_samples_ptr;
  imu_samples.clear();
  // read imu data
  std::string line;
  std::string item;
  double c[6];
  uint64_t t;
  bool set_offset_time = false;
  while (getline(imu_file,line)) {
    if (line[0] == '#')
      continue;
    std::istringstream is(line);
    int i = 0;
    while (getline(is, item, ',')) {
      std::stringstream ss;
      ss << item;
      if (i == 0)
        ss >> t;
      else
        ss >> c[i-1];
      i++;
    }
    if (!set_offset_time) {
      set_offset_time = true;
      offset_ts_ns = t;
    }
    t_last_imu = t;
    XP::ImuData imu_sample;
    float _t_100us = (t - offset_ts_ns)/1e5;
    imu_sample.time_stamp = _t_100us/1e4;
    imu_sample.ang_v(0) = c[0];
    imu_sample.ang_v(1) = c[1];
    imu_sample.ang_v(2) = c[2];
    imu_sample.accel(0) = c[3];
    imu_sample.accel(1) = c[4];
    imu_sample.accel(2) = c[5];

    VLOG(3) << "accel " << imu_sample.accel.transpose()
            << " gyro " << imu_sample.ang_v.transpose();
    imu_samples.push_back(imu_sample);
  }
  imu_file.close();
  LOG(INFO)<< "loaded " << imu_samples.size() << " imu samples";
  return imu_samples.size();
}

void load_asl_calib(const std::string &asl_path,
                    XP::DuoCalibParam &calib_param) {
  std::string cam0_yaml = asl_path + "/mav0/cam0/sensor.yaml";
//   std::string cam1_yaml = asl_path + "/mav0/cam1/sensor.yaml";
  std::string imu0_yaml = asl_path + "/mav0/imu0/sensor.yaml";
  YAML::Node cam0_calib = YAML::LoadFile(cam0_yaml);
//   YAML::Node cam1_calib = YAML::LoadFile(cam1_yaml);
  YAML::Node imu0_calib = YAML::LoadFile(imu0_yaml);
  // intrinsics
  std::vector<float> v_float = cam0_calib["intrinsics"].as<std::vector<float>>();
  calib_param.Camera.cv_camK_lr[0] << v_float[0], 0, v_float[2],
      0, v_float[1], v_float[3],
      0, 0, 1;
  calib_param.Camera.cameraK_lr[0] << v_float[0], 0, v_float[2],
      0, v_float[1], v_float[3],
      0, 0, 1;
//   v_float = cam1_calib["intrinsics"].as<std::vector<float>>();
  calib_param.Camera.cv_camK_lr[1] << v_float[0], 0, v_float[2],
      0, v_float[1], v_float[3],
      0, 0, 1;
  calib_param.Camera.cameraK_lr[1] << v_float[0], 0, v_float[2],
      0, v_float[1], v_float[3],
      0, 0, 1;
  // distortion_coefficients
  std::vector<double> v_double = cam0_calib["distortion_coefficients"].as<std::vector<double>>();
  calib_param.Camera.cv_dist_coeff_lr[0] = (cv::Mat_<float>(8, 1) << static_cast<float>(v_double[0]), static_cast<float>(v_double[1]),
      static_cast<float>(v_double[2]), static_cast<float>(v_double[3]), 0.0, 0.0, 0.0, 0.0);
//   v_double = cam1_calib["distortion_coefficients"].as<std::vector<double>>();
  calib_param.Camera.cv_dist_coeff_lr[1] = (cv::Mat_<float>(8, 1) << static_cast<float>(v_double[0]), static_cast<float>(v_double[1]),
      static_cast<float>(v_double[2]), static_cast<float>(v_double[3]), 0.0, 0.0, 0.0, 0.0);
  
  //TBS
  v_double = cam0_calib["T_BS"]["data"].as<std::vector<double>>();
  Eigen::Matrix4d b_t_c0 = Eigen::Map<Eigen::Matrix<double, 4, 4, Eigen::RowMajor>>(&v_double[0]);
//   v_double = cam1_calib["T_BS"]["data"].as<std::vector<double>>();
  Eigen::Matrix4d b_t_c1 = Eigen::Map<Eigen::Matrix<double, 4, 4, Eigen::RowMajor>>(&v_double[0]);
  v_double = imu0_calib["T_BS"]["data"].as<std::vector<double>>();
  Eigen::Matrix4d b_t_i = Eigen::Map<Eigen::Matrix<double, 4, 4, Eigen::RowMajor>>(&v_double[0]);
  // ASL {B}ody frame is the IMU
  // {D}evice frame is the left camera
  Eigen::Matrix4d d_t_cam0 = Eigen::Matrix4d::Identity();
  Eigen::Matrix4d d_t_b = d_t_cam0 * b_t_c0.inverse();
  Eigen::Matrix4d d_t_cam1 = d_t_b * b_t_c1;
  Eigen::Matrix4d d_t_imu = d_t_b * b_t_i;
  calib_param.Camera.D_T_C_lr[0] = Eigen::Matrix4f::Identity();
  calib_param.Camera.D_T_C_lr[1] = d_t_cam1.cast<float>();
  // Image size
  std::vector<int> v_int = cam0_calib["resolution"].as<std::vector<int>>();
  calib_param.Camera.img_size = cv::Size(v_int[0], v_int[1]);
  // IMU
  calib_param.Imu.accel_TK = Eigen::Matrix3f::Identity();
  calib_param.Imu.accel_bias = Eigen::Vector3f::Zero();
  calib_param.Imu.gyro_TK = Eigen::Matrix3f::Identity();
  calib_param.Imu.gyro_bias = Eigen::Vector3f::Zero();
  calib_param.Imu.accel_noise_var = Eigen::Vector3f{0.0016, 0.0016, 0.0016};
  calib_param.Imu.angv_noise_var = Eigen::Vector3f{0.0001, 0.0001, 0.0001};
  calib_param.Imu.D_T_I = d_t_imu.cast<float>();
  calib_param.device_id = "ASL";
  calib_param.sensor_type = XP::DuoCalibParam::SensorType::UNKNOWN;

  calib_param.initUndistortMap(calib_param.Camera.img_size);
  
  timedelay = cam0_calib["time-delay"].as<std::vector<double>>()[0];
  std::cout << "timedelay " << timedelay << std::endl;
}

float get_timestamp_from_img_name(const string& img_name,
                                  uint64_t offset_ns) {
  string ts_ns_string = fs::path(img_name).stem().string();
  int64_t offset_t = boost::lexical_cast<uint64_t>(ts_ns_string) - offset_ns;
  int64_t t = offset_t/1e5;
  return static_cast<float>(t)/1e4;
}

bool convert_to_asl_timestamp(const string& file_in,
                              const string& file_out,
                              uint64_t offset_ns) {
  FILE *fp_in = fopen(file_in.c_str(), "r");
  FILE *fp_out = fopen(file_out.c_str(), "w+");
  if (!fp_in || !fp_out) {
    LOG(ERROR) << "convert to asl timestamp error";
    return false;
  }
  float t;
  float x, y, z;
  float qx, qy, qz, qw;
  while (fscanf(fp_in, "%f %f %f %f %f %f %f %f", &t, &x, &y, &z, &qx, &qy, &qz, &qw) == 8) {
    double t_s = t + static_cast<double>(offset_ns*1e-9);
    fprintf(fp_out, "%lf %f %f %f %f %f %f %f\n", t_s, x, y, z, qx, qy, qz, qw);
  }
  fclose(fp_in);
  fclose(fp_out);
}

inline bool cmp_by_class_id(const cv::KeyPoint& lhs, const cv::KeyPoint& rhs)  {
  return lhs.class_id < rhs.class_id;
}

template <typename T>
void InitPOD(T& t) {
  memset(&t, 0, sizeof(t));
}

bool create_iba_frame(const vector<cv::KeyPoint>& kps_l,
                      const vector<cv::KeyPoint>& kps_r,
                      const vector<XP::ImuData>& imu_samples,
                      const float rig_time,
                      IBA::CurrentFrame* ptrCF, IBA::KeyFrame* ptrKF) {
  CHECK(std::is_sorted(kps_l.begin(), kps_l.end(), cmp_by_class_id));
  CHECK(std::is_sorted(kps_r.begin(), kps_r.end(), cmp_by_class_id));
  CHECK(std::includes(kps_l.begin(), kps_l.end(), kps_r.begin(), kps_r.end(), cmp_by_class_id));

  // IBA will handle *unknown* initial depth values
  IBA::Depth kUnknownDepth;
  kUnknownDepth.d = 0.0f;
  kUnknownDepth.s2 = 0.0f;
  static int last_added_point_id = -1;
  static int iba_iFrm = 0;
  auto kp_it_l = kps_l.cbegin(), kp_it_r = kps_r.cbegin();

  IBA::CurrentFrame& CF = *ptrCF;
  IBA::KeyFrame& KF = *ptrKF;

  CF.iFrm = iba_iFrm;
  InitPOD(CF.C); // needed to ensure the dumped frame deterministic even for unused field
  CF.C.C.R[0][0] = CF.C.v[0] = CF.C.ba[0] = CF.C.bw[0] = FLT_MAX;
  // MapPointMeasurement, process in ascending class id, left camera to right
  // Note the right keypoints is a subset of the left ones
  IBA::MapPointMeasurement mp_mea;
  InitPOD(mp_mea);
  mp_mea.x.S[0][0] = mp_mea.x.S[1][1] = 1.f;
  mp_mea.x.S[0][1] = mp_mea.x.S[1][0] = 0.f;
  for (; kp_it_l != kps_l.cend() && kp_it_l->class_id <= last_added_point_id; ++kp_it_l) {
    mp_mea.idx = kp_it_l->class_id;
    mp_mea.x.x[0] = kp_it_l->pt.x;
    mp_mea.x.x[1] = kp_it_l->pt.y;
    CF.zs.push_back(mp_mea);
    if (kp_it_r != kps_r.cend() && kp_it_r->class_id == kp_it_l->class_id) {
      mp_mea.x.x[0] = kp_it_r->pt.x;
      mp_mea.x.x[1] = kp_it_r->pt.y;
      CF.zs.push_back(mp_mea);
      ++kp_it_r;
    }
  }
  std::transform(imu_samples.begin(), imu_samples.end(), std::back_inserter(CF.us), XP::to_iba_imu);
  CF.t = rig_time;
  CF.d = kUnknownDepth;
  bool need_new_kf = std::distance(kp_it_l, kps_l.end()) >= 20 || CF.zs.size() < 20;
  if (std::distance(kp_it_l, kps_l.end()) == 0)
    need_new_kf = false;
  if (!need_new_kf) KF.iFrm = -1;
  else
    LOG(INFO) << "new keyframe " << KF.iFrm;

  if (!need_new_kf) {
    KF.iFrm = -1;
    //  to make it deterministic
    InitPOD(KF.C);
    InitPOD(KF.d);
  } else {
    KF.iFrm = CF.iFrm;
    KF.C = CF.C.C;
    // MapPointMeasurement, duplication of CF
    KF.zs = CF.zs;
    // MapPoint
    for(; kp_it_l != kps_l.cend(); ++kp_it_l) {
      IBA::MapPoint mp;
      InitPOD(mp.X);
      mp.X.idx = kp_it_l->class_id;
      mp.X.X[0] = FLT_MAX;
      mp_mea.iFrm = iba_iFrm;
      mp_mea.x.x[0] = kp_it_l->pt.x;
      mp_mea.x.x[1] = kp_it_l->pt.y;
      mp.zs.push_back(mp_mea);

      if (kp_it_r != kps_r.cend() && kp_it_r->class_id == kp_it_l->class_id) {
        mp_mea.x.x[0] = kp_it_r->pt.x;
        mp_mea.x.x[1] = kp_it_r->pt.y;
        mp.zs.push_back(mp_mea);
        kp_it_r++;
      } else {
        LOG(WARNING) << "add new feature point " << kp_it_l->class_id << " only found in left image";
      }
      KF.Xs.push_back(mp);
    }
    last_added_point_id = std::max(KF.Xs.back().X.idx, last_added_point_id);
    KF.d = kUnknownDepth;
  }
  ++iba_iFrm;
  return true;
}