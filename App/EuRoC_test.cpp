#include "EuRoC_test.h"

#define SHOW_POSE
int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InstallFailureSignalHandler();
  if (FLAGS_imgs_folder.empty()) {
    google::ShowUsageWithFlags(argv[0]);
    return -1;
  }
  
  // Load IMU samples to predict OF point locations
  std::list<XP::ImuData> imu_samples;
  std::string imu_file = FLAGS_imgs_folder + "/mav0/imu0/data.csv";
  uint64_t offset_ts_ns;
  if (load_imu_data(imu_file, &imu_samples, offset_ts_ns) > 0) {
    std::cout << "Load imu data. Enable OF prediciton with gyro\n";
  } else {
    std::cout << "Cannot load imu data.\n";
    return -1;
  }
  
  vector<string> img_file_paths;
  vector<string> iba_dat_file_paths;

  constexpr int reserve_num = 5000;
  img_file_paths.reserve(reserve_num);
  fs::path p(FLAGS_imgs_folder + "/mav0/cam0");
  if (!fs::is_directory(p)) {
    LOG(ERROR) << p << " is not a directory";
    return -1;
  }
  if (!fs::is_directory(FLAGS_imgs_folder + "/dat")) {
    fs::create_directories(FLAGS_imgs_folder + "/dat");
  }
  vector<string> limg_name, rimg_name;
  load_image_data(FLAGS_imgs_folder, limg_name, rimg_name);
  for (int i=0; i<limg_name.size(); i++) {
    string l_png = p.string() + "/data/" + limg_name[i];
    img_file_paths.push_back(l_png);
    iba_dat_file_paths.push_back(FLAGS_imgs_folder + "/dat/" + limg_name[i] + ".dat");
  }

  if (img_file_paths.size() == 0) {
    LOG(ERROR) << "No image files for detection";
    return -1;
  }
  XP::DuoCalibParam duo_calib_param;
  try {
    load_asl_calib(FLAGS_imgs_folder, duo_calib_param);
  } catch (...){
    LOG(ERROR) << "Load calibration file error";
    return -1;
  }
  // Create masks based on FOVs computed from intrinsics
  std::vector<cv::Mat_<uchar> > masks(2);
  for (int lr = 0; lr < 2; ++lr) {
    float fov;
    if (XP::generate_cam_mask(duo_calib_param.Camera.cv_camK_lr[lr],
                              duo_calib_param.Camera.cv_dist_coeff_lr[lr],
                              duo_calib_param.Camera.img_size,
                              &masks[lr],
                              &fov)) {
      std::cout << "camera " << lr << " fov: " << fov << " deg\n";
    }
  }
  
  // Adjust end image index for detection
  if (FLAGS_end_idx < 0 || FLAGS_end_idx > img_file_paths.size()) {
    FLAGS_end_idx = img_file_paths.size();
  }

  FLAGS_start_idx = std::max(0, FLAGS_start_idx);
  // remove all frames before the first IMU data
  while (FLAGS_start_idx < FLAGS_end_idx && get_timestamp_from_img_name(img_file_paths[FLAGS_start_idx], offset_ts_ns) <= imu_samples.front().time_stamp)
    FLAGS_start_idx++;

  XP::FeatureTrackDetector feat_track_detector(FLAGS_ft_len,
                                               FLAGS_ft_droprate,
                                               !FLAGS_not_use_fast,
                                               FLAGS_uniform_radius,
                                               duo_calib_param.Camera.img_size);
  const Eigen::Matrix4f T_Cl_Cr =
      duo_calib_param.Camera.D_T_C_lr[0].inverse() * duo_calib_param.Camera.D_T_C_lr[1];
      
#ifdef SHOW_POSE
  XP::PoseViewer pose_viewer;
  pose_viewer.set_clear_canvas_before_draw(true);
#endif //SHOW_POSE
  
  IBA::Solver solver;
  Eigen::Vector3f last_position = Eigen::Vector3f::Zero();
  float travel_dist = 0.f;
  if (FLAGS_save_feature) {
    IBA::SaveCalibration(FLAGS_imgs_folder + "/calibration.dat", to_iba_calibration(duo_calib_param));
  }
  solver.Create(to_iba_calibration(duo_calib_param),
                257,
                IBA_VERBOSE_NONE,
                IBA_DEBUG_NONE,
                257,
                FLAGS_iba_param_path,
                "" /* iba directory */);
#ifdef SHOW_POSE
  solver.SetCallbackLBA([&](const int iFrm, const float ts) {
#ifndef __DUO_VIO_TRACKER_NO_DEBUG__
    VLOG(1) << "===== start ibaCallback at ts = " << ts;
#endif
    // as we may be able to send out information directly in the callback arguments
    IBA::SlidingWindow sliding_window;
    solver.GetSlidingWindow(&sliding_window);
    const IBA::CameraIMUState& X = sliding_window.CsLF.back();
    const IBA::CameraPose& C = X.C;
    Eigen::Matrix4f W_vio_T_S = Eigen::Matrix4f::Identity();  // W_vio_T_S
    for (int i = 0; i < 3; ++i) {
      W_vio_T_S(i, 3) = C.p[i];
      for (int j = 0; j < 3; ++j) {
        W_vio_T_S(i, j) = C.R[j][i];  // C.R is actually R_SW
      }
    }
    Eigen::Matrix<float, 9, 1> speed_and_biases;
    for (int i = 0; i < 3; ++i) {
      speed_and_biases(i) = X.v[i];
      speed_and_biases(i + 3) = X.ba[i];
      speed_and_biases(i + 6) = X.bw[i];
    }
    Eigen::Vector3f cur_position = W_vio_T_S.topRightCorner(3, 1);
    travel_dist += (cur_position - last_position).norm();
    last_position = cur_position;
    pose_viewer.addPose(W_vio_T_S, speed_and_biases, travel_dist);
  });
#endif // SHOW_POSE
  solver.Start();
  
  float prev_time_stamp = 0.0f;
  // load previous image
  std::vector<cv::KeyPoint> pre_image_key_points;
  cv::Mat pre_image_features;
  for (int it_img = FLAGS_start_idx; it_img < FLAGS_end_idx; ++it_img) {
    VLOG(0) << " start detection at ts = " << fs::path(img_file_paths[it_img]).stem().string();
    auto read_img_start = std::chrono::high_resolution_clock::now();
    cv::Mat img_in_raw;
    img_in_raw = cv::imread(img_file_paths[it_img], CV_LOAD_IMAGE_GRAYSCALE);
    CHECK_EQ(img_in_raw.rows, duo_calib_param.Camera.img_size.height);
    CHECK_EQ(img_in_raw.cols, duo_calib_param.Camera.img_size.width);
    cv::Mat img_in_smooth;
    cv::blur(img_in_raw, img_in_smooth, cv::Size(3, 3));
    if (img_in_smooth.rows == 0) {
      LOG(ERROR) << "Cannot load " << img_file_paths[it_img];
      return -1;
    }
    // get timestamp from image file name (s)
    const float time_stamp = get_timestamp_from_img_name(img_file_paths[it_img], offset_ts_ns) + timedelay;
    std::vector<cv::KeyPoint> key_pnts;
    cv::Mat orb_feat;
    cv::Mat pre_img_in_smooth;
    // load slave image
    cv::Mat slave_img_smooth;  // for visualization later
    std::vector<cv::KeyPoint> key_pnts_slave;
    cv::Mat orb_feat_slave;
    std::vector<XP::ImuData> imu_meas;

    // Get the imu measurements within prev_time_stamp and time_stamp to compute old_R_new
    imu_meas.reserve(10);
    for (auto it_imu = imu_samples.begin(); it_imu != imu_samples.end(); ) {
      if (it_imu->time_stamp < time_stamp) {
        imu_meas.push_back(*it_imu);
        it_imu++;
        imu_samples.pop_front();
      } else {
        break;
      }
    }
    VLOG(1) << "imu_meas size = " << imu_meas.size();
    VLOG(1) << "img ts prev -> curr " << prev_time_stamp << " -> " << time_stamp;
    if (imu_meas.size() > 0) {
      VLOG(1) << "imu ts prev -> curr " << imu_meas.front().time_stamp
              << " -> " << imu_meas.back().time_stamp;
    }

    // use optical flow  from the 1st frame
    if (it_img != FLAGS_start_idx) {
      CHECK(it_img >= 1);
      VLOG(1) << "pre_image_key_points.size(): " << pre_image_key_points.size();
      const int request_feat_num = FLAGS_max_num_per_grid * FLAGS_grid_row_num * FLAGS_grid_col_num;
      feat_track_detector.build_img_pyramids(img_in_smooth,
                                             XP::FeatureTrackDetector::BUILD_TO_CURR);
      if (imu_meas.size() > 1) {
        // Here we simply the transformation chain to rotation only and assume zero translation
        cv::Matx33f old_R_new;
        XP::XpQuaternion I_new_q_I_old;  // The rotation between the new {I} and old {I}
        for (size_t i = 1; i < imu_meas.size(); ++i) {
          XP::XpQuaternion q_end;
          XP::IntegrateQuaternion(imu_meas[i - 1].ang_v,
                                  imu_meas[i].ang_v,
                                  I_new_q_I_old,
                                  imu_meas[i].time_stamp - imu_meas[i - 1].time_stamp,
                                  &q_end);
          I_new_q_I_old = q_end;
        }
        Eigen::Matrix3f I_new_R_I_old = I_new_q_I_old.ToRotationMatrix();
        Eigen::Matrix4f I_T_C =
            duo_calib_param.Imu.D_T_I.inverse() * duo_calib_param.Camera.D_T_C_lr[0];
        Eigen::Matrix3f I_R_C = I_T_C.topLeftCorner<3, 3>();
        Eigen::Matrix3f C_new_R_C_old = I_R_C.transpose() * I_new_R_I_old * I_R_C;
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            old_R_new(j, i) = C_new_R_C_old(i, j);
          }
        }

        if (VLOG_IS_ON(1)) {
          XP::XpQuaternion C_new_q_C_old;
          C_new_q_C_old.SetFromRotationMatrix(C_new_R_C_old);
          VLOG(1) << "C_new_R_C_old = \n" << C_new_R_C_old;
          VLOG(1) << "ea =\n" << C_new_q_C_old.ToEulerRadians() * 180 / M_PI;
        }
        feat_track_detector.optical_flow_and_detect(masks[0],
                                                    pre_image_features,
                                                    pre_image_key_points,
                                                    request_feat_num,
                                                    FLAGS_pyra_level,
                                                    FLAGS_fast_thresh,
                                                    &key_pnts,
                                                    &orb_feat,
                                                    cv::Vec2f(0, 0),  // shift init pixels
                                                    &duo_calib_param.Camera.cv_camK_lr[0],
                                                    &duo_calib_param.Camera.cv_dist_coeff_lr[0],
                                                    &old_R_new);
      } else {
        feat_track_detector.optical_flow_and_detect(masks[0],
                                                    pre_image_features,
                                                    pre_image_key_points,
                                                    request_feat_num,
                                                    FLAGS_pyra_level,
                                                    FLAGS_fast_thresh,
                                                    &key_pnts,
                                                    &orb_feat);
      }
      feat_track_detector.update_img_pyramids();
      VLOG(1) << "after OF key_pnts.size(): " << key_pnts.size() << " requested # "
              << FLAGS_max_num_per_grid * FLAGS_grid_row_num * FLAGS_grid_col_num;
    } else {
      // first frame
      feat_track_detector.detect(img_in_smooth,
                                 masks[0],
                                 FLAGS_max_num_per_grid * FLAGS_grid_row_num * FLAGS_grid_col_num,
                                 FLAGS_pyra_level,
                                 FLAGS_fast_thresh,
                                 &key_pnts,
                                 &orb_feat);
      feat_track_detector.build_img_pyramids(img_in_smooth,
                                             XP::FeatureTrackDetector::BUILD_TO_PREV);
    }

    std::sort(key_pnts.begin(), key_pnts.end(), cmp_by_class_id);
    std::sort(key_pnts_slave.begin(), key_pnts_slave.end(), cmp_by_class_id);
    // push to IBA
    IBA::CurrentFrame CF;
    IBA::KeyFrame KF;
    create_iba_frame(key_pnts, key_pnts_slave, imu_meas, time_stamp, &CF, &KF);
    CF.fileName = img_file_paths[it_img];
    solver.PushCurrentFrame(CF, KF.iFrm == -1 ? nullptr : &KF);
    if (FLAGS_save_feature) {
      IBA::SaveCurrentFrame(iba_dat_file_paths[it_img], CF, KF);
    }
    pre_image_key_points = key_pnts;
    pre_image_features = orb_feat.clone();
#ifdef SHOW_POSE
    // show pose
    pose_viewer.displayTo("trajectory");
#endif // SHOW_POSE
    cv::waitKey(1);
    prev_time_stamp = time_stamp;
    
  }
  std::string temp_file = "/tmp/" + std::to_string(offset_ts_ns) + ".txt";
  solver.SaveCamerasGBA(temp_file, false /* append */, true /* pose only */);
//   solver.SaveCamerasLBA(temp_file, false /* append */, true /* pose only */);
  solver.Stop();
  solver.Destroy();

  // for comparsion with asl groundtruth
  convert_to_asl_timestamp(temp_file, FLAGS_gba_camera_save_path, offset_ts_ns);
  return 0;
}
