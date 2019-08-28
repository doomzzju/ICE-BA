#include "ViewerPangolin.h"

void ViewerPangolin::Initialize(IBA::Solver* solver, const int img_width, const int img_height)
{
  MT::Thread::Initialize(1, 6, "ViewerPangolin");
  
  m_ViewpointX = 0;
  m_ViewpointY = -0.7;
  m_ViewpointZ = -1.8;
  m_ViewpointF = 500;

  m_KeyFrameSize = 0.05;
  m_KeyFrameLineWidth = 1.0;        
  m_CameraSize = 0.1;
  m_CameraLineWidth = 2.0;
  m_PointSize = 4;

  m_bFinished = false;

  m_img_width = img_width;
  m_img_height = img_height;

}

void ViewerPangolin::Reset()
{
  MT::Thread::Reset();
  // TODO
}

void ViewerPangolin::Run()
{
  m_bFinished = false;
  pangolin::CreateWindowAndBind("ICE-BA Viewer", 1920, 1440);

  // 3D Mouse handler requires depth testing to be enabled
  glEnable(GL_DEPTH_TEST);

  // Issue specific OpenGl we might need
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  pangolin::CreatePanel("menu").SetBounds(0.0, 1.0, 0.0, pangolin::Attach::Pix(175));
  pangolin::Var<bool> menuFollowCamera("menu.Follow Camera", false, true);
  pangolin::Var<bool> menuShowSlidingWindow("menu.Show Sliding Window", true, true);
  pangolin::Var<bool> menuShowLocalFrame("menu.Show Local Frame", true, true);
  pangolin::Var<bool> menuShowLocalPoints("menu.Show Local Points", true, true);
  pangolin::Var<bool> menuReset("menu.Reset", false, false);

  // Define Camera Render Object (for view / scene browsing)
  pangolin::OpenGlRenderState s_cam(
      pangolin::ProjectionMatrix(1920, 1440, m_ViewpointF, m_ViewpointF, 960, 720, 0.1, 1000),
	  pangolin::ModelViewLookAt(m_ViewpointX, m_ViewpointY, m_ViewpointZ, 0, 0, 0, 0.0, -1.0, 0.0)
  );

  // Add named OpenGL viewport to window and provide 3D Handler
  pangolin::View& d_cam = pangolin::CreateDisplay()
      .SetBounds(0.0, 1.0, pangolin::Attach::Pix(175), 1.0, -1920.0f / 1440.0f)
	  .SetHandler(new pangolin::Handler3D(s_cam));

  pangolin::GlTexture texVideo(m_img_width, m_img_height, GL_RGB, false, 0, GL_RGB, GL_UNSIGNED_BYTE);
  //pangolin::View& d_video = pangolin::Display("imgVideo").SetAspect(ROW / float(COL));
  pangolin::View& d_video = pangolin::Display("imgVideo").SetAspect(double(m_img_width) / m_img_height);
  pangolin::CreateDisplay()
	  .SetBounds(0.0, 0.3, pangolin::Attach::Pix(175), 1.0)
	  .SetLayout(pangolin::LayoutEqual)
	  .AddDisplay(d_video);

  pangolin::GlTexture coviVideo(600, 600, GL_RGB, false, 0, GL_RGB, GL_UNSIGNED_BYTE);
  //pangolin::View& d_video = pangolin::Display("imgVideo").SetAspect(ROW / float(COL));
  pangolin::View& covi_video = pangolin::Display("covisiblityVideo").SetAspect(double(600) / 600);
  pangolin::CreateDisplay()
	  .SetBounds(0.0, 0.3, 0.32, 0.5)
	  .SetLayout(pangolin::LayoutEqual)
	  .AddDisplay(covi_video);

  pangolin::OpenGlMatrix Twc;
  Twc.SetIdentity();

  bool bFollow = true;

  while (!m_bFinished)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(25));
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //get current pose in cv Mat
    getCurrentPose(Twc);

    //GetCurrentOpenGLCameraMatrix(Twc);

    if (menuFollowCamera && bFollow)
    {
      s_cam.Follow(Twc);
    }
    else if (menuFollowCamera && !bFollow)
    {
      s_cam.SetModelViewMatrix(pangolin::ModelViewLookAt(m_ViewpointX, m_ViewpointY, m_ViewpointZ, 0, 0, 0, 0.0, -1.0, 0.0));
      s_cam.Follow(Twc);
      bFollow = true;
    }
    else if (!menuFollowCamera && bFollow)
    {
      bFollow = false;
    }

    d_cam.Activate(s_cam);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    DrawCurrentCamera(Twc, Eigen::Vector3f(1.0f, 0.0f, 0.0f), m_CameraSize);

    if (menuShowKeyFrames)
    {
      DrawGroundTruth(menuShowKeyFrames);
      DrawKeyFrames(menuShowKF);
    DrawFilterWindowFrames(menuShowFilterWindow);
    }

    if (menuShowCG)
    DrawCovisibility();

    if (menuShowLF)
      DrawLocalizeFrames(menuShowLF);

    if (menuShowPoints)
      DrawMapPoints();
    if (menuShowTriangulatedPoints)
    DrawMSCKFTriangulatedPoints();
    if (menuShowPlanes)
      DrawPlanes();
    if (menuShowTracks)
    DrawTracks();
    if (menuShowAR)
      DrawARResult();
    if (menuShowInfo)
    {
      if (!ar_image.empty())
	      DrawInformation();
    }
    if(menuShowTracking && !bPauseRequested)
    {
      covi_video.Show(false);
      proj_video.Show(false);
      mapPointsView.Show(false);
      plotter_statescurve1.Show(true);
      plotter_statescurve2.Show(true);
      // Draw filter states
      plotter_statescurve1.Activate();
      mStatesLog1.Log(biasa_.x(),
	      biasa_.y(),
	      biasa_.z(),
	      radial_distortion_.x(),
	      radial_distortion_.y(),
	      tangential_distortion_.x(),
	      tangential_distortion_.y(),
	      time_delay_);

      plotter_statescurve2.Activate();
      double dfx = focal_point_calib_.x() - focal_point_.x();
      double dfy = focal_point_calib_.y() - focal_point_.y();
      double dcx = optical_center_calib_.x() - optical_center_.x();
      double dcy = optical_center_calib_.y() - optical_center_.y();
      mStatesLog2.Log(dfx, dfy, dcx, dcy);
    }
    else
    {
      plotter_statescurve1.Show(false);
      plotter_statescurve2.Show(false);
      covi_video.Show(true);
      proj_video.Show(true);
      mapPointsView.Show(true);
      if(!covisibility_map.empty())
	      coviVideo.Upload(covisibility_map.data, GL_BGR, GL_UNSIGNED_BYTE);
      covi_video.Activate();
      glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
      coviVideo.RenderToViewportFlipY();
      if (!projection_viewer.empty())
	      projVideo.Upload(projection_viewer.data, GL_BGR, GL_UNSIGNED_BYTE);
      proj_video.Activate();
      glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
      projVideo.RenderToViewportFlipY();
      if (!matMapPointsSize.empty())
	      mapPointsVideo.Upload(matMapPointsSize.data, GL_BGR, GL_UNSIGNED_BYTE);
      mapPointsView.Activate();
      glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
      mapPointsVideo.RenderToViewportFlipY();
    }

    //draw image
    if (!ar_image.empty())
      texVideo.Upload(ar_image.data, GL_BGR, GL_UNSIGNED_BYTE);
    d_video.Activate();
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    texVideo.RenderToViewportFlipY();

    pangolin::FinishFrame();

    if (menuClose)
    {
      SetFinish();
    }
    if (menuReset)
    {
      unique_lock<mutex> lk(reset_mtx);
      resetSystem = true;
      menuReset = false;
    }

  }
  mbFinished = false;
}

bool ViewerPangolin::BufferDataEmpty()
{
  return false;
}
