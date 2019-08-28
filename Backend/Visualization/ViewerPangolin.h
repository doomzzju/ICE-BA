#ifndef _VIEWER_PANGOLIN_H_
#define _VIEWER_PANGOLIN_H_

#include <opencv2/opencv.hpp>
#include <pangolin/pangolin.h>
#include <mutex>
#include <vector>
#include <memory>

#include "IBA.h"
#include "LocalMap.h"
#include "GlobalMap.h"

class ViewerPangolin : public MT::Thread {
public:
  virtual void Initialize(IBA::Solver *solver, const int img_width = 640, const int img_height = 480);
  virtual void Reset();
  virtual void Run();
  
protected:
  virtual bool BufferDataEmpty();
  
  void DrawCurrentCamera();
  
private:
  float m_ViewpointX, m_ViewpointY, m_ViewpointZ, m_ViewpointF;
  
  bool m_bFinished;
  
  //draw map params
  float m_KeyFrameSize;
  float m_KeyFrameLineWidth;
  float m_PointSize;
  float m_CameraSize;
  float m_CameraLineWidth;
  float m_boxSize;
  
  int m_img_width, m_img_height;
};

#endif // _VIEWER_PANGOLIN_H_