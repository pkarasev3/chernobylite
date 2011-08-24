#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>


#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

namespace chernobylite
{
using std::string;
using std::vector;

template<typename Type>
void fillWeightsGaussian(cv::Mat& weights, Type sigma_squared)
{
  for (int y = 0; y < weights.rows; y++)
  {
    for (int x = 0; x < weights.cols; x++)
    {
      Type y_h = ((Type)y) / (weights.rows - 1.0) - 0.5;
      Type x_h = ((Type)x) / (weights.cols - 1.0) - 0.5;
      x_h *= 2; //x goes from -1 to 1
      y_h *= 2; //y "" ""
      Type val = ((x_h * x_h) + (y_h*y_h))/(2*sigma_squared);
      val = (Type)std::exp(-val);
      weights.at<Type> (y, x) = val;
    }
  }
}

void fillWeightsGaussian32(cv::Mat& weights, float sigma_squared);
void fillWeightsGaussian64(cv::Mat& weights, double sigma_squared);

inline bool readKfromCalib(cv::Mat& K, cv::Mat& distortion, cv::Size & img_size, const std::string& calibfile)
{
  cv::FileStorage fs(calibfile, cv::FileStorage::READ);
  cv::Mat cameramat;
  cv::Mat cameradistortion;
  float width = 0, height = 0;
  if (fs.isOpened())
  {
    cv::read(fs["camera_matrix"], cameramat, cv::Mat());
    cv::read(fs["distortion_coefficients"], cameradistortion, cv::Mat());
    cv::read(fs["image_width"], width, 0);
    cv::read(fs["image_height"], height, 0);

    fs.release();

  }
  else
  {
    throw std::runtime_error("bad calibration!");
  }

  cv::Size _size(width, height);
  img_size = _size;

  cameramat.convertTo(K,CV_32F);
  distortion = cameradistortion;
  return true;
}

/** draw axes on top of image to show pose */
void poseDrawer(cv::Mat& drawImage, const cv::Mat& K,
                       const cv::Mat& w, const cv::Mat& t, 
                       const std::string scaleText = std::string(""), int lineThickness=4);


/**  given a "dir" as string and ending extension, put name of files
     into the vector string. vector is sorted lexicographically.
  */
void lsFilesOfType(const char * dir, const string& extension,
        vector<string>& files) {
    files.clear();
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(dir)) == NULL) {
        return;
    }

    while ((dirp = readdir(dp)) != NULL) {
        std::string name(dirp->d_name);
        size_t pos = name.find(extension);
        if (pos != std::string::npos) {
            files.push_back(name);
        }
    }
    closedir(dp);
    std::sort(files.begin(), files.end());
}

unsigned char*  mat2raw( const cv::Mat&  frame );

void raw2mat(unsigned char*  VideoData, cv::Mat& frame);

double pointMeanInOutDiffCost( const cv::Mat& img, const cv::Point2f& pt_xy, int rad );



void drawHistogram_internal(cv::Mat& img, const cv::Mat& query,
                            int npts, int lineHeight, const std::string& name, cv::Point2f xy_draw_start);

template<typename T>
void drawHistogramOfVectorOnImage( const std::vector<T>& vec, const std::vector<cv::Point2f>& xy,
                                   cv::Mat& img, cv::Point2f xy0, int npts = 64,
                                   int lineHeight=48, const std::string& name = std::string("") )
{
  assert( ! img.empty() ); // must be created already
  cv::Mat query(vec);
  drawHistogram_internal(img,query,npts,lineHeight,name,xy0);
}

void drawHistogramOfVectorOnImageF32( const std::vector<float>& vec, const std::vector<cv::Point2f>& xy,
                                   cv::Mat& img, cv::Point2f xy0, int npts = 64, int lineHeight=48,
                                   const std::string& name = std::string("") );

void computeHomography_NotNormalized( const cv::Mat& img0, const cv::Mat& img1,
                                      cv::Mat& H, cv::Mat& drawImg, double fscale=1.0);

void computeHomography_NotNormalized( const std::vector<cv::Point2f>& srcPts,
                                      const std::vector<cv::Point2f>& dstPts,
                                      cv::Mat& H, cv::Mat& drawImg);

}
