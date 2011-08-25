#include <chernobylite/opencv_helpers.hpp>
#include <chernobylite/math_util.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <fstream>
#include <iomanip>

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/tracking.hpp>

#include <ctype.h>
#include <vector>
#include <iostream>
#include <sstream>

using namespace cv;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using boost::lexical_cast;

namespace {

/** weights is an MxN matrix. Fill it in with Gaussian shape.
  * sigma_squared: smaller => more of the gaussian is near the center.
  * suggested: between 0.01-1.0 for sigma_squared
  */
template<typename T>
void fillWeightsGaussian(cv::Mat& weights, T sigma_squared)
{
  for (int y = 0; y < weights.rows; y++)
  {
    for (int x = 0; x < weights.cols; x++)
    {
      T y_h = ((T)y) / (weights.rows - 1.0) - 0.5;
      T x_h = ((T)x) / (weights.cols - 1.0) - 0.5;
      x_h *= 2; //x goes from -1 to 1
      y_h *= 2; //y "" ""
      T val = ((x_h * x_h) + (y_h*y_h))/(2*sigma_squared);
      val = (T)std::exp(-val);
      weights.at<T> (y, x) = val;
    }
  }
}

}

namespace chernobylite  {



/** does a *MALLOC*, you own the pointer it returns!
  */
unsigned char*  mat2raw( const Mat&  frame )
{
  //*ProcessFrameY8( short Width, short Height, unsigned char *VideoData )
  bool    isRGB   = (frame.channels() == 3);
  long int nBytes = frame.cols * frame.rows * (1+2*isRGB);
  unsigned char*  VideoData = (unsigned char*) malloc( nBytes );

  std::memcpy( VideoData, frame.ptr<unsigned char>(0), nBytes );
  return VideoData;
}

/** does a *COPY*, safe but inefficient!
  */
void raw2mat( unsigned char*  VideoData, Mat& frame)
{
  int     nchans  = frame.channels();
  bool    isRGB   = (nchans == 3);
  long int nBytes = frame.cols * frame.rows * (1+2*isRGB);

  assert( nBytes > 16 );
  std::memcpy( frame.ptr<unsigned char>(0), VideoData, nBytes );

}

double pointMeanInOutDiffCost( const cv::Mat& img, const cv::Point2f& pt_xy, int rad )
{
  cv::Point3f RGBi(0,0,0);
  cv::Point3f RGBo(0,0,0);
  int num_i = 1;
  int num_o = 1;

  int imin = std::max(0.0f,pt_xy.y - rad );
  int jmin = std::max(0.0f,pt_xy.x - rad );
  int imax = std::min(img.rows * 1.0f,pt_xy.y + rad );
  int jmax = std::min(img.cols * 1.0f,pt_xy.x + rad );
  for( int i = imin; i <= imax; i++ ) {
    for( int j = jmin; j <= jmax; j++ ) {
      double rad_ij = std::sqrt( std::pow( double(i - pt_xy.y),2.0) + std::pow( double(j - pt_xy.x),2.0) );
      cv::Point3f rgb_ij = img.at<cv::Point3f>(i,j);
      if( rad_ij < std::max(2.0,rad * 0.75) ) {
        RGBi += rgb_ij; // interior mean rgb
        num_i++;
      } else if( rad_ij >= rad ) {
        RGBo += rgb_ij; // exterior mean rgb
        num_o++;
      }
    }
  }
  double E = normL2(RGBi*(1.0/num_i) - RGBo*(1.0/num_o));
  return E;
}

void fillWeightsGaussian32(cv::Mat& weights, float sigma_squared)
{
  fillWeightsGaussian<float>(weights,sigma_squared);
}

void fillWeightsGaussian64(cv::Mat& weights, double sigma_squared)
{
  fillWeightsGaussian<double>(weights,sigma_squared);
}

namespace {

std::map<std::string,cv::Point2f>  range_of_named_plot;
std::map<std::string,cv::Mat>      cache_of_named_plot;
typedef std::map<std::string,cv::Point2f> name2xy;
typedef std::map<std::string,cv::Mat> name2mat;
cv::Mat hist;
}

void drawHistogram_internal(cv::Mat& img, const cv::Mat& query,
                            int npts, int lineHeight, const std::string& name,
                            cv::Point2f xy_draw_start )
{

  float ir[2];
  double dmin, dmax;
  int N = query.cols * query.rows;
  npts = (npts > N ? N : npts);

  cv::minMaxLoc(query,&dmin,&dmax);
  dmax = (dmax > dmin + 1e-3 ) ? dmax : dmin + 1e-3;

  name2xy::iterator iter = range_of_named_plot.find(name) ;
  if( range_of_named_plot.end() == iter ) {
    Mat new_hist = Mat::zeros(1,npts,CV_32F);
    cache_of_named_plot.insert(name2mat::value_type(name,new_hist));
    range_of_named_plot.insert(name2xy::value_type(name,cv::Point2f(dmin,dmax)));
  } else {
    dmin = std::min(dmin, (double) range_of_named_plot[name].x );
    dmax = std::max(dmax, (double) range_of_named_plot[name].y );
    range_of_named_plot[name] = cv::Point2f(dmin,dmax);
  }
  //cache_of_named_plot[name].copyTo(hist);

  ir[0] = dmin;
  ir[1] = dmax;
  const float* irange = ir;
  cv::Mat mask = cv::Mat::ones(query.size(),CV_8U);

  cv::calcHist(&query,1,0,mask,cache_of_named_plot[name],1,&npts,&irange,true,true);
  double hrange[2];
  cv::minMaxLoc(cache_of_named_plot[name],&(hrange[0]),&(hrange[1]));

  cv::Scalar color(200,220,255,100);
  float fscale = 0.3;
  int step     = 3;
  for( int k = 0; k < npts; k++ ) {
    float yval      = -lineHeight * cache_of_named_plot[name].at<float>(k) / (1e-9+hrange[1]);
    cv::Point2f pt1 = xy_draw_start + cv::Point2f(step*k,0);
    cv::Point2f pt2 = xy_draw_start + cv::Point2f(step*k,yval);
    cv::line(img,pt1,pt2,color,1,CV_AA);
  }

  std::stringstream ss1,ss2;
  ss1 << std::setprecision(2) << std::setiosflags(std::ios_base::scientific) << dmin;
  cv::Point2f minlabel = xy_draw_start+cv::Point2f(-12,8);
  cv::Point2f maxlabel = xy_draw_start+cv::Point2f( npts*step*0.9,8);
  cv::putText(img, ss1.str(), minlabel, CV_FONT_HERSHEY_SIMPLEX, fscale, color, 1, CV_AA, false);
  ss2 << std::setprecision(2) << dmax;
  cv::putText(img, ss2.str(), maxlabel, CV_FONT_HERSHEY_SIMPLEX, fscale, color, 1, CV_AA, false);

  if( !name.empty() ) {
    int baseline = 0;
    cv::Size sz   = cv::getTextSize(name.c_str(), CV_FONT_HERSHEY_SIMPLEX, 1, fscale, &baseline);
    cv::Point2f coordToPutName = (maxlabel+minlabel)*0.5+cv::Point2f(0,0);
    cv::putText(img, name, coordToPutName, CV_FONT_HERSHEY_SIMPLEX, fscale, color, 1, CV_AA, false);
  }
}


void drawHistogramOfVectorOnImageF32( const std::vector<float>& vec, const std::vector<cv::Point2f>& xy,
                                      cv::Mat& img, cv::Point2f xy0, int npts, int lineHeight,
                                      const std::string& name )
{
  drawHistogramOfVectorOnImage(vec,xy,img,xy0,npts,lineHeight,name);
}




namespace
{
// @ descriptor_extractor_matcher.cpp in opencv samples
// enforce symmetry of matching 1 to 2, 2 to 1
void crossCheckMatching( Ptr<DescriptorMatcher>& descriptorMatcher,
                         const Mat& descriptors1, const Mat& descriptors2,
                         vector<DMatch>& result_matches12, int knn=1 )
{
  result_matches12.clear();
  vector<vector<DMatch> > matches12, matches21;

  descriptorMatcher->knnMatch(descriptors1, descriptors2, matches12, knn );
  descriptorMatcher->knnMatch( descriptors2, descriptors1, matches21, knn );

  for( size_t m = 0; m < matches12.size(); m++ )
  {
    bool findCrossCheck = false;
    for( size_t fk = 0; fk < matches12[m].size(); fk++ )
    {
      DMatch forward = matches12[m][fk];

      for( size_t bk = 0; bk < matches21[forward.trainIdx].size(); bk++ )
      {
        DMatch backward = matches21[forward.trainIdx][bk];
        if( backward.trainIdx == forward.queryIdx )
        {
          result_matches12.push_back(forward);
          findCrossCheck = true;
          break;
        }
      }
      if( findCrossCheck ) break;
    }
  }
}
}

void poseDrawer(cv::Mat& drawImage, const cv::Mat& K,
                       const cv::Mat& w, const cv::Mat& t,
                       const std::string scaleText, int lineThickness)
{
  using namespace cv;
  Point3f z(0, 0, -0.25);
  Point3f x(0.25, 0, 0);
  Point3f y(0, -0.25, 0);
  Point3f o(0, 0, 0);
  vector<Point3f> op(4);
  op[1] = x, op[2] = y, op[3] = z, op[0] = o;
  vector<Point2f> ip;

  Mat D = Mat::zeros(4,1,CV_32F);
  projectPoints(Mat(op), w, t, K, D, ip);
  double axes_sz = drawImage.rows / 4.0;
  double zmin    = 5e-2;

  ip[1] = ip[0] + (ip[1]-ip[0] ) * ( axes_sz / norm( ip[1] - ip[0] ) );
  ip[2] = ip[0] + (ip[2]-ip[0] ) * ( axes_sz / norm( ip[2] - ip[0] ) );
  ip[3] = ip[0] + (ip[3]-ip[0] ) * ( (1.0/sqrt(2))*axes_sz / ( zmin + norm( ip[3] - ip[0] ) ) );

  // DRAW AXES LINES
  vector<Scalar> c(4); //colors
  c[0] = Scalar(255, 255, 255);
  c[1] = Scalar(205, 50, 50);//x
  c[2] = Scalar(100, 200, 0);//y
  c[3] = Scalar(200, 100, 205);//z
  line(drawImage, ip[0], ip[1], c[1],lineThickness,CV_AA);
  line(drawImage, ip[0], ip[2], c[2],lineThickness,CV_AA);
  line(drawImage, ip[0], ip[3], c[3],lineThickness,CV_AA);

  if( scaleText.size() > 1 )
  { // print some text on the image if desired
    int baseline = 0;
    Size sz = getTextSize(scaleText, CV_FONT_HERSHEY_SIMPLEX, 1, 1, &baseline);
    rectangle(drawImage, Point(10, 30 + 5),
              Point(10, 30) + Point(sz.width, -sz.height - 5), Scalar::all(0), -1);
    putText(drawImage, scaleText, Point(10, 30), CV_FONT_HERSHEY_SIMPLEX, 1.0, c[0], 1, CV_AA, false);
  }

  // DRAW LETTERS FOR AXES
  c[1] += Scalar(50,50,50);
  c[2] += Scalar(50,50,50);
  c[3] += Scalar(50,50,50);
  putText(drawImage, "Z", ip[3], CV_FONT_HERSHEY_SIMPLEX, 1.0, c[3], lineThickness, CV_AA, false);
  putText(drawImage, "Y", ip[2], CV_FONT_HERSHEY_SIMPLEX, 1.0, c[2], lineThickness, CV_AA, false);
  putText(drawImage, "X", ip[1], CV_FONT_HERSHEY_SIMPLEX, 1.0, c[1], lineThickness, CV_AA, false);

}

bool readKfromCalib(cv::Mat& K, cv::Mat& distortion, cv::Size & img_size, const std::string& calibfile)
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

void computeHomography_NotNormalized(const cv::Mat &img0, const cv::Mat &img1,
                                     cv::Mat &H, cv::Mat& drawImg,std::vector<double>& WT_result,
                                      double fscale )
{

  Ptr<FeatureDetector>     detector  = FeatureDetector::create( "ORB" /*"DynamicSURF"*/ );
  Ptr<DescriptorExtractor> extractor = DescriptorExtractor::create( "BRIEF"/*"SURF"*/);
  Ptr<DescriptorMatcher>   matcher   = DescriptorMatcher::create("BruteForce-Hamming" /*"BruteForce"*/);

  assert( detector != NULL && extractor != NULL && matcher != NULL );

  // First Image
  cout << endl << "< Extracting keypoints from first image..." << endl;
  vector<KeyPoint> keypoints0;
  detector->detect( img0, keypoints0 );
  cout << keypoints0.size() << " points" << endl << ">" << endl;
  cout << "< Computing descriptors for keypoints from first image..." << endl;
  Mat descriptors0;
  extractor->compute( img0, keypoints0, descriptors0 );
  cout << ">" << endl;

  // Second Image
  cout << endl << "< Extracting keypoints from second image..." << endl;
  vector<KeyPoint> keypoints1;
  detector->detect( img1, keypoints1 );
  cout << keypoints1.size() << " points" << endl << ">" << endl;
  cout << "< Computing descriptors for keypoints from second image..." << endl;
  Mat descriptors1;
  extractor->compute( img1, keypoints1, descriptors1 );
  cout << ">" << endl;

  vector<DMatch> result_matches;
  int knn = 1;
  cout << "type0: " << descriptors0.type()
       << ", type1: " << descriptors1.type() << endl;
  crossCheckMatching( matcher, descriptors0, descriptors1, result_matches, knn );

  vector<int> queryIdxs( result_matches.size() );
  vector<int> trainIdxs( result_matches.size() );
  for( size_t i = 0; i < result_matches.size(); i++ )
  {
    queryIdxs[i] = result_matches[i].queryIdx;
    trainIdxs[i] = result_matches[i].trainIdx;
  }

  cout << "< Computing homography (RANSAC)..." << endl;
  vector<Point2f> points0; KeyPoint::convert(keypoints0, points0, queryIdxs);
  vector<Point2f> points1; KeyPoint::convert(keypoints1, points1, trainIdxs);

  vector<Point3f> points0xyz(points0.size());

  double imgW = img0.cols;
  double imgH = img0.rows;
  double f0   = (fscale > 1.0 ) ? fscale : sqrt(imgW*imgH);

  for( size_t i1 = 0; i1 < points1.size(); i1++ )
  {
    points0xyz[i1].x = (points0[i1].x - (imgW-1)/2.0) / f0;
    points0xyz[i1].y = (points0[i1].y - (imgH-1)/2.0) / f0;
    points0xyz[i1].z = 0.0;
  }


  // struct AdaptiveRansacParams { ... };
  double reprojThresh = 1.0/fscale;
  int num_matched = 0;
  int min_N       = 8;
  int max_N       = min_N + 5;
  int maxIters    = 64;
  int iters       = 0;
  double stepVerh = 1.5;
  double stepVniz = 0.95;
  min_N           = std::max( min_N, (int) (points0.size() * 0.01) );
  max_N           = std::max( max_N, (int) (points0.size() * 0.10) );
  vector<char> matchesMask( result_matches.size(), 0 );
  Mat H01;


  vector<Point3f> points0xyz_keep; // best fit in the homography solving
  vector<Point2f> points1_keep;    // best fit in the homography solving
  points0xyz_keep.reserve(max_N);
  points1_keep.reserve(max_N);

  while( ( (num_matched < min_N) || (num_matched > max_N) ) && (iters < maxIters) )
  {// relax the threshold until we have a good set
    iters++;
    points0xyz_keep.clear();
    points1_keep.clear();
    num_matched = 0;

    H01 = findHomography( Mat(points0), Mat(points1), CV_RANSAC, reprojThresh );
    cout << ">" << endl;

    Mat points0_mapped;
    perspectiveTransform(Mat(points0), points0_mapped, H01);

    for( size_t i1 = 0; i1 < points1.size(); i1++ )
    { // record the ones that reproject within threshold
      double reproj_error = normL2(points1[i1] - points0_mapped.at<Point2f>((int)i1,0));
      if( reproj_error <= reprojThresh ) { // inlier
        matchesMask[i1] = 1;
        num_matched++;
        points0xyz_keep.push_back(points0xyz[i1]);
        points1_keep.push_back(points1[i1]);
      } else { // fail this index
        matchesMask[i1] = 0;
      }
    }
    cout << "iters" << iters << ", matched " << num_matched << " of "
         << points1.size() << ", thresh = " << reprojThresh << endl;
    reprojThresh *= (num_matched < min_N ) ? stepVerh : stepVniz;
    if( 0 == (iters % 32 ) ) {
      stepVerh = (1.0 + stepVerh)/2.0;
      stepVniz = (1.0 + stepVniz)/2.0;
    }
  }
  H01.convertTo(H,CV_32F);
  drawImg = cv::Mat();

  // draw inliers
  cv::drawMatches( img0, keypoints0, img1, keypoints1, result_matches, drawImg,
                   CV_RGB(0, 255, 0), CV_RGB(0, 0, 255), matchesMask );
  Mat tvec,rvec;
  Mat K = Mat::eye(3,3,CV_32F);
  K.at<float>(0,0) = f0;
  K.at<float>(1,1) = f0;
  K.at<float>(0,2) = (imgW-1.0)/2.0;
  K.at<float>(1,2) = (imgH-1.0)/2.0;
  if( 1 ) {
    cv::solvePnP(Mat(points0xyz_keep),Mat(points1_keep),K,Mat()/*distortion*/,
                       rvec,tvec,false/*init-guess*/);
    cout << "rvec = " << rvec << ", tvec = " << tvec << endl;
    vector<double> WT(6);
    assert( rvec.type() == CV_64F && tvec.type() == CV_64F );
    for( int k=0; k < 3; k++ ) {
      WT[k]   = rvec.at<double>(k,0);
      WT[k+3] = tvec.at<double>(k,0);
    }
    WT_result = WT;

  }

}



void write_RT_to_csv( const string& csvFile, const vector<double>& WT ) {
  cout << "attempting to write csv file in " << csvFile << endl;
  std::ofstream  data(csvFile.c_str());
  stringstream ss;
  // ss << "wt_out= " << Mat(WT) << ";"; matlab-format (?)
  ss <<WT[0]<< ","<< WT[1]<< ","<<WT[2]<<","<<WT[3]<<","<<WT[4]<<","<<WT[5];
  data << ss.str() << endl;
  data.close();
}

void load_RT_from_csv( const string& csvFile, vector<float>& WT ) {
  cout << "attempting to open csv file in " << csvFile << endl;
  std::ifstream  data(csvFile.c_str());
  std::string line;
  while(std::getline(data,line))
  {
    std::stringstream  lineStream(line);
    std::string        cell;
    int column_index = 0;
    while(std::getline(lineStream,cell,',')) {
      if(column_index > 5) { cerr << "Oh snap, bogus CSV file!" << endl; exit(1); }
      WT[column_index++] = lexical_cast<float>(cell);
    }
  }
  data.close();
}



}

