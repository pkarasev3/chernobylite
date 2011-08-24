
/** pkarasev3  18 july 2011
  *   test accuracy / validity of getting rigid pose from features
  *   detected on template image and an observation, along with the
  *   homography H matrix.
  */

#include "chernobylite/opencv_helpers.hpp"
#include <boost/lexical_cast.hpp>

using namespace cv;
using namespace std;
using namespace chernobylite;

int main( int ac, char* av [] )
{
  if( ac < 3 ) {
    cerr << "usage: " << av[0] << " img1.jpg  img2.jpg [fscale](optional)" << endl;
    exit(1);
  }
  string name0 = av[1];
  string name1 = av[2];
  Mat img0rgb = imread( name0 );
  Mat img1rgb = imread( name1 );
  assert( !img0rgb.empty() && !img1rgb.empty() );

  double f0 = img0rgb.cols * 1.00; // default approximation for focal-scale
  if( ac > 3 ) {
    f0 = boost::lexical_cast<double>(av[3]);
  }
  cout << "using f-scale " << f0 << endl;
  Mat H, drawImg;
  vector<double> WT;
  computeHomography_NotNormalized(img0rgb,img1rgb,H,drawImg, WT, f0);
  string winName = "correspondences. q to exit.";

  string outFile = "wt_feature_points.out";
  write_RT_to_csv(outFile,WT);
  cout << Mat(WT) << " \n ...  wrote to " << outFile << endl;

  char key = '@';
  while( key != 'q' ) {
    imshow( winName, drawImg );
    key = waitKey(100);
  }

  // 1) try the PnP
  // 2) "look for modified f, enforce origin-to-origin,
  //                   orthonormal R columns" after first-pass

  return 0;
}
