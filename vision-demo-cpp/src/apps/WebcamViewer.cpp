#include <sstream>
#include <fstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread/mutex.hpp>
//#include <boost/thread/thread.hpp>

using namespace cv;
using namespace std;

namespace po = boost::program_options;
using boost::lexical_cast;

namespace {

struct Options {
  string videoN;        // /dev/videoN
  int    frameWidth;
  int    frameHeight;
};

int options(int ac, char ** av, Options& opts) {
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "Produce help message.")
      ("video,V", po::value<string>(&opts.videoN)->default_value("0"),
                      "video source integer or string. use 300 for 1394")
      ("width,w", po::value<int>(&opts.frameWidth)->default_value(320),"frame width")
      ("height,h", po::value<int>(&opts.frameHeight)->default_value(240),"frame height");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }

  return 0;

}

struct MouseEvent
{
  MouseEvent( )  {
    click_point.x = 0;
    click_point.y = 0;
    is_clicked = false;
  }

  Point2f click_point;
  bool is_clicked;

  void drawClick( Mat& outimg ) {
    if( is_clicked ) { // draw the click point, if it exists
      cv::circle(outimg,click_point,5,Scalar(255,50,100),1,CV_AA);
      cv::circle(outimg,click_point,4,Scalar(255,255,255),1,CV_AA);
    }
  }

  void operator()(int event, int x, int y, int flags)  {
    if (flags & CV_EVENT_FLAG_LBUTTON )    {
      is_clicked = true;
      click_point = Point2f(x, y);
      cout << "clicked point: " << Mat(click_point) << endl;
    }
  }
};

static void onMouse(int event, int x, int y, int flags, void* userdata)
{
  MouseEvent* data = (MouseEvent*)userdata;
  (*data)(event, x, y, flags);
}

Mat frame, outimg;
vector<Mat> image_container(1);
bool isClickedToStart = false;
string window_name    = "anonymous";

} // end app-local namespace


int main(int argc, char** argv) {
  Options opts;
  if (options(argc, argv, opts)) // failed to parse!
    return 1;

//  boost::shared_ptr<KatyushaTracker> prob( new KatyushaTracker );
//  PFTracker  pft( prob );

  int vid_int = boost::lexical_cast<int>(opts.videoN);
  VideoCapture capture( vid_int ); // at /dev/videoN
  if( !capture.isOpened() ) {
    cout << "opening as integer video code " << vid_int
         << "failed, trying as string" << endl;
    capture.open(opts.videoN);
  }
  if (!capture.isOpened()) {
    cerr << "unable to open video device " << opts.videoN << endl;
    return 1;
  }
  capture.set(CV_CAP_PROP_FRAME_WIDTH,opts.frameWidth);
  capture.set(CV_CAP_PROP_FRAME_HEIGHT,opts.frameHeight);

  namedWindow(window_name);
  MouseEvent mousetrap;
  setMouseCallback(window_name, onMouse, &mousetrap);
  vector<Mat> channels(3);
  Point3f rgbPoint;
  int saveCount = 0;
  for (;;)   // Main Loop
  {
    capture >> frame; // grab frame data from webcam
    if (frame.empty())
      continue;
    if (outimg.empty())
      frame.copyTo(outimg);

    { // Display
      outimg = frame * 0.1 + outimg * 0.9;
      mousetrap.drawClick(outimg);
      imshow( window_name, outimg );
    }

    if( mousetrap.is_clicked && !isClickedToStart )
    { // initialize tracker on the click point
      isClickedToStart = true;
    }

    if( isClickedToStart ) {
      int px = mousetrap.click_point.x;
      int py = mousetrap.click_point.y;
      split(outimg,channels); assert(frame.type()==CV_8UC3);
      rgbPoint.x = channels[0].at<u_int8_t>(py,px);
      rgbPoint.y = channels[1].at<u_int8_t>(py,px);
      rgbPoint.z = channels[2].at<u_int8_t>(py,px);
      cout << "RGB point:  " << Mat(rgbPoint) << "\r";
      cout.flush();
    }

    char key = waitKey(20);
    if( 'q' == key ) // quit if we hit q key
      break;

    if( 's' == key ) {
      stringstream ss;
      ss << argv[0] << "_"
         << std::setw(4) << std::setfill('0') << saveCount << ".png";
      imwrite(ss.str(),outimg); cout << "saved " << ss.str() << endl;
      saveCount++;
    }
  }
  return 0;



}
