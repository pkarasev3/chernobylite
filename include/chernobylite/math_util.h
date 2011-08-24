#include <vector>
#include <cmath>
#include <vector>
#include <opencv2/core/core.hpp>

/** fast norm functions. profiling code indicates that cv::norm
  * wastes a lot of CPU, very bad for deeply nested loops!
  */
inline double normL1(const cv::Point3f& p)
{
  return abs(p.x) + abs(p.y) + abs(p.z);
}
inline double normINF(const cv::Point3f& p)
{
  return std::max(std::max(abs(p.x), abs(p.y)), abs(p.z));
}

inline double normL2(const cv::Point3f& p)
{
  return std::sqrt(p.dot(p));
}
inline double normL1(const cv::Point2f& p)
{
  return abs(p.x) + abs(p.y);
}
inline double normINF(const cv::Point2f& p)
{
  return std::max(abs(p.x), abs(p.y));
}
inline double normL2(const cv::Point2f& p)
{
  return std::sqrt(p.dot(p));
}
