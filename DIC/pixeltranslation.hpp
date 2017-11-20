#ifndef PIXEL_TRANSLATION
#define PIXEL_TRANSLATION

#include <opencv2/opencv.hpp>
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculatePixelTranslation(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset);
/*--------------------------------------------------------------------------*/
#endif
