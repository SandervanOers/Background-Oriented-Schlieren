#ifndef PIXEL_TRANSLATION
#define PIXEL_TRANSLATION

#include <opencv2/opencv.hpp>
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculatePixelTranslation(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset);
/*--------------------------------------------------------------------------*/
extern std::vector<std::vector<double>> calculatePixelTranslationRandom(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const std::vector<double> &X_positions, const std::vector<double> &Y_positions);
/*--------------------------------------------------------------------------*/
extern std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const double &Indexi, const double &Indexj);
/*--------------------------------------------------------------------------*/
#endif
