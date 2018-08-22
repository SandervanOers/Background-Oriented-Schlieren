#ifndef PIXEL_TRANSLATION
#define PIXEL_TRANSLATION

#include <opencv2/opencv.hpp>
#include <random>
#include <thread>
#include <mutex>
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculatePixelTranslation(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset);
/*--------------------------------------------------------------------------*/
extern std::vector<std::vector<double>> calculatePixelTranslationRandom(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const std::vector<double> &X_positions, const std::vector<double> &Y_positions);
/*--------------------------------------------------------------------------*/
extern std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const double &Indexi, const double &Indexj);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads);
/*--------------------------------------------------------------------------*/
void calculateInitialGuess_Thread(const unsigned int &tid, const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &CorrelationCoefficietnt, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr);
/*--------------------------------------------------------------------------*/
#endif
