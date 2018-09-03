#ifndef PIXEL_TRANSLATION
#define PIXEL_TRANSLATION

#include <opencv2/opencv.hpp>
#include <random>
#include <thread>
#include "nonlineariteration.hpp"
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculatePixelTranslation(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset);
/*--------------------------------------------------------------------------*/
extern std::vector<std::vector<double>> calculatePixelTranslationRandom(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const std::vector<double> &X_positions, const std::vector<double> &Y_positions);
/*--------------------------------------------------------------------------*/
extern std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const double &Indexi, const double &Indexj);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads, const unsigned int &MaxPixelYVertical);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess_Iteration(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction);
/*--------------------------------------------------------------------------*/
void calculateInitialGuess_Thread(const unsigned int &tid, const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &CorrelationCoefficient, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr, const unsigned int &MaxPixelYVertical);
/*--------------------------------------------------------------------------*/
void calculateInitialGuess_Thread_Iteration(const unsigned int &tid, const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction);
/*--------------------------------------------------------------------------*/
#endif
