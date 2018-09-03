
#ifndef ITERATION
#define ITERATION

#include <opencv2/opencv.hpp>
#include <math.h>
#include <thread>
#include <numeric>
#include "PointsWithValue.hpp"
extern "C" {
#include "interpol.h"
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> iteration(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction);
/*--------------------------------------------------------------------------*/
void iteration_thread(std::vector<std::vector<double>> &SolutionNeighbours, const unsigned int k, const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction);
/*--------------------------------------------------------------------------*/
void iteration_Thread(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, std::vector<Points_With_Value> &Locations_Best_Correlation);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Point> get_valid_Neighbours(const cv::Mat &M_valid_points, const cv::Mat &Computed_Points, const unsigned int &x, const unsigned int &y, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned &offset);
/*--------------------------------------------------------------------------*/
#endif
