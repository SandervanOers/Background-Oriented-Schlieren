
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
extern std::vector<cv::Point> get_valid_Neighbours(const cv::Mat &M_valid_points, const cv::Mat &Computed_Points, const unsigned int &x, const unsigned int &y, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned &offset);
/*--------------------------------------------------------------------------*/
static std::vector<double> iteration_rigid_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
static std::vector<double> iteration_affine_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
static std::vector<double> iteration_irregular_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
static std::vector<double> iteration_quadratic_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
static double Correlation_Coefficient_ZNSSD(const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &sum_g_minus_gm_squared, const std::vector<double> &f_values, const std::vector<double> &g_values);
/*--------------------------------------------------------------------------*/
static void calculate_Hessian_Jacobian_quadratic(cv::Mat &Hessian, cv::Mat &Jacobian, const cv::Mat &img, float *fptr_img1, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &Indexi, const unsigned int &Indexj, const std::vector<double> &f_values, const double &fm, const double &sum_f_minus_fm_squared, const std::vector<double> &g_values, const double &gm, const double &sum_g_minus_gm_squared, const double &lambda);
/*--------------------------------------------------------------------------*/
//static double getDerivativeValue(float *fptr_img1, const unsigned int &cols, const unsigned int &rows, const double &x, const double &y, const unsigned int &SplineDegree, const unsigned int &direction);
/*--------------------------------------------------------------------------*/
#endif
