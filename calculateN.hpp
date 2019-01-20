#ifndef DetN
#define DetN

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <numeric>

#include "PositionDirection.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}
using namespace cv;
/*--------------------------------------------------------------------------*/
static double calculateMean(const cv::Mat &Mat);
/*--------------------------------------------------------------------------*/
extern double calculateLf(const double &focal_length, const double &Lm);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateDirectionCosines(const cv::Mat &GridX, const cv::Mat &GridY, const double &L_f, const double &Distance_From_Pixels_To_Meters);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateIntersectionPlaneLine(const std::vector<cv::Mat> &InitialPosition, const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLaw(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLawMat(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const cv::Mat &n_0, const cv::Mat &n_1);
/*--------------------------------------------------------------------------*/
extern PositionDirection calculateIntersectionRungeKutta(const std::vector<cv::Mat> &InitialPosition, const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const cv::Mat &n_field, const double &n_0, const double &n_2, const unsigned int &SplineDegree, const double &L_t, const unsigned int &Number_Of_Steps);
/*--------------------------------------------------------------------------*/
extern PositionDirection calculateIntersectionConstantRefraction(const PositionDirection &InitialPositionDirection, const std::vector<double> &PlaneDefinition, const double &n_initial, const double &n_final);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> ForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &focal_length, const double &Lm, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const cv::Mat &n_field, const unsigned int &SplineDegree, const unsigned int &Number_Of_Steps);
/*--------------------------------------------------------------------------*/
static double getDerivativeValue(float *fptr_img1, const unsigned int &cols, const unsigned int &rows, const double &x, const double &y, const unsigned int &SplineDegree, const unsigned int &direction);
/*--------------------------------------------------------------------------*/
cv::Mat calculateTransformationMatrix(const std::vector<double> &PlaneDefinition);
/*--------------------------------------------------------------------------*/
#endif