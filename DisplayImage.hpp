#ifndef DImage
#define DImage

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <cppad/cppad.hpp>

#include <math.h>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <random>
#include <functional>
#include <sys/types.h>
#include <sys/stat.h>
#include <thread>
#include <future>
//#include <iostream>
//#include <fstream>

#include "InputOut.hpp"
#include "pixeltranslation.hpp"
#include "nonlineariteration.hpp"
#include "PointsWithValue.hpp"
#include "PositionDirection.hpp"
#include "calculateN.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}
using namespace cv;
/*--------------------------------------------------------------------------*/
static void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored);
/*--------------------------------------------------------------------------*/
static bool sort_by_C_value (const Points_With_Value &lhs, const Points_With_Value &rhs);
/*--------------------------------------------------------------------------*/
static void compute_Save_GridX_Y(const cv::Size &Size, const unsigned int &xStart_ROI, const unsigned int &yStart_ROI, const unsigned int &GridLength, const std::string path);
/*--------------------------------------------------------------------------*/
#endif