#ifndef H_DIC
#define H_DIC

#include "pixeltranslation.hpp"
#include "nonlineariteration.hpp"
#include "PointsWithValue.hpp"
#include "InputVariables.hpp"
#include "InputOut.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}

/*--------------------------------------------------------------------------*/
extern void DIC(const cv::Mat &img, const cv::Mat &img1, const InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
static void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored);
/*--------------------------------------------------------------------------*/
static bool sort_by_C_value (const Points_With_Value &lhs, const Points_With_Value &rhs);
/*--------------------------------------------------------------------------*/
static void compute_Save_GridX_Y(const cv::Size &Size, const InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
#endif