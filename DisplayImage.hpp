#ifndef DImage
#define DImage

#include <stdio.h>
#include <opencv2/opencv.hpp>

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

#endif