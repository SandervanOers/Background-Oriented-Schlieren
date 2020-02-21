#ifndef CSVImage
#define CSVImage
#include <opencv2/opencv.hpp>
/*--------------------------------------------------------------------------*/
extern cv::Mat load_matrix(std::string path, std::string filename, const int &skiplines);
/*--------------------------------------------------------------------------*/
#endif
