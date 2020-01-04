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
static cv::Mat load_matrix(std::string path, std::string filename, const int &skiplines);
/*--------------------------------------------------------------------------*/
static bool sort_by_C_value (const Points_With_Value &lhs, const Points_With_Value &rhs);
/*--------------------------------------------------------------------------*/
static void compute_Save_GridX_Y(const cv::Size &Size, const unsigned int &xStart_ROI, const unsigned int &yStart_ROI, const unsigned int &GridLength, const std::string path);
/*--------------------------------------------------------------------------*/
class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, ','))
            {
                m_data.push_back(cell);
            }
            // This checks for a trailing comma with no data after it.
            if (!lineStream && cell.empty())
            {
                // If there was a trailing comma then add an empty element.
                m_data.push_back("");
            }
        }
    private:
        std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   
/*--------------------------------------------------------------------------*/
class CSVIterator
{   
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};
/*--------------------------------------------------------------------------*/
void removeCharsFromString( std::string &str, char* charsToRemove ) {
   for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
      str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   }
}
/*--------------------------------------------------------------------------*/
/*std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str)
{
    std::vector<std::string>   result;
    std::string                line;
    std::getline(str,line);

    std::stringstream          lineStream(line);
    std::string                cell;

    while(std::getline(lineStream,cell, ','))
    {
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}*/
/*--------------------------------------------------------------------------*/
#endif
