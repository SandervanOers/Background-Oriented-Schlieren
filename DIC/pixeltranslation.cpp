# include "pixeltranslation.hpp"
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculatePixelTranslation(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset)
    {
    unsigned int horx_ROI = Reference.cols-SubsetLength/2-SubsetLength/2-2*offset;
    unsigned int very_ROI = Reference.rows-SubsetLength/2-SubsetLength/2-2*offset;
    cv::Mat DispX((very_ROI)/GridLength+1, (horx_ROI)/GridLength+1, CV_64F);
    cv::Mat DispY((very_ROI)/GridLength+1, (horx_ROI)/GridLength+1, CV_64F);
    cv::Mat CorrelationCoefficient((very_ROI)/GridLength+1, (horx_ROI)/GridLength+1, CV_64F);

    for (unsigned int i = 0; i < DispX.cols; i++)
    {
        for (unsigned int j = 0; j < DispX.rows; j++)
        {
            unsigned int Indexi = offset+SubsetLength/2 + i * GridLength;
            unsigned int Indexj = offset+SubsetLength/2 + j * GridLength;
            cv::Mat temp = Reference(cv::Range(Indexj-SubsetLength/2,Indexj+SubsetLength/2+1), cv::Range(Indexi-SubsetLength/2, Indexi+SubsetLength/2+1));
            cv::Mat result;
            double minVal; double maxVal;
            cv::Point minLoc; cv::Point maxLoc;
            cv::Point matchLoc;
            // Compute match
            cv::matchTemplate(Deformed, temp, result, CV_TM_CCOEFF_NORMED);
            // Localizing the best match with minMaxLoc
            cv::minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );
            matchLoc = maxLoc;
            double tempi;
            tempi = (double)matchLoc.y + (SubsetLength/2) - (double)Indexj;
            DispY.at<double>(j,i) = tempi;
            tempi = (double)matchLoc.x+ SubsetLength/2 - (double)Indexi;
            DispX.at<double>(j,i) = tempi;
            CorrelationCoefficient.at<double>(j,i) = maxVal;
        }
        std::cout << "Percentage Done: " << (double)i/(DispX.cols-1)*100 << "%" << std::endl;
    }
    std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(CorrelationCoefficient);
    return ReturnVector;
    }
