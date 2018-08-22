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
/*--------------------------------------------------------------------------*/
extern std::vector<std::vector<double>>  calculatePixelTranslationRandom(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const std::vector<double> &X_positions, const std::vector<double> &Y_positions)
{
    std::vector<double> DispX, DispY, CorrelationCoefficient;
    for (unsigned int k = 0; k<X_positions.size(); k++)
    {
        unsigned int Indexi = X_positions[k];
        unsigned int Indexj = Y_positions[k];

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
        DispY.push_back(tempi);
        tempi = (double)matchLoc.x+ SubsetLength/2 - (double)Indexi;
        DispX.push_back(tempi);
        CorrelationCoefficient.push_back(maxVal);

        std::cout << "Percentage Done: " << (double)k/(X_positions.size())*100 << "%" << std::endl;
    }
    std::vector<std::vector<double>>  ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(CorrelationCoefficient);
    return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const double &Indexi, const double &Indexj)
{
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
        double DispX, DispY, CorrelationCoefficient;
        DispY = (double)matchLoc.y + (SubsetLength/2) - (double)Indexj;
        DispX = (double)matchLoc.x+ SubsetLength/2 - (double)Indexi;
        CorrelationCoefficient = maxVal;

    std::vector<double> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(CorrelationCoefficient);
    return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads)
{
	std::cout << "Computing Initial Guess" << std::endl;
	cv::Mat DispX(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat DispY(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat CorrelationCoefficient(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);

	std::vector<std::thread> threads;
	for (unsigned int l = 0; l < Number_Of_Threads; l++)
	{
		unsigned int xl = 1+round(0.1*horx_ROI/GridLength);
		unsigned int xr = horx_ROI/GridLength-round(0.1*horx_ROI/GridLength);
		unsigned int yl = 1+round(0.1*very_ROI/GridLength) + (double)l/Number_Of_Threads * (very_ROI/GridLength-round(0.1*very_ROI/GridLength)-(1+round(0.1*very_ROI/GridLength)));
		unsigned int yr = 1+round(0.1*very_ROI/GridLength) + (l+1.0)/Number_Of_Threads * (very_ROI/GridLength-round(0.1*very_ROI/GridLength)-(1+round(0.1*very_ROI/GridLength)))-1;

		threads.push_back(std::thread(calculateInitialGuess_Thread, l, Number_Of_Threads, Reference, Deformed, std::ref(DispX), std::ref(DispY), std::ref(CorrelationCoefficient), SubsetLength, GridLength, offset, xl, xr, yl, yr));
	}  
	
	for (auto& th : threads) 
		th.join();
    
	std::cout << "Computation Initial Guess Completed " << std::endl;
	std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(CorrelationCoefficient);
    return ReturnVector;
	
}
/*--------------------------------------------------------------------------*/
void calculateInitialGuess_Thread(const unsigned int &tid, const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &CorrelationCoefficient, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr)
{
	thread_local std::uniform_int_distribution<unsigned> ux(xl, xr);
    thread_local std::uniform_int_distribution<unsigned> uy(yl, yr);
	std::random_device rdx{};
    std::random_device rdy{};
    std::default_random_engine ex{rdx()};
    std::default_random_engine ey{rdy()};
    auto dicex = std::bind(ux, ex);
    auto dicey = std::bind(uy, ey);	
	
	//std::cout << "Total Points: " << (xr-xl+1)*(yr-yl+1) << std::endl;
    for (unsigned int k = 0; k < 100/Number_Of_Threads; k++) // 1 promille of total points // std::max(100., 0.01*(xr-xl)*(yr-yl))
    {
    // Pick X_positions and Y_positions randomly on Grid (from GridLength)
    // Get random i in range (1,horx_ROI/GridLength) => unsigned int Indexi = offset+SubsetLength/2 + i * GridLength; => X_positions(k) = Indexi;
    // Get random j in range (1,very_ROI/GridLength) => unsigned int Indexj = offset+SubsetLength/2 + j * GridLength; => Y_posiitons(k) = Indexj;
    // Note: Removed boundaries from range => Can give weird results
    // Do calculatePixelTranslation for these randomly chosen points
		unsigned int i = dicex();
        unsigned int j = dicey();
        unsigned int Indexi = offset+SubsetLength/2 + i * GridLength;
        unsigned int Indexj = offset+SubsetLength/2 + j * GridLength;
        std::vector <double> Displ;
        Displ = calculatePixelTranslationRandom_SinglePoint(Reference, Deformed, SubsetLength, Indexi, Indexj);

        DispX.at<double>(j,i) = Displ[0];
        DispY.at<double>(j,i) = Displ[1];
        CorrelationCoefficient.at<double>(j,i) = Displ[2];
	}
}
/*--------------------------------------------------------------------------*/