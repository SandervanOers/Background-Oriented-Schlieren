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
extern std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const double &Indexi, const double &Indexj, const unsigned int &MaxPixelYVertical)
{
	cv::Mat temp = Reference(cv::Range(Indexj-SubsetLength/2,Indexj+SubsetLength/2+1), cv::Range(Indexi-SubsetLength/2, Indexi+SubsetLength/2+1));
	
	double DispX=0, DispY=0, CorrelationCoefficient=-1;
	
	// check wheter the reference image is (nearly) uniform. 
	cv::Mat     mean;
	cv::Mat     stddev;
	cv::meanStdDev (temp, mean, stddev );
	double       mean_pxl = mean.data[0];
	double       stddev_pxl = stddev.data[0];
	if (stddev_pxl>=1)
	{	
		unsigned int dyl, dyr;
		if (Indexj>MaxPixelYVertical)
		{
			dyl = Indexj-MaxPixelYVertical;
		}
		else
		{
			dyl = 0;
		}		
		if (Indexj+MaxPixelYVertical+1 < Deformed.rows)
		{
			dyr = Indexj+MaxPixelYVertical+1;
		}
		else
		{
			dyr = Deformed.rows;
		}
		cv::Mat Deformed2 = Deformed(cv::Range(dyl,dyr), cv::Range::all());
		
		cv::Mat result;
		double minVal; double maxVal;
		cv::Point minLoc; cv::Point maxLoc;
		cv::Point matchLoc;
		cv::matchTemplate(Deformed2, temp, result, CV_TM_CCOEFF_NORMED);
		// Localizing the best match with minMaxLoc
		cv::minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );
		matchLoc = maxLoc;
		DispY = (double)matchLoc.y + (SubsetLength/2) + dyl-(double)Indexj;//- MaxPixelYVertical;// - (double)Indexj;
		DispX = (double)matchLoc.x + SubsetLength/2 - (double)Indexi;
		CorrelationCoefficient = maxVal;
	}
    std::vector<double> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(CorrelationCoefficient);
    return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads, const unsigned int &MaxPixelYVertical)
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

		threads.push_back(std::thread(calculateInitialGuess_Thread, l, Number_Of_Threads, Reference, Deformed, std::ref(DispX), std::ref(DispY), std::ref(CorrelationCoefficient), SubsetLength, GridLength, offset, xl, xr, yl, yr, MaxPixelYVertical));
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
void calculateInitialGuess_Thread(const unsigned int &tid, const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &CorrelationCoefficient, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr, const unsigned int &MaxPixelYVertical)
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
    for (unsigned int k = 0; k < ceil(100/Number_Of_Threads); k++) // 1 promille of total points // std::max(100., 0.01*(xr-xl)*(yr-yl))
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
        Displ = calculatePixelTranslationRandom_SinglePoint(Reference, Deformed, SubsetLength, Indexi, Indexj, MaxPixelYVertical);
        //Displ = calculatePixelTranslationRandom_SinglePoint(Reference(cv::Range(Indexj-SubsetLength/2,Indexj+SubsetLength/2+1), cv::Range(Indexi-SubsetLength/2, Indexi+SubsetLength/2+1)), Deformed, SubsetLength, Indexi, Indexj);

        DispX.at<double>(j,i) = Displ[0];
        DispY.at<double>(j,i) = Displ[1];
        CorrelationCoefficient.at<double>(j,i) = Displ[2];
	}
}
/*--------------------------------------------------------------------------*/
void calculateInitialGuess_Thread_Iteration(const unsigned int &tid, const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction)
{
	thread_local std::uniform_int_distribution<unsigned> ux(xl, xr);
    thread_local std::uniform_int_distribution<unsigned> uy(yl, yr);
	std::random_device rdx{};
    std::random_device rdy{};
    std::default_random_engine ex{rdx()};
    std::default_random_engine ey{rdy()};
    auto dicex = std::bind(ux, ex);
    auto dicey = std::bind(uy, ey);
	
	// 0.8 is the minimum correlation coefficient we find acceptable
	double DX=0, DY=0, CC=0.8, I=0, J=0;
    for (unsigned int k = 0; k < ceil(100/Number_Of_Threads); k++) 
    {
    // Pick X_positions and Y_positions randomly on Grid (from GridLength)
    // Get random i in range (1,horx_ROI/GridLength) => unsigned int Indexi = offset+SubsetLength/2 + i * GridLength; => X_positions(k) = Indexi;
    // Get random j in range (1,very_ROI/GridLength) => unsigned int Indexj = offset+SubsetLength/2 + j * GridLength; => Y_posiitons(k) = Indexj;
    // Do calculatePixelTranslation for these randomly chosen points
		unsigned int i = dicex();
        unsigned int j = dicey();
        unsigned int Indexi = offset+SubsetLength/2 + i * GridLength;
        unsigned int Indexj = offset+SubsetLength/2 + j * GridLength;
        std::vector <double> Displ;
        Displ = calculatePixelTranslationRandom_SinglePoint(Reference, Deformed, SubsetLength, Indexi, Indexj, MaxPixelYVertical);
        //Displ = calculatePixelTranslationRandom_SinglePoint(Reference(cv::Range(Indexj-SubsetLength/2,Indexj+SubsetLength/2+1), cv::Range(Indexi-SubsetLength/2, Indexi+SubsetLength/2+1)), Deformed, SubsetLength, Indexi, Indexj);
		
		if (Displ[2] > CC)
		{
			DX = Displ[0];
			DY = Displ[1];
			CC = Displ[2];
			I = i;
			J = j;
		}
	}
	// Use the Initial Guess to Calculate the Iterated Solution 
	std::vector<double> InitialCondition = {DX, DY, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<double> point1 = iteration(Reference, fptr_img1, I, J, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
	
	if (point1.back()>CC)
	{
		DispX.at<double>(J,I) = point1[0];
		DispY.at<double>(J,I)= point1[1];
		Ux.at<double>(J,I) = point1[2];
		Vx.at<double>(J,I) = point1[3];
		Uy.at<double>(J,I) = point1[4];
		Vy.at<double>(J,I) = point1[5];
		Uxy.at<double>(J,I) = point1[6];
		Vxy.at<double>(J,I) = point1[7];
		Uxx.at<double>(J,I) = point1[8];
		Vxx.at<double>(J,I)= point1[9];
		Uyy.at<double>(J,I) = point1[10];
		Vyy.at<double>(J,I) = point1[11];
		CorrelationCoefficient.at<double>(J,I) = point1[12];
		Computed_Points.at<uchar>(J,I) = 1;
	}
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess_Iteration(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction)
{
	std::cout << "Computing Initial Guess" << std::endl;
	cv::Mat DispX(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat DispY(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat CorrelationCoefficient(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat Computed_Points(DispX.size(), CV_8UC1, 0.0);
    cv::Mat Ux(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vx(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uxy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vxy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uxx(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vxx(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uyy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vyy(DispX.size(), CV_64FC1, 0.0);
	
	std::vector<std::thread> threads;
	for (unsigned int l = 0; l < Number_Of_Threads; l++)
	{
		unsigned int xl = 1+round(0.1*horx_ROI/GridLength);
		unsigned int xr = horx_ROI/GridLength-round(0.1*horx_ROI/GridLength);
		unsigned int yl = 1+round(0.1*very_ROI/GridLength) + (double)l/Number_Of_Threads * (very_ROI/GridLength-round(0.1*very_ROI/GridLength)-(1+round(0.1*very_ROI/GridLength)));
		unsigned int yr = 1+round(0.1*very_ROI/GridLength) + (l+1.0)/Number_Of_Threads * (very_ROI/GridLength-round(0.1*very_ROI/GridLength)-(1+round(0.1*very_ROI/GridLength)))-1;
		// std::cout << l << " " << xl << " " << xr << " " << yl << " " << yr << std::endl; 
		// Possible Error when Number_Of_Threads is large and GridLength is large: yl == yr or even yl>yr.
		threads.push_back(std::thread(calculateInitialGuess_Thread_Iteration, l, Number_Of_Threads, Reference, Deformed, fptr_img1, std::ref(DispX), std::ref(DispY), std::ref(Ux), std::ref(Vx), std::ref(Uy), std::ref(Vy), std::ref(Uxy), std::ref(Vxy), std::ref(Uxx), std::ref(Vxx), std::ref(Uyy), std::ref(Vyy), std::ref(CorrelationCoefficient), std::ref(Computed_Points), SplineDegree, SubsetLength, GridLength, offset, xl, xr, yl, yr, MaxPixelYVertical, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction));
	
	}  
	for (auto& th : threads) 
		th.join();
	
	std::cout << "Computation Initial Guess Completed " << std::endl;
	std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(Ux);
    ReturnVector.push_back(Vx);
    ReturnVector.push_back(Uy);
    ReturnVector.push_back(Vy);
    ReturnVector.push_back(Uxy);
    ReturnVector.push_back(Vxy);
    ReturnVector.push_back(Uxx);
    ReturnVector.push_back(Vxx);
    ReturnVector.push_back(Uyy);
    ReturnVector.push_back(Vyy);
    ReturnVector.push_back(CorrelationCoefficient);
    ReturnVector.push_back(Computed_Points);
    return ReturnVector;
}