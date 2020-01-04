#include "DisplayImage.hpp"

int main(int argc, char** argv )
{
	auto tr1 = std::chrono::high_resolution_clock::now();

	/*--------------------------------------------------------------------------*/
	// Check and Process Input
	std::string pathname_string;
	unsigned int SplineDegree, SubsetLength, GridLength, ShapeFunction, propagation_function, ordering, xStart, xEnd, yStart, yEnd, offset, Number_Of_Threads, MaxPixelYVertical;
	double abs_tolerance_threshold, rel_tolerance_threshold, minimum_corrcoeff_IG;
	int returnvalue = checkinput(argc, argv, pathname_string, SplineDegree, SubsetLength, GridLength, ShapeFunction, propagation_function, ordering, xStart, xEnd, yStart, yEnd, offset, Number_Of_Threads, MaxPixelYVertical, abs_tolerance_threshold, rel_tolerance_threshold, minimum_corrcoeff_IG);
	if (returnvalue < 0)
	{
		return -1;
	}
	std::cout << std::endl << "\033[1;32mInput Parsed\033[0m\n" << std::endl;
	cv::String path(pathname_string);
	/*--------------------------------------------------------------------------*/
	// Read Images from Disk
	cv::Mat img, img1;
	unsigned int xStart_ROI, yStart_ROI, horx_ROI, very_ROI;
	returnvalue = readImageDataFromFile(path, img, img1, xStart, xEnd, yStart, yEnd, SubsetLength, offset, xStart_ROI, yStart_ROI, horx_ROI, very_ROI, ordering);
	if (returnvalue < 0)
	{
		return -1;
	}
	std::cout << "\033[1;32mImages Loaded\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
    //cv::Mat M_valid_points = get_valid_points(img, SubsetLength, offset);
    cv::Mat M_valid_points(img.size(), CV_8U, Scalar(0));
    for (unsigned int i = 0; i < static_cast<unsigned int>(img.rows); i++)
    {
        for (unsigned int j = 0; j < static_cast<unsigned int>(img.cols); j++)
        {
            if (i >= SubsetLength/2+offset && i <= static_cast<unsigned int>(img.rows)-SubsetLength/2-offset && j >= SubsetLength/2+offset && j <= static_cast<unsigned int>(img.cols)-SubsetLength/2-offset)
            {
                M_valid_points.at<uchar>(i,j) = 1;
            }
        }
    }
	/*--------------------------------------------------------------------------*/
    // Calculate interpolation coefficients for g
    cv::Mat prova_img1= img1.clone();
    float *fptr_img1 = prova_img1.ptr<float>(0);
    SamplesToCoefficients(fptr_img1, img1.cols, img1.rows, SplineDegree);
	std::cout << "\033[1;32mB-splines coefficients calculated\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
	// TO DO: Need to rewrite nonlineariteration.cpp to use img1.cols and img1.rows instead of img.cols and img.rows since interpol.c needs width and height of image.
	std::vector<cv::Mat> IG = calculateInitialGuess_Iteration(img, img1, fptr_img1, SplineDegree, SubsetLength, GridLength, horx_ROI, very_ROI, offset, Number_Of_Threads, MaxPixelYVertical, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction, minimum_corrcoeff_IG);
	cv::Mat DispX = IG[0].clone();
	cv::Mat DispY = IG[1].clone();
	cv::Mat Ux =  	IG[2].clone();
	cv::Mat Vx =   	IG[3].clone();
	cv::Mat Uy =    IG[4].clone();
	cv::Mat Vy =    IG[5].clone();
	cv::Mat Uxy =   IG[6].clone();
	cv::Mat Vxy =   IG[7].clone();
	cv::Mat Uxx =   IG[8].clone();
	cv::Mat Vxx =   IG[9].clone();
	cv::Mat Uyy =   IG[10].clone();
	cv::Mat Vyy =   IG[11].clone();
	cv::Mat CorrelationCoefficient =  IG[12].clone();
	cv::Mat Computed_Points =  IG[13].clone();

	compute_Save_GridX_Y(DispX.size(), xStart_ROI, yStart_ROI, GridLength, pathname_string);
	std::cout << "\033[1;32mInitial Points Computed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
    std::vector<Points_With_Value> Locations_Best_Correlation;
    cv::Mat nonZeroCoordinates;
    cv::findNonZero(CorrelationCoefficient>0, nonZeroCoordinates);
    for (unsigned int i = 0; i < nonZeroCoordinates.total(); i++ )
	{
		Locations_Best_Correlation.push_back(Points_With_Value(CorrelationCoefficient.at<double>(nonZeroCoordinates.at<Point>(i)), nonZeroCoordinates.at<Point>(i)));
    }

    // Sort first elements
    sort(Locations_Best_Correlation.begin(), Locations_Best_Correlation.end(), sort_by_C_value);
    for (auto i = Locations_Best_Correlation.begin(); i < Locations_Best_Correlation.end(); i++)
	{
		std::cout << (*i).Loc << ": " << (*i).C_value << std::endl;
	}
	std::cout << std::endl << "\033[1;32mList of Points with Values Sorted\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
    const bool plotting = 0;
    int ct = 0;
    std::string name = "CC_";
    std::string type = ".png";
    std::string folderName = "Output/";
    if (plotting)
    {
        std::stringstream ss;
        ss<<name<<(ct)<<type;
        std::string filename = ss.str();
        ss.str("");
        cv::Mat Copy_CCF = CorrelationCoefficient.clone();
        Copy_CCF.convertTo(Copy_CCF, CV_8UC1, 255.0);
        imwrite(filename, Copy_CCF);
    }
	//auto tr5 = std::chrono::high_resolution_clock::now();

	//std::chrono::duration<double> TR = tr5-tr5;
	//auto T = TR.count();
	/*--------------------------------------------------------------------------*/

	cv::Ptr<cv::Formatter> fmt = cv::Formatter::get(cv::Formatter::FMT_DEFAULT);
	fmt->set64fPrecision(2);
	fmt->set32fPrecision(2);

	int outputcount = 1;
	while(static_cast<unsigned int>(cv::sum(Computed_Points).val[0]) < Computed_Points.total())
    {
        // Get best point in queue
        cv::Point matchLoc = Locations_Best_Correlation[0].Loc;
        // Use Initial Guess from (neighbouring) initial point
        std::vector<double> InitialCondition = {DispX.at<double>(matchLoc), DispY.at<double>(matchLoc), Ux.at<double>(matchLoc), Vx.at<double>(matchLoc), Uy.at<double>(matchLoc), Vy.at<double>(matchLoc), Uxy.at<double>(matchLoc), Vxy.at<double>(matchLoc), Uxx.at<double>(matchLoc), Vxx.at<double>(matchLoc), Uyy.at<double>(matchLoc), Vyy.at<double>(matchLoc)};
        std::vector<cv::Point> Neighbours = get_valid_Neighbours(M_valid_points, Computed_Points, matchLoc.x, matchLoc.y, SubsetLength, GridLength, offset);
        for (auto i = Neighbours.begin(); i < Neighbours.end(); i++)
        {
            if (propagation_function==1 && GridLength < SubsetLength/2)
            {
                //std::cout << "InitialCondition old: " << InitialCondition[0] << ", " << InitialCondition[1] << std::endl;
                InitialCondition[0] += Ux.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength + Uy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                InitialCondition[1] += Vx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength + Vy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                //std::cout << "InitialCondition new: " << InitialCondition[0] << ", " << InitialCondition[1] << std::endl;
                ///InitialCondition[2] += Uxx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;// + Uxy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                //InitialCondition[3] += Vxx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;// + Vxy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                //InitialCondition[4] += Uyy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;// + Uxy.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;
                ///InitialCondition[5] += Vyy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;// + Vxy.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;
            }
			//auto tr3 = std::chrono::high_resolution_clock::now();
			//std::cout << "before" << std::endl;
            std::vector<double> point2 = iteration(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
			//std::cout << "after" << std::endl;
			//auto tr4 = std::chrono::high_resolution_clock::now();
			//std::chrono::duration<double> elapsed_seconds = tr4-tr3;
			//T += elapsed_seconds.count();
            DispX.at<double>(*i) = point2[0];
            DispY.at<double>(*i) = point2[1];
            Ux.at<double>(*i) = point2[2];
            Vx.at<double>(*i) = point2[3];
            Uy.at<double>(*i) = point2[4];
            Vy.at<double>(*i) = point2[5];
            Uxy.at<double>(*i) = point2[6];
            Vxy.at<double>(*i) = point2[7];
            Uxx.at<double>(*i) = point2[8];
            Vxx.at<double>(*i) = point2[9];
            Uyy.at<double>(*i) = point2[10];
            Vyy.at<double>(*i) = point2[11];
			//std::cout << "11" << std::endl;
            CorrelationCoefficient.at<double>(*i) = point2.back();
			//std::cout << "12: " << point2.back() << std::endl;
            Computed_Points.at<uchar>(*i) = 1;
			//std::cout << "13" << std::endl;

			//std::cout << fmt->format(CorrelationCoefficient) << std::endl;
			//std::cout << *i << std::endl << std::endl;
			//std::cout << std::fixed << std::setprecision(2) << CorrelationCoefficient << std::endl;
            Locations_Best_Correlation.push_back(Points_With_Value(point2.back(), *i));

        }
		//std::cout << "Before Loc"<< std::endl;
        Locations_Best_Correlation.erase(Locations_Best_Correlation.begin());
        sort(Locations_Best_Correlation.begin(), Locations_Best_Correlation.end(), sort_by_C_value);

    //for (auto i = Locations_Best_Correlation.begin(); i < Locations_Best_Correlation.end(); i++)
	//{
	//	std::cout << (*i).Loc << ": " << (*i).C_value << std::endl;
	//}
	//std::cout << "After Loc"<< std::endl;

        if (Neighbours.empty())
        {
            //std::cout << "no valid or uncomputed neighbours" << std::endl;
        }
        else
        {
            if (plotting)
            {
                ct++;
                std::stringstream ss;
                ss<<name<<(ct)<<type;
                std::string filename = ss.str();
                ss.str("");
                cv::Mat Copy_CCF = CorrelationCoefficient.clone();
                Copy_CCF.convertTo(Copy_CCF, CV_8UC1, 255.0);
                imwrite(filename, Copy_CCF);
            }
            else
            {
                double computed = cv::sum(Computed_Points).val[0]/static_cast<double>(Computed_Points.total())*100.0;
                if (computed>outputcount)
                {
                    outputcount++;
                    outputcount++;
                    outputcount++;
                    std::cout << "Points Computed: " << std::setprecision(2) << computed << "%" << std::endl;
                }
            }
        }
    }
	std::cout << std::endl << "\033[1;32mAll Points Computed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
	/*
	        // Get best point in queue
			auto tr3 = std::chrono::high_resolution_clock::now();
			std::vector<std::future<std::vector<double>>> futures;

        cv::Point matchLoc = Locations_Best_Correlation[0].Loc;
        // Use Initial Guess from (neighbouring) initial point
        std::vector<double> InitialCondition = {DispX.at<double>(matchLoc), DispY.at<double>(matchLoc), Ux.at<double>(matchLoc), Vx.at<double>(matchLoc), Uy.at<double>(matchLoc), Vy.at<double>(matchLoc), Uxy.at<double>(matchLoc), Vxy.at<double>(matchLoc), Uxx.at<double>(matchLoc), Vxx.at<double>(matchLoc), Uyy.at<double>(matchLoc), Vyy.at<double>(matchLoc)};
        std::vector<cv::Point> Neighbours = get_valid_Neighbours(M_valid_points, Computed_Points, matchLoc.x, matchLoc.y, SubsetLength, GridLength, offset);
		std::vector<std::vector<double>> List_Neighbour_Solution;

        for (auto i = Neighbours.begin(); i < Neighbours.end(); i++)
       // for (auto i = Neighbours_list.begin(); i < Neighbours_list.end(); i++)
        {
            if (propagation_function==1 && GridLength < SubsetLength/2)
            {
                InitialCondition[0] += Ux.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength + Uy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                InitialCondition[1] += Vx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength + Vy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
            }
			//futures.emplace_back(std::async(std::launch::async, iteration, img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction));
			futures.emplace_back(std::async(iteration, img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction));
			//std::vector<double> point2 = iteration(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
			//List_Neighbour_Solution.push_back(point2);
		}

			Locations_Best_Correlation.erase(Locations_Best_Correlation.begin());



		for (auto &future : futures) {
		   //List_Neighbour_Solution.push_back(future.get());
		//}

		//for (auto i = 0; i < Neighbours.size(); i++)
		//{
			//std::vector<double> point2 = List_Neighbour_Solution[i];
			std::vector<double> point2 = future.get();
			cv::Point Loc(point2[13], point2[14]);
			DispX.at<double>(Loc) = point2[0];
			DispY.at<double>(Loc) = point2[1];
			Ux.at<double>(Loc) = point2[2];
			Vx.at<double>(Loc) = point2[3];
			Uy.at<double>(Loc) = point2[4];
			Vy.at<double>(Loc) = point2[5];
			Uxy.at<double>(Loc) = point2[6];
			Vxy.at<double>(Loc) = point2[7];
			Uxx.at<double>(Loc) = point2[8];
			Vxx.at<double>(Loc) = point2[9];
			Uyy.at<double>(Loc) = point2[10];
			Vyy.at<double>(Loc) = point2[11];
			CorrelationCoefficient.at<double>(Loc) = point2[12];
			Computed_Points.at<uchar>(Loc) = 1;
			Locations_Best_Correlation.push_back(Points_With_Value(point2[12], Loc));
		}
		auto tr4 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = tr4-tr3;
		//std::cout << "Time iteration = " << elapsed_seconds.count() << std::endl;
		//tr5 = tr5 + elapsed_seconds.count();
		T += elapsed_seconds.count();
        sort(Locations_Best_Correlation.begin(), Locations_Best_Correlation.end(), sort_by_C_value);

        //if (Neighbours.empty())
        //{
        //    //std::cout << "no valid or uncomputed neighbours" << std::endl;
        //}
       // else
       // {
       //     if (plotting)
       //     {
        //        ct++;
        //        std::stringstream ss;
        //        ss<<name<<(ct)<<type;
        //        std::string filename = ss.str();
        //        ss.str("");
        //        cv::Mat Copy_CCF = CorrelationCoefficient.clone();
        //        Copy_CCF.convertTo(Copy_CCF, CV_8UC1, 255.0);
        //        imwrite(filename, Copy_CCF);
         //   }
        //    else
         //   {
                double computed = cv::sum(Computed_Points).val[0]/Computed_Points.total()*100;
                if (computed>outputcount)
                {
                    outputcount++;
                    outputcount++;
                    std::cout << "Points Computed: " << std::setprecision(2) << computed << "%" << std::endl;
                }
            //}
        //}
    }
clTabCtrl	*/
	std::cout << "Reliability = " << std::setprecision(5) << cv::sum(CorrelationCoefficient).val[0]/static_cast<double>(CorrelationCoefficient.total()) << std::endl;

    if (plotting==1)
    {
        std::stringstream ss;
        ss<<folderName<<"CC_PF"<<(propagation_function)<<"_Bn"<<(SplineDegree)<<"_SF"<<(ShapeFunction)<<"_GL"<<(GridLength)<<"_SL"<<(SubsetLength)<<type;
        std::string filename = ss.str();
        ss.str("");
        cv::Mat Copy_CCF_16 = CorrelationCoefficient.clone();
        Copy_CCF_16.convertTo(Copy_CCF_16, CV_16U, 255.0*256.0);
        imwrite(filename, Copy_CCF_16);
    }
	
	/*--------------------------------------------------------------------------*/
    store_matrix(path,"U0", DispX);
    store_matrix(path,"V0", DispY);
    store_matrix(path,"Ux", Ux);
    store_matrix(path,"Vx", Vx);
    store_matrix(path,"Uy", Uy);
    store_matrix(path,"Vy", Vy);
    store_matrix(path,"Uxy", Uxy);
    store_matrix(path,"Vxy", Vxy);
    store_matrix(path,"Uxx", Uxx);
    store_matrix(path,"Uyy", Uyy);
    store_matrix(path,"Vxx", Vxx);
    store_matrix(path,"Vyy", Vyy);
    store_matrix(path,"CorrelationCoefficient", CorrelationCoefficient);
	/*--------------------------------------------------------------------------*/
	std::string filename = "U0";
	cv::Mat In = load_matrix(path, filename, 1);
	std::cout << In << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << DispX - In << std::endl;
	std::cout << std::endl;
	/*--------------------------------------------------------------------------*/
	
	//cv::Mat Copy_CCF_16 = CorrelationCoefficient.clone();
	//Copy_CCF_16.convertTo(Copy_CCF_16, CV_16U, 255.0*256.0);
	//imwrite("Images/CC.png", Copy_CCF_16);
	//
	//flip(Copy_CCF_16, Copy_CCF_16, 0);
	//imwrite("Images/CCflipped.png", Copy_CCF_16);

	auto tr2= std::chrono::high_resolution_clock::now();
	std::cout << "DIC took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr2-tr1).count()
	<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr2-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr2-tr1).count() << " minutes\n"<<std::endl;
	// std::cout << "Iterations took "<< T << " seconds" << std::endl;

	/*--------------------------------------------------------------------------*/
	if (ordering==1)
	{
		DispX = - DispX;
		DispY = - DispY;
	}
	// Median Filter
	/*
	DispX.convertTo(DispX, CV_32F);
	DispY.convertTo(DispY, CV_32F);
	medianBlur ( DispX, DispX, 3 );
	medianBlur ( DispY, DispY, 3 );
	DispX.convertTo(DispX, CV_64F);
	DispY.convertTo(DispY, CV_64F);
	 * */
	// Gaussian Filter
	//GaussianBlur(DispX, DispX, Size(5, 5), 0);
	//GaussianBlur(DispY, DispY, Size(5, 5), 0);
	// Store Again?
	 
	/*--------------------------------------------------------------------------*/
	auto tr3= std::chrono::high_resolution_clock::now();
    cv::Mat GridX(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat GridY(DispX.size(), CV_64FC1, Scalar(0));
	for (unsigned int i = 0; i < static_cast<unsigned int>(GridX.cols); i++)
	{
		for (unsigned int j = 0; j < static_cast<unsigned int>(GridX.rows); j++)
		{
			GridX.at<double>(j,i) = xStart_ROI + i*GridLength;
			GridY.at<double>(j,i) = yStart_ROI + j*GridLength;
		}
	}
	// Camera
	double focal_length = 50.0/1000.0;
	double Distance_From_Pixels_To_Meters = 6.45e-6;
	double n_0 = 1.0003;
	double n_1 = 1.52;
	double n = 1.333;
	// Lengths Small Tank
	double L_c = 10000;
	double L_t = 0.191;//0.168;//0.168;//
	double L_g = 0.004; //0.01
	double L_s = 0.0;
	// Lengths Full Tank
	/*
	double L_c = 1.60;
	double L_g = 5.8/1000;
	double L_t = 13.8/100-2*L_g;
	double L_s = 0.0;//0.523;
	 * */
	std::vector<double> Lengths{L_c, L_g, L_t, L_s};
	double corr_cut_off = 0.98;
	//CalibrationFigures(GridX, GridY, DispX, DispY, CorrelationCoefficient, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, path);
	/*--------------------------------------------------------------------------*/
    /*
	 * cv::Mat GX(1, 21, CV_64FC1, Scalar(0));
    cv::Mat GY(1, 21, CV_64FC1, Scalar(0));
	for (unsigned int i = 0; i < static_cast<unsigned int>(GX.cols); i++)
	{
		for (unsigned int j = 0; j < static_cast<unsigned int>(GX.rows); j++)
		{
			GX.at<double>(j,i) = i*20;
			GY.at<double>(j,i) = j*20;
		}
	}
	//std::cout << GX << std::endl;
	
	double a, b, c;
	double L_m = 1.778;
	std::vector<double> PlaneDefinition{0,0,-1,L_m};
	Lengths[0] = L_m-2*L_g-L_t-L_s;
	//calculate_Displacements(GX, GY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n); 
	std::cout << "L_c = " << L_m-2*L_g-L_t-L_s << std::endl;*/
	/*--------------------------------------------------------------------------*/
    double a, b, c, L_m, meanGridX, meanGridY;
	/*
	//std::vector<double> CalibrationNumbers = Calibration(GridX, GridY, DispX, DispY, CorrelationCoefficient, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, path, corr_cut_off);
	std::vector<double> CalibrationNumbers = Calibration2(GridX, GridY, DispX, DispY, CorrelationCoefficient, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, path, corr_cut_off);
	//std::vector<double> CalibrationNumbers = CalibrationExtended(GridX, GridY, DispX, DispY, CorrelationCoefficient, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, path);
	a = CalibrationNumbers[0];
	b = CalibrationNumbers[1];
	c = CalibrationNumbers[2];
	L_m = CalibrationNumbers[3];
	meanGridX = CalibrationNumbers[4];
	meanGridY = CalibrationNumbers[5];
	//L_t = CalibrationNumbers[4];
	*/
	std::cout << std::endl << "\033[1;32mCalibration Completed\033[0m\n" << std::endl;   
	
	a = -0.12082091 ;
	b = 0.12082091 ;
	c =  -0.43362714 ;
	L_m =  2.1608463 ;
	meanGridX =  555.04703 ;
	meanGridY =  -765.25704;
	 
	double normPD = calculateNorm(a, b, c);
	double d = -c/normPD*L_m;
	//double pd[] = {a/normPD, b/normPD, c/normPD, d};
	std::vector<double> PlaneDefinition{0, 0, 0, 0};
	PlaneDefinition[0] = a/normPD; //{a/normPD, b/normPD, c/normPD, d};
	PlaneDefinition[1] = b/normPD;
	PlaneDefinition[2] = c/normPD;
	PlaneDefinition[3] = d;
	//PlaneDefinition.assign(pd, pd+4);
	L_c = c/(a+b+c)*L_m-L_s-2.0*L_g-L_t;
	Lengths[0] = L_c;
	Lengths[2] = L_t;
	//calculateNFigures(GridX, GridY, DispX, DispY, CorrelationCoefficient, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, path);

    
	auto tr4= std::chrono::high_resolution_clock::now();
	std::cout << "Calibration took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr4-tr3).count()
	<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr4-tr3).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr4-tr3).count() << " minutes"<<std::endl;
	std::cout << std::endl << "\033[1;32mN Calculation Completed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
	//cv::Mat nfield = CalculateN(GridX, GridY, DispX, DispY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, path);
	cv::Mat nfield = CalculateN2(GridX, GridY, meanGridX, meanGridY, DispX, DispY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, path);
	store_matrix(path,"nfield", nfield);
	auto tr5= std::chrono::high_resolution_clock::now();
	std::cout << "n Calculation took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr5-tr4).count()
	<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr5-tr4).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr5-tr4).count() << " minutes"<<std::endl;
	std::cout << std::endl << "\033[1;32mN Calculation Completed\033[0m\n" << std::endl;	

	/*--------------------------------------------------------------------------*/
	/*
	auto tr3= std::chrono::high_resolution_clock::now();

	unsigned int SizeGrid = 6;
    //unsigned int offsetSplineDegree = 2*((SplineDegree+1)/2) > 5 ? 2*((SplineDegree+1)/2) : 5;
	//std::cout << "offsetSplineDegree = " << offsetSplineDegree << std::endl;
	cv::Mat GX(cv::Size(SizeGrid,SizeGrid), CV_64FC1, Scalar(0));
	cv::Mat GY(cv::Size(SizeGrid,SizeGrid), CV_64FC1, Scalar(0));
	// Camera
	double focal_length = 50.0/1000.0;
	double Distance_From_Pixels_To_Meters = 6.45e-6;
	double n_0 = 1;
	double n_1 = 1.5;
	// Lengths
	double L_t = 0.168;
	double L_g = 0.01;
	double L_s = 0;
	double Lm = 2.0;
	// Plane Definition
	double a = 1.0;
	double b = 1.0;
	double c = 5.0;
	double d = -c*Lm;
	// Normalize normal vector
	double normPD = sqrt(a*a+b*b+c*c);
	std::vector<double> PlaneDefinition{a/normPD, b/normPD, c/normPD, 1};
	PlaneDefinition[3] = -c/normPD*Lm;
	double L_c = c/normPD/(a/normPD+b/normPD+c/normPD)*Lm-L_s-2*L_g-L_t;
	std::cout << "L_c = " << L_c << std::endl;
	std::cout << "L_total = " << L_c+2*L_g+L_t+L_s << std::endl;
	unsigned int Number_Of_Steps = 1000;
	std::vector<double> Lengths{L_c, L_g, L_t, L_s};

	cv::Mat Dx(GX.size(), CV_64FC1, Scalar(0));
	cv::Mat Dy(GX.size(), CV_64FC1, Scalar(0));
	cv::Mat n_field(GX.size(), CV_64FC1);
	for (unsigned int i = 0; i < static_cast<unsigned int>(GX.cols); i++)
	{
		for (unsigned int j = 0; j < static_cast<unsigned int>(GX.rows); j++)
		{
			GX.at<double>(i,j) = i*200.0;
			GY.at<double>(i,j) = j*200.0;
			Dx.at<double>(i,j) = 10+i;
			Dy.at<double>(i,j) = j;
			n_field.at<double>(i,j) = 1.333;//1.0+i/(double)SizeGrid/2.0;
		}
	}
	std::cout << GX << std::endl << std::endl;
	std::cout << GY << std::endl << std::endl;
	cv::Ptr<cv::Formatter> formatMat=Formatter::get(cv::Formatter::FMT_DEFAULT);
	formatMat->set64fPrecision(4);

	std::vector<cv::Mat> D6 = ForwardModel(GX, GY, Dx, Dy, focal_length, Lm, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n_field, SplineDegree, Number_Of_Steps);

	std::cout << "X6 = " << std::endl << D6[0] << std::endl<< std::endl;
	std::cout << "Z6 = " << std::endl << D6[2] << std::endl<< std::endl;

	auto tr4= std::chrono::high_resolution_clock::now();

	std::cout << "Total took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr4-tr3).count()
	<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr4-tr3).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr4-tr3).count() << " minutes"<<std::endl;
	 }

	 std::cout << "1" << std::endl;
	 int pp = poly_test();
	 std::cout << "2" << std::endl;
*/
	auto trend = std::chrono::high_resolution_clock::now();
	std::cout << "Total took " << std::chrono::duration_cast<std::chrono::milliseconds>(trend-tr1).count()
	<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(trend-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(trend-tr1).count() << " minutes"<<std::endl;
    return 0;
}
/*--------------------------------------------------------------------------*/
void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored)
{
    std::ofstream myfile;
    myfile.open(path+"/"+filename+".csv");
    myfile << format(Matrix_To_Be_Stored,cv::Formatter::FMT_MATLAB);
    myfile.close();
}
/*--------------------------------------------------------------------------*/
cv::Mat load_matrix(std::string path, std::string filename, const int &skiplines)
{

    std::ifstream file(path+"/"+filename+".csv");
    std::ifstream file1(path+"/"+filename+".csv");
	int line = 0;
	int numberofrows, numberofcols;
    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
		if (line>=skiplines)
		{
			numberofrows = (*loop).size();
		}
		line++;
    }
	numberofcols = line-skiplines;
	
	cv::Mat In(numberofcols, numberofrows, CV_64FC1);
	line = 0;
    for(CSVIterator loop(file1); loop != CSVIterator(); ++loop)
    {
		if (line>=skiplines)
		{
			for (unsigned int i = 0; i < (*loop).size(); i++ )
			{
				std::string s = (*loop)[i];
				removeCharsFromString( s, ";" );
				In.at<double>(line-skiplines,i) =  std::stod(s);
			}
		}
		line++;		
    }
	return In;
}
/*--------------------------------------------------------------------------*/
bool sort_by_C_value (const  Points_With_Value &lhs, const Points_With_Value &rhs)
{
    return lhs.C_value > rhs.C_value;
}
/*--------------------------------------------------------------------------*/
static void compute_Save_GridX_Y(const cv::Size &Size, const unsigned int &xStart_ROI, const unsigned int &yStart_ROI, const unsigned int &GridLength, const std::string path)
{
    cv::Mat GridX(Size, CV_64FC1, Scalar(0));
    cv::Mat GridY(Size, CV_64FC1, Scalar(0));
	for (unsigned int i = 0; i < static_cast<unsigned int>(GridX.cols); i++)
	{
		for (unsigned int j = 0; j < static_cast<unsigned int>(GridX.rows); j++)
		{
			GridX.at<double>(j,i) = xStart_ROI + i*GridLength;
			GridY.at<double>(j,i) = yStart_ROI + j*GridLength;
		}
	}
    store_matrix(path,"GridX", GridX);
	store_matrix(path,"GridY", GridY);
}
/*--------------------------------------------------------------------------*/

