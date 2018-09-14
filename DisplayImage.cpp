#include "DisplayImage.hpp"
void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored)
{
    std::ofstream myfile;
    myfile.open(path+"/"+filename+".csv");
    myfile << format(Matrix_To_Be_Stored,cv::Formatter::FMT_MATLAB);
    myfile.close();
}
bool sort_by_C_value (const  Points_With_Value &lhs, const Points_With_Value &rhs) 
{
    return lhs.C_value > rhs.C_value;
}
int main(int argc, char** argv )
{
	auto tr1 = std::chrono::high_resolution_clock::now();
	for (unsigned int kk=0; kk <1; kk++)
	{
    if ( argc != 15)
    {
        printf("usage: DisplayImage.out <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction OrderingImages xStart xEnd yStart yEnd NumberOfThreads MaxPixelYVertical Tolerance\n");
        return -1;
    }
	
	const char *pathname;
	pathname = argv[1];
	struct stat info;

	if( stat( pathname, &info ) != 0 )
	{
		printf( "cannot access %s\n", pathname );
		return -1;
	}
	else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows 
	{
		//printf( "%s is a directory\n", pathname );
	}
	else
	{
		printf( "%s is no directory\n", pathname );
		return -1;
	}
	
    unsigned int SplineDegree;
    SplineDegree = atoi(argv[2]);
    if (!(SplineDegree > 2 && SplineDegree < 9))
    {
        printf("Spline degree must be between 2 and 9 \n");
        return -1;
    }
    unsigned int SubsetLength;
    unsigned int GridLength;
    SubsetLength = atoi(argv[3]);
    if (!(SubsetLength > 1))
    {
        printf("Subset Length must be larger than 1 \n");
        return -1;
    }
    GridLength = atoi(argv[4]);
    if (!(GridLength > 0))
    {
        printf("GridLength must be larger than 0 \n");
        return -1;
    }
    unsigned int ShapeFunction = atoi(argv[5]);
    if (!(ShapeFunction==0 || ShapeFunction==1 || ShapeFunction==2 || ShapeFunction==3))
    {
        printf("Shape function must be 0, 1, 2 or 3 \n");
        return -1;
    }
    if (SubsetLength < GridLength)
    {
        printf("Warning: Subset Length smaller than GridLength: you are not using all available data. \n");
    }
    unsigned int propagation_function = atoi(argv[6]);
    if (propagation_function != 0 && propagation_function != 1)
    {
        printf("Propagation Function must be on (1) or off (0)");
        return -1;
    }
    unsigned int ordering = atoi(argv[7]);
    if (ordering != 0 && ordering != 1)
    {
        printf("Ordering must be natural (0) or reverse (1)");
        return -1;
    }
	unsigned int xStart = atoi(argv[8]);
	unsigned int xEnd = atoi(argv[9]);
	unsigned int yStart = atoi(argv[10]);
	unsigned int yEnd = atoi(argv[11]);
	if (xStart > xEnd)
	{
		std::cout << "xStart larger than xEnd" << std::endl;
		return -1;
	}
	if (yStart > yEnd)
	{
		std::cout << "yStart larger than yEnd" << std::endl;
		return -1;
	}
    unsigned int offset = 2*((SplineDegree+1)/2) > 5 ? 2*((SplineDegree+1)/2) : 5;
    offset = GridLength > (SubsetLength/2+offset) ? GridLength-SubsetLength/2 : offset; 
	if ((yEnd-yStart) < SubsetLength/2+offset)
	{
		std::cout << "Vertical Range of Image is too small with this Subset" << std::endl;
	}	
	if ((xEnd-xStart) < SubsetLength/2+offset)
	{
		std::cout << "Horizontal Range of Image is too small with this Subset" << std::endl;
	}
	unsigned int Number_Of_Threads = atoi(argv[12]);
	if (Number_Of_Threads > std::thread::hardware_concurrency())
	{
		std::cout << "Specified Number of Threads larger than maximum available on this system. Using the system's maximum" << std::endl;
		Number_Of_Threads = std::thread::hardware_concurrency();
	}	
	unsigned int MaxPixelYVertical = atoi(argv[13]); 
	double tolerance = atof(argv[14]);
	if (tolerance < 0)
	{
		std::cout << "Tolerance is negative"<< std::endl;
		return -1;
	}
	if (tolerance > 1)
	{
		std::cout << "Tolerance is too large" << std::endl;
		return -1;
	}
	// Minimum Acceptable Correlation Coefficient for Initial Guess
	double minimum_corrcoeff = 0.8;
    // Stopping Criterion
    double abs_tolerance_threshold = tolerance;
    double rel_tolerance_threshold = tolerance;
    
	cv::String path(pathname);
    std::vector<cv::String> filenames;
    std::vector<cv::Mat> data;
	cv::String path1;
	path1 = path+"/*.tif";
    cv::glob(path1,filenames,true); // recurse
	// Read First Image to Determine Size 
	{
		cv::Mat im = cv::imread(filenames[0] , IMREAD_GRAYSCALE );
		if (im.empty()) 
		{
			std::cout << "No Images Found" << std::endl;
			return -1;
		}
		else
		{
			if (xEnd > static_cast<unsigned int>(im.cols) || yEnd > static_cast<unsigned int>(im.rows))
			{
				std::cout << "Region of Interest (ROI) larger than Image" << std::endl;
				std::cout << "Image Size = " << im.size() << std::endl;
				return -1;
			}	
		}
	}
	// Read all tif files in folder	
    for (size_t k=0; k<filenames.size(); ++k)
    {
         cv::Mat im = cv::imread(filenames[k] , IMREAD_GRAYSCALE );
         if (im.empty()) continue; //only proceed if sucsessful
         im = im(Range(yStart,yEnd),Range(xStart,xEnd));
         data.push_back(im);
         std::cout << filenames[k] << std::endl;
    }
			  
    cv::Mat img, img1;
    // Lyon
    //unsigned int xStart = 1950; //2080
    //unsigned int xEnd   = 3210;
    //unsigned int yStart = 270;
    //unsigned int yEnd   = 2180;
	
    //unsigned int horx = xEnd - xStart;
   // unsigned int very = yEnd - yStart;
    unsigned int xStart_ROI = xStart+SubsetLength/2+offset;
    unsigned int xEnd_ROI   = xEnd-SubsetLength/2-offset;
    unsigned int yStart_ROI = yStart+SubsetLength/2+offset;
    unsigned int yEnd_ROI   = yEnd-SubsetLength/2-offset;
    unsigned int horx_ROI = xEnd_ROI - xStart_ROI;
    unsigned int very_ROI = yEnd_ROI - yStart_ROI;

    std::cout << "xStart_ROI = " << xStart_ROI << std::endl;
    std::cout << "yStart_ROI = " << yStart_ROI << std::endl;
    std::cout << "horx_ROI = " << horx_ROI << std::endl;
    std::cout << "very_ROI = " << very_ROI << std::endl;
	
    if (ordering==0)
    {
        (data.at(0)).copyTo(img);
        (data.at(1)).copyTo(img1);
    }
    else if (ordering==1)
    {
        (data.at(1)).copyTo(img);
        (data.at(0)).copyTo(img1);
    }
    imwrite("show.tif", img);
    imwrite("show1.tif", img1);
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
    img.convertTo(img, CV_32F);
    img1.convertTo(img1, CV_32F);

    // Calculate interpolation coefficients for g
    cv::Mat prova_img1= img1.clone();
    float *fptr_img1 = prova_img1.ptr<float>(0);
    SamplesToCoefficients(fptr_img1, img1.cols, img1.rows, SplineDegree);
	
	std::vector<cv::Mat> IG = calculateInitialGuess_Iteration(img, img1, fptr_img1, SplineDegree, SubsetLength, GridLength, horx_ROI, very_ROI, offset, Number_Of_Threads, MaxPixelYVertical, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction, minimum_corrcoeff);
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
        std::cout << (*i).Loc << ": " << (*i).C_value << std::endl;

    const bool plotting = 0;
    int ct = 0;
    int outputcount = 1;
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
	auto tr5 = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> TR = tr5-tr5;
	auto T = TR.count();
	
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
			auto tr3 = std::chrono::high_resolution_clock::now();
            std::vector<double> point2 = iteration(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
			auto tr4 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = tr4-tr3;	
			T += elapsed_seconds.count();
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
            CorrelationCoefficient.at<double>(*i) = point2.back();
            Computed_Points.at<uchar>(*i) = 1;
            Locations_Best_Correlation.push_back(Points_With_Value(point2.back(), *i));
			
        }
		 
        Locations_Best_Correlation.erase(Locations_Best_Correlation.begin());
        sort(Locations_Best_Correlation.begin(), Locations_Best_Correlation.end(), sort_by_C_value);

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
                    std::cout << "Points Computed: " << std::setprecision(2) << computed << "%" << std::endl;
                }
            }
        } 
    }		/*
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
	*/
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
	
    store_matrix(path,"GridX", GridX);
	store_matrix(path,"GridY", GridY);
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
	
        cv::Mat Copy_CCF_16 = CorrelationCoefficient.clone();
        Copy_CCF_16.convertTo(Copy_CCF_16, CV_16U, 255.0*256.0);
        imwrite("Images/CC.png", Copy_CCF_16);
		
		flip(Copy_CCF_16, Copy_CCF_16, 0);
		imwrite("Images/CCflipped.png", Copy_CCF_16);
		
	//Mat Temp;
	//DispX.convertTo(Temp, CV_8U);
	//Mat DX_Median, DY_Median;
	//DX_Median = Temp.clone();
	//DispY.convertTo(Temp, CV_8U);
	//DY_Median = Temp.clone();
	////Apply median filter
    //medianBlur(Temp, DX_Median, 3 );
    //medianBlur(DispY, DY_Median, 3 );
	//store_matrix(path,"U0filter", DX_Median);
	// store_matrix(path,"V0filter", DY_Median);
	 
	auto tr2= std::chrono::high_resolution_clock::now();
	std::cout << "Total took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr2-tr1).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr2-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr2-tr1).count() << " minutes"<<std::endl;
	 	
		//std::cout << "Iterations took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr5).count()
		//<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr5-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr5-tr1).count() << " minutes"<<std::endl;
	 	 std::cout << "Iterations took "<< T << " seconds" << std::endl;
	 }
    return 0;
}