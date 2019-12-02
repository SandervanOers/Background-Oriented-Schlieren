#include "InputOut.hpp"

/*--------------------------------------------------------------------------*/
extern int checkinput(unsigned int argc, char *argv[], std::string &pathname_string, unsigned int &SplineDegree, unsigned int &SubsetLength, unsigned int &GridLength, unsigned int &ShapeFunction, unsigned int &propagationfunction, unsigned int &ordering, unsigned int &xStart, unsigned int &xEnd, unsigned int &yStart, unsigned int &yEnd, unsigned int &offset, unsigned int &Number_Of_Threads, unsigned int &MaxPixelYVertical, double &abs_tolerance_threshold, double &rel_tolerance_threshold, double &minimum_corrcoeff_IG )
{
	// Checks number of Arguments
    if ( argc != 16)
    {
		printf("Too few arguments\n");
        printf("usage: DisplayImage.out <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction OrderingImages xStart xEnd yStart yEnd NumberOfThreads MaxPixelYVertical Tolerance MimCorrCoeffIG\n");
		return -1;
    }

	// Checks pathname
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
	pathname_string = std::string(pathname);

	// Saves Command
	std::ofstream myfile;
	myfile.open (pathname_string+"/commandline.txt");
	myfile << "usage: DisplayImage.out <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction OrderingImages xStart xEnd yStart yEnd NumberOfThreads MaxPixelYVertical Tolerance\n";
	for (unsigned int i=0;i<15;i++)
	{
		myfile << argv[i] << " ";
	}
	myfile.close();

	// Checks input
    SplineDegree = atoi(argv[2]);
    if (!(SplineDegree > 2 && SplineDegree < 9))
    {
        printf("Spline degree must be between 2 and 9 \n");
        return -1;
    }
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
    ShapeFunction = atoi(argv[5]);
    if (!(ShapeFunction==0 || ShapeFunction==1 || ShapeFunction==2 || ShapeFunction==3))
    {
        printf("Shape function must be 0, 1, 2 or 3 \n");
        return -1;
    }
    if (SubsetLength < GridLength)
    {
        printf("Warning: Subset Length smaller than GridLength: you are not using all available data. \n");
    }
    propagationfunction = atoi(argv[6]);
    if (propagationfunction != 0 && propagationfunction != 1)
    {
        printf("Propagation Function must be on (1) or off (0)");
        return -1;
    }
	ordering = atoi(argv[7]);
    if (ordering != 0 && ordering != 1)
    {
        printf("Ordering must be natural (0) or reverse (1)");
        return -1;
    }
	xStart = atoi(argv[8]);
	xEnd = atoi(argv[9]);
	yStart = atoi(argv[10]);
	yEnd = atoi(argv[11]);
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
    offset = 2*((SplineDegree+1)/2) > 5 ? 2*((SplineDegree+1)/2) : 5;
	std::cout << " offset = " << offset << std::endl;
    offset = GridLength > (SubsetLength/2+offset) ? GridLength-SubsetLength/2 : offset;
    std::cout << "offset = " << offset << std::endl;
    offset++;
    std::cout << "offset = " << offset << std::endl;
    std::cout << " SubsetLength/2 = "  << SubsetLength/2 << std::endl;
	if ((yEnd-yStart) < SubsetLength/2+offset)
	{
		std::cout << "Vertical Range of Image is too small with this Subset" << std::endl;
		return -1;
	}
	if ((xEnd-xStart) < SubsetLength/2+offset)
	{
		std::cout << "Horizontal Range of Image is too small with this Subset" << std::endl;
		return -1;
	}
	Number_Of_Threads = atoi(argv[12]);
	if (Number_Of_Threads > std::thread::hardware_concurrency())
	{
		std::cout << "Specified Number of Threads larger than maximum available on this system. Using the system's maximum" << std::endl;
		Number_Of_Threads = std::thread::hardware_concurrency();
	}
	MaxPixelYVertical = atoi(argv[13]);
	if (MaxPixelYVertical > (yEnd-yStart))
	{
		std::cout << "Maximum Vertical Pixel Allowed is larger than Vertical Image Size"<< std::endl;
		std::cout << "TO DO: check if this is an error"<< std::endl;
		return -1;
	}
	double tolerance = atof(argv[14]);
	if (tolerance < 0)
	{
		std::cout << "Tolerance is negative"<< std::endl;
		return -1;
	}
	if (tolerance >= 1)
	{
		std::cout << "Tolerance is too large" << std::endl;
		return -1;
	}
    // Stopping Criterion
    abs_tolerance_threshold = tolerance;
    rel_tolerance_threshold = tolerance;
	// Minimum Acceptable Correlation Coefficient for Initial Guess
	minimum_corrcoeff_IG =  atoi(argv[15]);

	return 1;
}
/*--------------------------------------------------------------------------*/
extern int readImageDataFromFile(const cv::String &path, cv::Mat &img, cv::Mat &img1, const unsigned int &xStart, const unsigned int &xEnd, const unsigned int &yStart, const unsigned int &yEnd, const unsigned int &SubsetLength, const unsigned int &offset, unsigned int &xStart_ROI, unsigned int &yStart_ROI, unsigned int &horx_ROI, unsigned int &very_ROI, const unsigned int &ordering)
{
    std::vector<cv::String> filenames;
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
	std::vector<cv::Mat> data;
	// Read all tif files in folder
    for (size_t k=0; k<filenames.size(); ++k)
    {
         cv::Mat im = cv::imread(filenames[k] , IMREAD_GRAYSCALE );
         if (im.empty()) continue; //only proceed if sucsessful
		 // Change Here
         im = im(Range(yStart,yEnd),Range(xStart,xEnd));
         data.push_back(im);
         std::cout << filenames[k] << std::endl;
    }
	std::cout << std::endl;
    xStart_ROI = xStart+SubsetLength/2+offset;
    unsigned int xEnd_ROI   = xEnd-SubsetLength/2-offset;
    yStart_ROI = yStart+SubsetLength/2+offset;
    unsigned int yEnd_ROI   = yEnd-SubsetLength/2-offset;
    horx_ROI = xEnd_ROI - xStart_ROI;
    very_ROI = yEnd_ROI - yStart_ROI;

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
	/*
	// Restrict Template to specified Region
	img = img(Range(yStart,yEnd),Range(xStart,xEnd));
	// We do not want to restrict the Deformed Image (we might miss data)
	// TO DO: Need to rewrite nonlineariteration.cpp to use img1.cols and img1.rows instead of img.cols and img.rows since interpol.c needs width and height of image.
	img1 = img1(Range(yStart,yEnd),Range(xStart,xEnd));
	 * */
    imwrite(path+"/Template.png", img);
    imwrite(path+"/Deformed.png", img1);


    img.convertTo(img, CV_32F);
    img1.convertTo(img1, CV_32F);
	return 1;
}
/*--------------------------------------------------------------------------*/
