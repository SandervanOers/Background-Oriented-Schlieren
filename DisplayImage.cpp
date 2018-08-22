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
//#include <pthread.h>
#include <thread>

#include "pixeltranslation.hpp"
#include "nonlineariteration.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}
using namespace cv;

void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored)
{
    std::ofstream myfile;
    myfile.open(path+"/"+filename+".csv");
    myfile << format(Matrix_To_Be_Stored,cv::Formatter::FMT_MATLAB);
    myfile.close();
};
struct Points_With_Value {
    double C_value;
    cv::Point  Loc;
    Points_With_Value(double k, cv::Point s) : C_value(k), Loc(s) {}
};
bool sort_by_C_value (const  Points_With_Value &lhs, const Points_With_Value &rhs) {
    /*if (isnan(lhs.C_value) || isnan(rhs.C_value))
    {
        std::cout << "lhs.C_value = " << lhs.C_value <<", rhs.C_value = " << rhs.C_value << std::endl;
        std::cout << "lhs.C_value > rhs.C_value : " << (lhs.C_value > rhs.C_value) << std::endl;
    }*/
    return lhs.C_value > rhs.C_value;
    }
void *worker_thread(void *arg)
{
        printf("This is worker_thread #%ld\n", (long)arg);
        pthread_exit(NULL);
}
int main(int argc, char** argv )
{
    if ( argc != 13)
    {
        printf("usage: DisplayImage.out <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction OrderingImages xStart xEnd yStart yEnd NumberOfThreads\n");
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
    if (!(SplineDegree > 1 && SplineDegree < 9))
    {
        printf("Spline degree must be between 1 and 9 \n");
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
	
	unsigned int NumberOfThreads = atoi(argv[12]);
	if (NumberOfThreads > std::thread::hardware_concurrency())
	{
		std::cout << "Specified Number of Threads larger than maximum available on this system. Using the system's maximum" << std::endl;
		NumberOfThreads = std::thread::hardware_concurrency();
	}
	
    // Stopping Criterion
    double abs_tolerance_threshold = 1e-8;
    double rel_tolerance_threshold = 1e-8;

    unsigned int offset = 2*((SplineDegree+1)/2) > 5 ? 2*((SplineDegree+1)/2) : 5;
    offset = GridLength > (SubsetLength/2+offset) ? GridLength-SubsetLength/2 : offset; //
    //std::cout << "offset = " << offset << std::endl;
    
	cv::String path(pathname);
    std::vector<cv::String> filenames;
    std::vector<cv::Mat> data;
	cv::String path1("Images/*.tif");
	std::cout << path1 << std::endl;
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
			if (xEnd > im.cols || yEnd > im.rows)
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
    for (unsigned int i = 0; i < img.rows; i++)
    {
        for (unsigned int j = 0; j < img.cols; j++)
        {
            if (i >= SubsetLength/2+offset && i <= img.rows-SubsetLength/2-offset && j >= SubsetLength/2+offset && j <= img.cols-SubsetLength/2-offset)
            {
                M_valid_points.at<uchar>(i,j) = 1;
            }
        }
    }
    img.convertTo(img, CV_32F);
    img1.convertTo(img1, CV_32F);

    std::vector<cv::Mat> Displacements;
    std::cout << "very_ROI = " << very_ROI << ", horx_ROI = " << horx_ROI << std::endl;
	std::vector<cv::Mat> IG = calculateInitialGuess(img, img1, SubsetLength, GridLength, horx_ROI, very_ROI, offset, NumberOfThreads); 
	cv::Mat DispX =  IG[0].clone();
	cv::Mat DispY =  IG[1].clone();
	cv::Mat CorrelationCoefficient =  IG[2].clone();

    /// Construct an fm-table and df-table
    // fm = mean of grey-level values of all points in a subset
    // for each grid point, take the average of the points the subset
    // df = sqrt( sum ( (f - fm) ^ 2 )


	auto tr1 = std::chrono::high_resolution_clock::now();
    // Calculate interpolation coefficients for g
    cv::Mat prova_img1= img1.clone();
    float *fptr_img1 = prova_img1.ptr<float>(0);
    SamplesToCoefficients(fptr_img1, img1.cols, img1.rows, SplineDegree);
    //std::cout << gcoeff - gcoeff1 << std::endl;

    cv::Mat GridX(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat GridY(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Computed_Points(DispX.size(), CV_8UC1, Scalar(0));
    cv::Mat Ux(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Vx(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Uy(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Vy(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Uxy(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Vxy(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Uxx(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Vxx(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Uyy(DispX.size(), CV_64FC1, Scalar(0));
    cv::Mat Vyy(DispX.size(), CV_64FC1, Scalar(0));

    // Find the point (seed point) which is best described by (pixel-integer) rigid translation
    // We do this for both the upper half and lower half for the two layer case
    std::vector<Points_With_Value> Locations_Best_Correlation;
    for (unsigned int k = 0; k < 2; k++)
    {
        double minVal; double maxVal;
        cv::Point minLoc; cv::Point maxLoc;
        cv::Point matchLoc;
        cv::minMaxLoc(CorrelationCoefficient(Range(k*(CorrelationCoefficient.rows/2+1),k*(CorrelationCoefficient.rows/2)+CorrelationCoefficient.rows/2),Range::all()), &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );
        matchLoc = maxLoc;
        matchLoc.y += k*(CorrelationCoefficient.rows/2+1);
        std::cout << "Maximum found at " << matchLoc << ", Value = " << maxVal << ", Initial Guess U0 = " << DispX.at<double>(matchLoc) << ", Initial Guess V0 = " << DispY.at<double>(matchLoc) << std::endl;
        while (cv::abs(DispY.at<double>(matchLoc))>50 )
        {
            std::cout << "Maximum Rejected" << std::endl;
            CorrelationCoefficient.at<double>(matchLoc) = 0;
            cv::minMaxLoc(CorrelationCoefficient(Range(k*(CorrelationCoefficient.rows/2+1),k*(CorrelationCoefficient.rows/2)+CorrelationCoefficient.rows/2),Range::all()), &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );
            matchLoc = maxLoc;
            matchLoc.y += k*(CorrelationCoefficient.rows/2+1);
            std::cout << "New Maximum found at " << matchLoc << ", Value = " << maxVal << ", Initial Guess U0 = " << DispX.at<double>(matchLoc) << ", Initial Guess V0 = " << DispY.at<double>(matchLoc) << std::endl;
        }
        // Do a full solve
        std::vector<double> InitialCondition = {DispX.at<double>(matchLoc), DispY.at<double>(matchLoc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<double> point1 = iteration(img, fptr_img1, matchLoc.x, matchLoc.y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
        GridX.at<double>(matchLoc) = xStart_ROI+matchLoc.x*GridLength;
        GridY.at<double>(matchLoc) = yStart_ROI+matchLoc.y*GridLength;
        DispX.at<double>(matchLoc) = point1[0];
        DispY.at<double>(matchLoc) = point1[1];
        Ux.at<double>(matchLoc) = point1[2];
        Vx.at<double>(matchLoc) = point1[3];
        Uy.at<double>(matchLoc) = point1[4];
        Vy.at<double>(matchLoc) = point1[5];
        Uxy.at<double>(matchLoc) = point1[6];
        Vxy.at<double>(matchLoc) = point1[7];
        Uxx.at<double>(matchLoc) = point1[8];
        Vxx.at<double>(matchLoc) = point1[9];
        Uyy.at<double>(matchLoc) = point1[10];
        Vyy.at<double>(matchLoc) = point1[11];
        CorrelationCoefficient.at<double>(matchLoc) = point1.back();
        Computed_Points.at<uchar>(matchLoc) = 1;

        std::cout << "Solution: " ;
        for (auto i = point1.begin(); i < point1.end(); i++)
            std::cout << *i << " ";
        std::cout << std::endl;
        InitialCondition = {DispX.at<double>(matchLoc), DispY.at<double>(matchLoc), Ux.at<double>(matchLoc), Vx.at<double>(matchLoc), Uy.at<double>(matchLoc), Vy.at<double>(matchLoc), Uxy.at<double>(matchLoc), Vxy.at<double>(matchLoc), Uxx.at<double>(matchLoc), Vxx.at<double>(matchLoc), Uyy.at<double>(matchLoc), Vyy.at<double>(matchLoc)};
		std::vector<cv::Point> Neighbours1 = get_valid_Neighbours(M_valid_points, Computed_Points, matchLoc.x, matchLoc.y, SubsetLength, GridLength, offset);
		for (auto i = Neighbours1.begin(); i < Neighbours1.end(); i++)
        {
            std::cout << *i << std::endl;
                    std::vector<double> point2 = iteration(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
                    GridX.at<double>(*i) = xStart_ROI+(*i).x*GridLength;
                    GridY.at<double>(*i) = yStart_ROI+(*i).y*GridLength;
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
            std::cout << "Solution: " ;
            for (auto s = point2.begin(); s < point2.end(); s++)
                std::cout << *s << " ";
            std::cout << std::endl;
            // Add Location to queue
            Locations_Best_Correlation.push_back(Points_With_Value(point2.back(), *i));
        }
    }

    std::cout << "Computed Neighbours of Initial Points" << std::endl;

    // Sort first (eight) elements
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
		auto tr2 = std::chrono::high_resolution_clock::now();
	    std::cout << "Function IG took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(tr2-tr1).count()
              << " milliseconds" << std::endl;
    std::cout << "Start Iteration" << std::endl;	
   while(cv::sum(Computed_Points).val[0] < Computed_Points.total())
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
            //std::cout << *i << std::endl;
            std::vector<double> point2 = iteration(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
            GridX.at<double>(*i) = xStart_ROI+(*i).x*GridLength;
            GridY.at<double>(*i) = yStart_ROI+(*i).y*GridLength;
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
        /*
        for (auto l = Locations_Best_Correlation.begin(); l < Locations_Best_Correlation.end();l++)
        {
            if (isnan((*l).C_value))
            {
                    std::cout << "New List " << std::endl;
                    for (auto k = Locations_Best_Correlation.begin(); k < Locations_Best_Correlation.end();k++)
                           std::cout << (*k).Loc << ": " << (*k).C_value << std::endl;
                    std::cout << "End List " << std::endl;

            }
        }
        */
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
                ///cv::Mat Copy_CCF = CorrelationCoefficient.clone();
                ///Copy_CCF.convertTo(Copy_CCF, CV_8UC1, 255.0);
                ///imwrite(filename, Copy_CCF);
            }
            else
            {
                double computed = cv::sum(Computed_Points).val[0]/Computed_Points.total()*100;
                if (computed>outputcount)
                {
                    outputcount++;
                    outputcount++;
                    std::cout << "Points Computed: " << std::setprecision(2) << computed << "%" << std::endl;
                }
            }
        }
    }
	
	std::cout << "Reliability = " << std::setprecision(5) << cv::sum(CorrelationCoefficient).val[0]/CorrelationCoefficient.total() << std::endl;
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
	
    std::cout << Ux.size() << std::endl;
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
	 
    return 0;
}