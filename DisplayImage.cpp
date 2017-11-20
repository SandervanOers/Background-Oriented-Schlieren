#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <math.h>
#include <numeric>
#include <algorithm>
#include <sstream>


#include "pixeltranslation.hpp"
#include "nonlineariteration.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}
using namespace cv;

struct Points_With_Value {
    double C_value;
    cv::Point  Loc;
    Points_With_Value(double k, cv::Point s) : C_value(k), Loc(s) {}
};
bool sort_by_C_value (const  Points_With_Value &lhs, const Points_With_Value &rhs) {
    return lhs.C_value > rhs.C_value;
    }

int main(int argc, char** argv )
{
    if ( argc != 7 )
    {
        printf("usage: DisplayImage.out <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction\n");
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
    // Stopping Criterion
    double abs_tolerance_threshold = 1e-8;
    double rel_tolerance_threshold = 1e-8;


    unsigned int offset = 2*((SplineDegree+1)/2) > 5 ? 2*((SplineDegree+1)/2) : 5;
    offset = GridLength > (SubsetLength/2+offset) ? GridLength-SubsetLength/2 : offset; //
    std::cout << "offset = " << offset << std::endl;
    /// Read all tif files in folder
    unsigned int xStart = 775;//750;
    unsigned int xEnd   = 1220;//1360;
    unsigned int yStart = 1;
    unsigned int yEnd   = 515;//550;
    unsigned int horx = xEnd - xStart;
    unsigned int very = yEnd - yStart;

    unsigned int xStart_ROI = xStart+SubsetLength/2+offset;
    unsigned int xEnd_ROI   = xEnd-SubsetLength/2-offset;
    unsigned int yStart_ROI = yStart+SubsetLength/2+offset;
    unsigned int yEnd_ROI   = yEnd-SubsetLength/2-offset;
    unsigned int horx_ROI = xEnd_ROI - xStart_ROI;
    unsigned int very_ROI = yEnd_ROI - yStart_ROI;

    cv::String path("./*.tif");
    std::vector<cv::String> filenames;
    std::vector<cv::Mat> data;

    cv::glob(path,filenames,true); // recurse
    for (size_t k=0; k<filenames.size(); ++k)
    {
         cv::Mat im = cv::imread(filenames[k] , IMREAD_GRAYSCALE );
         if (im.empty()) continue; //only proceed if sucsessful
         im = im(Range(yStart,yEnd),Range(xStart,xEnd));
         data.push_back(im);
    }

    cv::Mat img, img1;
    (data.at(0)).copyTo(img);

    cv::Mat M_valid_points(img.size(), CV_8U, Scalar(0));
    std::cout << "img.rows = " << img.rows << std::endl;
    std::cout << "img.cols = " << img.cols << std::endl;
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

    //(data.at(0)).copyTo(img);
    //(data.at(2)).copyTo(img1);
    (data.at(2)).copyTo(img);
    (data.at(0)).copyTo(img1);

    img.convertTo(img, CV_32F);
    img1.convertTo(img1, CV_32F);

    std::vector<cv::Mat> Displacements;
    unsigned int GridLengthRough = 80;
    std::cout << horx_ROI << ", " << img.cols-SubsetLength/2-SubsetLength/2-offset*2 << std::endl;
    std::cout << very_ROI << ", " << img.rows-SubsetLength/2-SubsetLength/2-offset*2 << std::endl;
    GridLengthRough = GridLength*(GridLengthRough/GridLength);
    Displacements = calculatePixelTranslation(img, img1, SubsetLength, GridLengthRough, offset);

    std::cout << "very_ROI = " << very_ROI << ", horx_ROI = " << horx_ROI << std::endl;
    cv::Mat DispX(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, Scalar(0));
    cv::Mat DispY(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, Scalar(0));
    cv::Mat CorrelationCoefficient(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, Scalar(0));

    std::cout << Displacements[0].size() << std::endl;

    std::cout << Displacements[0] << std::endl;
    std::cout << Displacements[1] << std::endl;
    std::cout << Displacements[2] << std::endl;
    for (unsigned int i = 0; i < (Displacements.at(0)).cols; i++)
    {
        for (unsigned int j = 0; j < (Displacements.at(0)).rows; j++)
        {
            DispX.at<double>(j*GridLengthRough/GridLength,i*GridLengthRough/GridLength) = (Displacements.at(0)).at<double>(j,i);
            DispY.at<double>(j*GridLengthRough/GridLength,i*GridLengthRough/GridLength) = (Displacements.at(1)).at<double>(j,i);
            CorrelationCoefficient.at<double>(j*GridLengthRough/GridLength,i*GridLengthRough/GridLength) = (Displacements.at(2)).at<double>(j,i);
        }
    }

    /// Construct an fm-table and df-table
    // fm = mean of grey-level values of all points in a subset
    // for each grid point, take the average of the points the subset
    // df = sqrt( sum ( (f - fm) ^ 2 )

    // Calculate interpolation coefficients for g
    cv::Mat prova_img1= img1.clone();
    float *fptr_img1 = prova_img1.ptr<float>(0);
    //Mat gcoeff = Mat(img1.size(), CV_32F, (float*)fptr_img1).clone();
    SamplesToCoefficients(fptr_img1, img1.cols, img1.rows, SplineDegree);
    //Mat gcoeff1 = Mat(img1.size(), CV_32F, (float*)fptr_img1).clone();

    //std::cout << gcoeff - gcoeff1 << std::endl;


    // Matrix to keep track of which points we have computed
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
        // Do a full solve
        std::vector<double> InitialCondition = {DispX.at<double>(matchLoc), DispY.at<double>(matchLoc), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<double> point1 = iteration(img, fptr_img1, matchLoc.x, matchLoc.y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
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
        /*
        std::vector<double> point1;
        switch (ShapeFunction)
        {
        case 0:
            point1 = iteration_rigid_LM(img, fptr_img1, matchLoc.x, matchLoc.y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
            break;
        case 1:
            point1 = iteration_affine_LM(img, fptr_img1, matchLoc.x, matchLoc.y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
            Ux.at<double>(matchLoc) = point1[2];
            Vx.at<double>(matchLoc) = point1[3];
            Uy.at<double>(matchLoc) = point1[4];
            Vy.at<double>(matchLoc) = point1[5];
            break;
        case 2:
            point1 = iteration_irregular_LM(img, fptr_img1, matchLoc.x, matchLoc.y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
            Ux.at<double>(matchLoc) = point1[2];
            Vx.at<double>(matchLoc) = point1[3];
            Uy.at<double>(matchLoc) = point1[4];
            Vy.at<double>(matchLoc) = point1[5];
            Uxy.at<double>(matchLoc) = point1[6];
            Vxy.at<double>(matchLoc) = point1[7];
            break;
        case 3:
            point1 = iteration_quadratic_LM(img, fptr_img1, matchLoc.x, matchLoc.y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
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
            break;
        default:
            std::cout << "Shape Function Incorrect" << std::endl;
            break;
        }
        DispX.at<double>(matchLoc) = point1[0];
        DispY.at<double>(matchLoc) = point1[1];
        CorrelationCoefficient.at<double>(matchLoc) = point1.back();
        */
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

            /*
            std::vector<double> point2;
            // Use Initial Guess from (neighbouring) initial point

            switch (ShapeFunction)
            {
            case 0:
                point2 = iteration_rigid_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
                break;
            case 1:
                point2 = iteration_affine_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);

                Ux.at<double>(*i) = point2[2];
                Vx.at<double>(*i) = point2[3];
                Uy.at<double>(*i) = point2[4];
                Vy.at<double>(*i) = point2[5];
                break;
            case 2:
                point2 = iteration_irregular_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
                Ux.at<double>(*i) = point2[2];
                Vx.at<double>(*i) = point2[3];
                Uy.at<double>(*i) = point2[4];
                Vy.at<double>(*i) = point2[5];
                Uxy.at<double>(*i) = point2[6];
                Vxy.at<double>(*i) = point2[7];
                break;
            case 3:
                point2 = iteration_quadratic_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
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
                break;
            default:
                std::cout << "Shape Function Incorrect" << std::endl;
                break;
            }
            DispX.at<double>(*i) = point2[0];
            DispY.at<double>(*i) = point2[1];
            CorrelationCoefficient.at<double>(*i) = point2.back();
            */
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

    int64 t0 = cv::getTickCount();

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
            std::vector<double> point2 = iteration(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction);
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
            /*
            std::vector<double> point2;
            switch (ShapeFunction)
            {
            case 0:
                point2 = point2 = iteration_rigid_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
                break;
            case 1:
                point2 = iteration_affine_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
                Ux.at<double>(*i) = point2[2];
                Vx.at<double>(*i) = point2[3];
                Uy.at<double>(*i) = point2[4];
                Vy.at<double>(*i) = point2[5];
                break;
            case 2:
                point2 = iteration_irregular_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
                Ux.at<double>(*i) = point2[2];
                Vx.at<double>(*i) = point2[3];
                Uy.at<double>(*i) = point2[4];
                Vy.at<double>(*i) = point2[5];
                Uxy.at<double>(*i) = point2[6];
                Vxy.at<double>(*i) = point2[7];
                break;
            case 3:
                point2 = iteration_quadratic_LM(img, fptr_img1, (*i).x, (*i).y, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
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
                break;
            default:
                std::cout << "Shape Function Incorrect" << std::endl;
                break;
            }

            DispX.at<double>(*i) = point2[0];
            DispY.at<double>(*i) = point2[1];
            CorrelationCoefficient.at<double>(*i) = point2.back();
            */
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
            std::cout << "Points Computed: " << cv::sum(Computed_Points).val[0]/Computed_Points.total()*100 << "%" << std::endl;
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
        }
    }


    unsigned int XStart = 0;
    if (GridLength==3) XStart = 45;
    if (GridLength==5) XStart = 30;

    std::cout << "Reliability = " << cv::sum(CorrelationCoefficient(Range::all(), Range(XStart, CorrelationCoefficient.cols))).val[0]/CorrelationCoefficient(Range::all(), Range(XStart, CorrelationCoefficient.cols)).total() << std::endl;
    std::cout << "Reliability = " << cv::sum(CorrelationCoefficient).val[0]/CorrelationCoefficient.total() << std::endl;

    //std::cout << "Reliability Small = " << cv::sum(CorrelationCoefficient(Range::all(), Range(0, 5.0/7.0*CorrelationCoefficient.cols))).val[0]/CorrelationCoefficient(Range::all(), Range(0, 5.0/7.0*CorrelationCoefficient.cols)).total() << std::endl;

    int64 t1 = cv::getTickCount();
    double secs = (t1-t0)/cv::getTickFrequency();

    std::cout << "Time taken = " << secs << " seconds" << std::endl;
//    if (1)
//    {
//        std::stringstream ss;
//        ss<<"CC_final_8_"<<(ShapeFunction)<<type;
//        std::string filename = ss.str();
//        ss.str("");
//        cv::Mat Copy_CCF_8 = CorrelationCoefficient_Final.clone();
//        Copy_CCF_8.convertTo(Copy_CCF_8, CV_8U, 255.0);
//        imwrite(filename, Copy_CCF_8);
//    }
    if (1)
    {
        std::stringstream ss;
        ss<<folderName<<"CC_PF"<<(propagation_function)<<"_Bn"<<(SplineDegree)<<"_SF"<<(ShapeFunction)<<"_GL"<<(GridLength)<<"_SL"<<(SubsetLength)<<type;
        std::string filename = ss.str();
        ss.str("");
        cv::Mat Copy_CCF_16 = CorrelationCoefficient.clone();
        Copy_CCF_16.convertTo(Copy_CCF_16, CV_16U, 255.0*256.0);
        imwrite(filename, Copy_CCF_16);
    }

















}
