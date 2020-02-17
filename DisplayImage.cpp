#include "DisplayImage.hpp"

int main(int argc, char** argv )
{
	auto tr1 = std::chrono::high_resolution_clock::now();

	InputVariables inputvariables;
	int returnvalue = readinput(argc, argv, inputvariables);
	if (returnvalue < 0)
	{
		return -1;
	}
	std::cout << std::endl << "\033[1;32mInput Parsed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
	if (inputvariables.DICNeeded == 1)
	{
		// Read Images from Disk
		cv::Mat img, img1;
		returnvalue = readImageDataFromFile(img, img1, inputvariables);
		if (returnvalue < 0)
		{
			return -1;
		}
		std::cout << std::endl << "\033[1;32mImages Loaded\033[0m\n" << std::endl;
		/*--------------------------------------------------------------------------*/
		DIC(img, img1, inputvariables);
		auto tr2= std::chrono::high_resolution_clock::now();
		std::cout << "DIC took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr2-tr1).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr2-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr2-tr1).count() << " minutes\n"<<std::endl;
	}
	/*--------------------------------------------------------------------------*/
	cv::Mat DX = load_matrix(inputvariables.path, "U0", 1);
	cv::Mat DY = load_matrix(inputvariables.path, "V0", 1);
	cv::Mat GridX = load_matrix(inputvariables.path, "GridX", 1);
	cv::Mat GridY = load_matrix(inputvariables.path, "GridY", 1);
	cv::Mat CC = load_matrix(inputvariables.path, "CorrelationCoefficient", 1);
	std::cout << std::endl << "\033[1;32mLoading Matrices Completed\033[0m\n" << std::endl;   
	/*--------------------------------------------------------------------------*/
	if (inputvariables.ordering==1)
	{
		DX = - DX;
		DY = - DY;
	}
	DX = DX - 0.2;
	DY = DY + 0.3;
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
	GaussianBlur(DX, DX, Size(5, 5), 0);
	GaussianBlur(DY, DY, Size(5, 5), 0);
	// Store Again?
	std::cout << std::endl << "\033[1;32mFiltering Completed\033[0m\n" << std::endl;  
	/*--------------------------------------------------------------------------*/
	auto tr3= std::chrono::high_resolution_clock::now();
	ExperimentalSetupVariables experimentalsetupvariables;
	// Camera
	experimentalsetupvariables.focal_length = 50.0e-3;
	experimentalsetupvariables.Distance_From_Pixels_To_Meters = 6.45e-6;
	experimentalsetupvariables.n_0 = 1.0003;
	experimentalsetupvariables.n_1 = 1.52;
	experimentalsetupvariables.n = 1.333;
	// Lengths Small Tank
	experimentalsetupvariables.L_c = 1e4;
	experimentalsetupvariables.L_g = 0.8/100;//0.004; //0.01
	experimentalsetupvariables.L_t = 0.168;//0.191;//15.0/100-2*experimentalsetupvariables.L_g;//13.8/100-2*experimentalsetupvariables.L_g;//0.168;//
	experimentalsetupvariables.L_s = 0.0;//0.523;
	//experimentalsetupvariables.Lengths{experimentalsetupvariables.L_c, experimentalsetupvariables.L_t, experimentalsetupvariables.L_g, experimentalsetupvariables.L_s};
	// Lengths Full Tank
	/*
	double L_c = 1.60;
	double L_g = 5.8/1000;
	double L_t = 13.8/100-2*L_g;
	double L_s = 0.0;//0.523;
	 * */
	std::vector<double> Lengths{experimentalsetupvariables.L_c, experimentalsetupvariables.L_g, experimentalsetupvariables.L_t, experimentalsetupvariables.L_s};
	double corr_cut_off = 0.8;//
	/*--------------------------------------------------------------------------*/
	//CalibrationFigures2(GridX, GridY, DX, DY, CC, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, experimentalsetupvariables.n, inputvariables.path);
	/*--------------------------------------------------------------------------*/
	CalibrationValues calibrationValues;
	if (inputvariables.CalibrationNeeded==1)
	{
		std::vector<double> CalibrationNumbers = Calibration2(GridX, GridY, DX, DY, CC, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, experimentalsetupvariables.n, inputvariables.path, corr_cut_off);
		//std::vector<double> CalibrationNumbers = Calibration(GridX, GridY, DX, DY, CC, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, experimentalsetupvariables.n, inputvariables.path, corr_cut_off);
		calibrationValues.a = CalibrationNumbers[0];
		calibrationValues.b = CalibrationNumbers[1];
		calibrationValues.c = CalibrationNumbers[2];
		calibrationValues.L_m = CalibrationNumbers[3];
		calibrationValues.meanGridX = CalibrationNumbers[4];
		calibrationValues.meanGridY = CalibrationNumbers[5];
		std::cout << std::endl << "\033[1;32mCalibration Completed\033[0m\n" << std::endl;   

		for (auto i = CalibrationNumbers.begin(); i != CalibrationNumbers.end(); ++i)
			std::cout << *i << ' ';
			
		std::cout << std::endl<< std::endl<< std::endl;
		auto tr4= std::chrono::high_resolution_clock::now();
		std::cout << "Calibration took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr4-tr3).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr4-tr3).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr4-tr3).count() << " minutes"<<std::endl;
	}
	else
	{
		// Load Numbers
		cv::Mat CalibrationNumbers = load_matrix(inputvariables.path, "CalibrationFigures", 1);
		std::cout << "CalibrationNumbers = " << CalibrationNumbers << std::endl;
		calibrationValues.a = CalibrationNumbers.at<double>(0);
		calibrationValues.b = CalibrationNumbers.at<double>(1);
		calibrationValues.c = CalibrationNumbers.at<double>(2);
		calibrationValues.L_m = CalibrationNumbers.at<double>(3);
		calibrationValues.meanGridX = CalibrationNumbers.at<double>(4);
		calibrationValues.meanGridY = CalibrationNumbers.at<double>(5);
		std::cout << std::endl << "\033[1;32mReading Calibration Completed\033[0m\n" << std::endl;   
	}
	/*
	//std::vector<double> CalibrationNumbers = Calibration(GridX, GridY, DX, DY, CC, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, inputvariables.path, corr_cut_off);
	std::vector<double> CalibrationNumbers = Calibration2(GridX, GridY, DX, DY, CC, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, inputvariables.path, corr_cut_off);
	//std::vector<double> CalibrationNumbers = CalibrationExtended(GridX, GridY, DX, DY, CC, focal_length, Lengths, Distance_From_Pixels_To_Meters, n_0, n_1, n, inputvariables.path);
	*/


	/*--------------------------------------------------------------------------*/
	if (inputvariables.CalculateRefractionIndex == 1)
	{
	
		double normPD = calculateNorm(calibrationValues.a, calibrationValues.b, calibrationValues.c);
		double d = -calibrationValues.c/normPD*calibrationValues.L_m;
		std::vector<double> PlaneDefinition{0, 0, 0, 0};
		PlaneDefinition[0] = calibrationValues.a/normPD;
		PlaneDefinition[1] = calibrationValues.b/normPD;
		PlaneDefinition[2] = calibrationValues.c/normPD;
		PlaneDefinition[3] = d;
		//PlaneDefinition.assign(pd, pd+4);
		experimentalsetupvariables.L_c = calibrationValues.c/(calibrationValues.a+calibrationValues.b+calibrationValues.c)*calibrationValues.L_m-experimentalsetupvariables.L_s-2.0*experimentalsetupvariables.L_g-experimentalsetupvariables.L_t;
		Lengths[0] = experimentalsetupvariables.L_c;
		Lengths[2] = experimentalsetupvariables.L_t;
		//calculateNFigures(GridX, GridY, DX, DY, CC, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, inputvariables.path);
		
		auto tr5= std::chrono::high_resolution_clock::now();
		//cv::Mat nfield = CalculateN(GridX, GridY, DX, DY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, inputvariables.path);
		cv::Mat nfield = CalculateN2(GridX, GridY, calibrationValues.meanGridX, calibrationValues.meanGridY, DX, DY, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, PlaneDefinition, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, inputvariables.path);
		store_matrix(inputvariables.path,"nfield", nfield);
		auto tr6= std::chrono::high_resolution_clock::now();
		std::cout << "n Calculation took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr6-tr5).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr6-tr5).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr6-tr5).count() << " minutes"<<std::endl;
		std::cout << std::endl << "\033[1;32mN Calculation Completed\033[0m\n" << std::endl;
	}

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
	unsigned int numberofrows, numberofcols;
    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
		if (line>=skiplines)
		{
			numberofrows = static_cast<unsigned int>((*loop).size());
		}
		line++;
    }
	numberofcols = static_cast<unsigned int>(line-skiplines);
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

