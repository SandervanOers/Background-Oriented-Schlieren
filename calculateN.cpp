#include "calculateN.hpp"

/*--------------------------------------------------------------------------*/
static cv::Mat signum(cv::Mat src)
{
    cv::Mat dst = (src >= 0) & 1;
    dst.convertTo(dst,CV_64F, 2.0, -1.0);
    return dst;
}
/*--------------------------------------------------------------------------*/
static double calculateMean(const cv::Mat &Mat)
{
	Scalar meanMatrix = mean(Mat);
	return meanMatrix[0];
}/*--------------------------------------------------------------------------*/
extern double calculateLf(const double &focal_length, const double &Lm)
{
	return 1.0/(1.0/focal_length-1.0/Lm);
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateDirectionCosines(const cv::Mat &GridX, const cv::Mat &GridY, const double &L_f, const double &Distance_From_Pixels_To_Meters)
{
	cv::Mat L_F(GridX.size(), CV_64FC1, L_f);
	
	cv::Mat X = GridX-calculateMean(GridX);
	cv::Mat Y = GridY-calculateMean(GridY);
	// Scale X and Y from pixels to meters
	// Depends on Camera Properties
	X = X*Distance_From_Pixels_To_Meters;
	Y = Y*Distance_From_Pixels_To_Meters;
	
	cv::Mat normII = X.mul(X)+Y.mul(Y)+L_f*L_f;
	cv::sqrt(normII, normII);
	//std::cout << normII << std::endl<< std::endl;
	/*for (unsigned int i = 0; i < GridX.rows; i++) 
	{
		for (unsigned int j = 0; j < GridX.cols; j++)
		{
			// Direction vector I = (-X_i, -Y_j, L_f)
			// Norm of Direction Vector ||I|| = sqrt(X_i^2+Y_j^2+L_f^2);
			double normI = sqrt(X.at<double>(i,j)*X.at<double>(i,j)+Y.at<double>(i,j)*Y.at<double>(i,j)+L_f*L_f);
			
			alpha.at<double>(i,j) = -X.at<double>(i,j)/normI;
			beta.at<double>(i,j) = -Y.at<double>(i,j)/normI;
			gamma.at<double>(i,j) = L_f/normI;
			// Check:  Result should always be one
			//std::cout << alpha.at<double>(i,j)*alpha.at<double>(i,j) + beta.at<double>(i,j)*beta.at<double>(i,j)+gamma.at<double>(i,j)*gamma.at<double>(i,j) << std::endl;
		}
	}*/
	cv::Mat alpha = X.mul(-1.0/normII);
	cv::Mat beta = Y.mul(-1.0/normII);
	cv::Mat gamma = L_F.mul(1.0/normII);
	std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(alpha);
    ReturnVector.push_back(beta);
    ReturnVector.push_back(gamma);
    return ReturnVector;	
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateIntersectionPlaneLine(const std::vector<cv::Mat> &InitialPosition, const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition)
{
	// Equation Plane
	// a x + b y + c z + d = 0
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	double d = PlaneDefinition[3];
	// Equation Line
	// P = S + I L
	// P = Intersection Point (unknown)
	// S = Initial Position
	// I = Direction Cosines
	// L = Length (unknown)
	cv::Mat S_x = InitialPosition.at(0);
	cv::Mat S_y = InitialPosition.at(1);
	cv::Mat S_z = InitialPosition.at(2);
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);
	
	cv::Mat L(S_x.size(), CV_64FC1, Scalar(0));
	// First Solve for L
	//for (unsigned int i = 0; i < S_x.rows; i++)
	//{
	//	for (unsigned int j = 0; j < S_x.cols; j++)
	//	{
	//		L.at<double>(i,j) = - (d+a*S_x.at<double>(i,j)+b*S_y.at<double>(i,j)+c*S_z.at<double>(i,j))/(a*alpha.at<double>(i,j)+b*beta.at<double>(i,j)+c*gamma.at<double>(i,j));
	//	}
	//}
	L = (d+a*S_x+b*S_y+c*S_z).mul(-1.0/(a*alpha+b*beta+c*gamma));
	
	// Calculate Intersection Point
	//cv::Mat P_x(S_x.size(), CV_64FC1, Scalar(0));
	//cv::Mat P_y(S_x.size(), CV_64FC1, Scalar(0));
	//cv::Mat P_z(S_x.size(), CV_64FC1, Scalar(0));
	cv::Mat P_x = S_x+alpha.mul(L);
	cv::Mat P_y = S_y+beta.mul(L);
	cv::Mat P_z = S_z+gamma.mul(L);
	
	std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(P_x);
    ReturnVector.push_back(P_y);
    ReturnVector.push_back(P_z);
    return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLaw(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1)
{
	// Normal vector of Plane
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);
	
	//cv::Mat alphasquare(alpha.size(), CV_64FC1);
	//cv::Mat betasquare(alpha.size(), CV_64FC1);
	//cv::Mat gammasquare(alpha.size(), CV_64FC1);
	//cv::Mat total(alpha.size(), CV_64FC1);
	//cv::pow(alpha, 2, alphasquare);
	//cv::pow(beta, 2, betasquare);
	//cv::pow(gamma, 2, gammasquare);
	//std::cout << "Total = " << alphasquare + betasquare + gammasquare << std::endl << std::endl;
	
	// Angle of Incidence
	cv::Mat cosIncidence = -(alpha*a+beta*b+gamma*c);
	//std::cout << "cosIncidence: " << std::endl << cosIncidence << std::endl << std::endl;
	//std::cout << "Sign cosIncidence: " << std::endl << signum(cosIncidence) << std::endl << std::endl;
	//cv::Mat Sign(alpha.size(), CV_64FC1);
	//Sign = signum(cosIncidence);
	cv::Mat Sign = signum(cosIncidence);
	//std::cout << "Ratio = " << n_0/n_1 << std::endl<< std::endl;
	// Snell's law
	cv::Mat lhsSnell;
	cv::sqrt(1.0-cosIncidence.mul(cosIncidence),lhsSnell);
	cv::Mat sinRefracted =  n_0/n_1*lhsSnell;
	cv::Mat cosRefracted;
	cv::sqrt(1.0-sinRefracted.mul(sinRefracted),cosRefracted);
	//std::cout << "cosRefracted: " << std::endl << cosRefracted << std::endl << std::endl;
	
	//std::cout << "cosRefracted Type " << cosRefracted.type() << std::endl;
	//std::cout << "Signum Type " << Sign.type() << std::endl;
	// Refracted Ray
	cv::Mat T_x = n_0/n_1*alpha+(n_0/n_1*cosIncidence-Sign.mul(cosRefracted))*a;
	cv::Mat T_y = n_0/n_1*beta+(n_0/n_1*cosIncidence-Sign.mul(cosRefracted))*b;
	cv::Mat T_z = n_0/n_1*gamma+(n_0/n_1*cosIncidence-Sign.mul(cosRefracted))*c;
	//cv::pow(T_x, 2, alphasquare);
	//cv::pow(T_y, 2, betasquare);
	//cv::pow(T_z, 2, gammasquare);
	//std::cout << "New Total = " << alphasquare + betasquare + gammasquare << std::endl << std::endl;
	
	/*
	// New Method
	cosIncidence = (alpha*a+beta*b+gamma*c);
	cv::Mat k1(alpha.size(), CV_64FC1);
	cv::Mat k2(alpha.size(), CV_64FC1);
	cv::Mat rhs(alpha.size(), CV_64FC1);
	cv::sqrt(cosIncidence.mul(cosIncidence)-(1.0-(n_0/n_1)*(n_0/n_1)), rhs);
	k1 = - cosIncidence + rhs;
	k2 = - cosIncidence - rhs;
	cv::Mat F1x(alpha.size(), CV_64FC1);
	cv::Mat F1y(alpha.size(), CV_64FC1);
	cv::Mat F1z(alpha.size(), CV_64FC1);
	cv::Mat F2x(alpha.size(), CV_64FC1);
	cv::Mat F2y(alpha.size(), CV_64FC1);
	cv::Mat F2z(alpha.size(), CV_64FC1);
	
	F1x = alpha+k1*a;
	F1y = beta+k1*b;
	F1z = gamma+k1*c;
	F2x = alpha+k2*a;
	F2y = beta+k2*b;
	F2z = gamma+k2*c;
	
	cv::Mat Incidence1 = (alpha*F1x+beta*F1y+gamma*F1z);
	cv::Mat Incidence2 = (alpha*F2x+beta*F2y+gamma*F2z);
	
	std::cout << std::endl << "New Method " << std::endl << std::endl;;
	std::cout << "Incidence 1 = " << Incidence1 << std::endl << std::endl;;
	std::cout << "Incidence 2 = " << Incidence2 << std::endl << std::endl;;
	
	//std::cout << "F1 X Direction: " << std::endl << F1x << std::endl << std::endl;
	//std::cout << "F1 Z Direction: " << std::endl << F1z << std::endl << std::endl;	
	//std::cout << "F2 X Direction: " << std::endl << F2x << std::endl << std::endl;
	//std::cout << "F2 Z Direction: " << std::endl << F2z << std::endl << std::endl;

	cv::sqrt(F1x.mul(F1x)+F1y.mul(F1y)+F1z.mul(F1z), total);
	F1x = F1x.mul(1.0/total);
	F1y = F1y.mul(1.0/total);
	F1z = F1z.mul(1.0/total);	
	std::cout << "F1 X Direction: " << std::endl << F1x << std::endl << std::endl;
	std::cout << "F1 Z Direction: " << std::endl << F1z << std::endl << std::endl;

	cv::pow(F1x, 2, alphasquare);
	cv::pow(F1y, 2, betasquare);
	cv::pow(F1z, 2, gammasquare);	
	std::cout << "F1 Total = " << alphasquare + betasquare + gammasquare << std::endl << std::endl;
	cv::sqrt(F2x.mul(F2x)+F2y.mul(F2y)+F2z.mul(F2z), total);
	F2x = F2x.mul(1.0/total);
	F2y = F2y.mul(1.0/total);
	F2z = F2z.mul(1.0/total);
	std::cout << "F2 X Direction: " << std::endl << F2x << std::endl << std::endl;
	std::cout << "F2 Z Direction: " << std::endl << F2z << std::endl << std::endl;
	cv::pow(F2x, 2, alphasquare);
	cv::pow(F2y, 2, betasquare);
	cv::pow(F2z, 2, gammasquare);	
	std::cout << "F2 Total = " << alphasquare + betasquare + gammasquare << std::endl << std::endl;
	 * */
	std::vector<cv::Mat> ReturnVector;
	ReturnVector.push_back(T_x);
	ReturnVector.push_back(T_y);
	ReturnVector.push_back(T_z);
	return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLawMat(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const cv::Mat &n_0, const cv::Mat &n_1)
{
	// Normal vector of Plane
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);
	
	//cv::Mat alphasquare(alpha.size(), CV_64FC1);
	//cv::Mat betasquare(alpha.size(), CV_64FC1);
	//cv::Mat gammasquare(alpha.size(), CV_64FC1);
	//cv::pow(alpha, 2, alphasquare);
	//cv::pow(beta, 2, betasquare);
	//cv::pow(gamma, 2, gammasquare);
	//std::cout << "Total = " << alphasquare + betasquare + gammasquare << std::endl << std::endl;
	// Angle of Incidence
	cv::Mat cosIncidence = -(alpha*a+beta*b+gamma*c);
	cv::Mat Sign = signum(cosIncidence);
	//std::cout << "gamma = " << gamma << std::endl << std::endl;
	//std::cout << "cosIncidence = " << cosIncidence << std::endl<< std::endl;
	cv::Mat Ratio = n_0.mul(1.0/n_1);
	//std::cout << "n_0 = " << n_0 << std::endl << std::endl;
	//std::cout << "n_1 = " << n_1 << std::endl << std::endl;
	//std::cout << "Ratio = " << Ratio << std::endl<< std::endl;
	// Snell's law
	cv::Mat lhsSnell;
	cv::sqrt(1.0-cosIncidence.mul(cosIncidence),lhsSnell);
	//std::cout << lhsSnell << std::endl<< std::endl;
	cv::Mat sinRefracted = Ratio.mul(lhsSnell);
	//std::cout << "sinRefracted = " << sinRefracted << std::endl<< std::endl;
	cv::Mat cosRefracted;
	cv::sqrt(1.0-sinRefracted.mul(sinRefracted),cosRefracted);
	//std::cout << "cosRefracted = " << cosRefracted << std::endl<< std::endl;
	
	// Refracted Ray
	cv::Mat T_x = Ratio.mul(alpha)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*a;
	cv::Mat T_y = Ratio.mul(beta)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*b;
	cv::Mat T_z = Ratio.mul(gamma)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*c;
	//cv::pow(T_x, 2, alphasquare);
	//cv::pow(T_y, 2, betasquare);
	//cv::pow(T_z, 2, gammasquare);
	//std::cout << "New Total = " << alphasquare + betasquare + gammasquare << std::endl << std::endl;
	
	std::vector<cv::Mat> ReturnVector;
	ReturnVector.push_back(T_x);
	ReturnVector.push_back(T_y);
	ReturnVector.push_back(T_z);
	return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern PositionDirection calculateIntersectionRungeKutta(const std::vector<cv::Mat> &InitialPosition, const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const cv::Mat &n_field, const double &n_0, const double &n_2, const unsigned int &SplineDegree, const double &L_t, const unsigned int &Number_Of_Steps)
{
	cv::Mat S_x = InitialPosition.at(0);
	cv::Mat S_y = InitialPosition.at(1);
	cv::Mat S_z = InitialPosition.at(2);
	unsigned int rows = S_x.rows;
	unsigned int cols = S_x.cols;
	
	cv::Mat n_0_mat(S_x.size(), CV_64FC1, n_0);
	cv::Mat n_1_mat(S_x.size(), CV_64FC1);
	// Step 1: Determine Coordinate Transformation
	cv::Mat Q = calculateTransformationMatrix(PlaneDefinition);
	cv::Mat Qt;
	cv:transpose(Q, Qt);
	//std::cout << std::endl << "Q = " << Q << std::endl << std::endl; 
	//std::cout << std::endl << "Qt = " << Qt << std::endl << std::endl; 
	// Step 2: Transform Initial Positions to new coordinates u,v 
	// Old unit vectors
	cv::Vec3d unit_x (1,0,0);
	cv::Vec3d unit_y (0,1,0);
	cv::Vec3d unit_z (0,0,1);
	// New Unit Vectors
	cv::Mat unit_u = Q*cv::Mat(unit_x,false);
	cv::Mat unit_v = Q*cv::Mat(unit_y,false);
	cv::Mat unit_w = Q*cv::Mat(unit_z,false);
	std::cout << std::endl;
	//std::cout << "unit_u = " << unit_u << std::endl;
	//std::cout << "unit_v = " << unit_v << std::endl;
	//std::cout << "unit_w = " << unit_w << std::endl;
	std::cout << "n = (" << PlaneDefinition[0] << ", " << PlaneDefinition[1] << ", " << PlaneDefinition[2] << ")"<< std::endl;
	//std::cout << std::endl << "S_x = " << std::endl << S_x << std::endl << std::endl;
	//std::cout << std::endl << "S_y = " << std::endl << S_y << std::endl << std::endl;
	//std::cout << std::endl << "S_z = " << std::endl << S_z << std::endl << std::endl;
	// Apply Coordinate Transformation to each Position
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			cv::Vec3d position(S_x.at<double>(i,j), S_y.at<double>(i,j), S_z.at<double>(i,j));
			double S_u = position.dot(unit_u);
			double S_v = position.dot(unit_v);
			double S_w = position.dot(unit_w);
			S_x.at<double>(i,j) = S_u;
			S_y.at<double>(i,j) = S_v;
			S_z.at<double>(i,j) = S_w; //Should all be equal
		}
	}
	//std::cout << std::endl << "S_x = " << std::endl << S_x << std::endl << std::endl;
	//std::cout << std::endl << "S_y = " << std::endl << S_y << std::endl << std::endl;
	//std::cout << std::endl << "S_z = " << std::endl << S_z << std::endl << std::endl;	
	// Step 3: Initiate n_1 in new coordinates u, v
	
	// First calculate B-splines coefficients for the n_field
    cv::Mat prova_img2 = n_field.clone();
	prova_img2.convertTo(prova_img2, CV_32F);
    float *fptr_img2 = prova_img2.ptr<float>(0);
    SamplesToCoefficients(fptr_img2, cols, rows, SplineDegree);
	// Mapping from indices for n_field to coordinates in S_x and S_y:
	// x = min(S_x) + i *(max(S_x)-min(S_x))/length(S_x);
	// y = min(S_y) + j *(max(S_y)-min(S_y))/length(S_y);
	// i = length(S_x) * (x - min(S_x))/(max(S_x)-min(S_x));
	// j = length(S_y) * (y - min(S_y))/(max(S_y)-min(S_y));
	double minSx, maxSx, minSy, maxSy;
	cv::minMaxLoc(S_x, &minSx, &maxSx);	
	cv::minMaxLoc(S_y, &minSy, &maxSy);	
	//std::cout << "minSx = " << minSx << ", maxSx = " << maxSx << std::endl;
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			double X = S_x.at<double>(i,j);
			double Y = S_y.at<double>(i,j);
			double I = cols*(X-minSx)/(maxSx-minSx);
			double J = rows*(Y-minSy)/(maxSy-minSy);
			n_1_mat.at<double>(j,i) = InterpolatedValue(fptr_img2,  cols, rows, I, J, SplineDegree);
		}
	}
	//std::cout << std::endl << "n_field = " << n_field << std::endl << std::endl;
	//std::cout << std::endl << "n_1_mat = " << n_1_mat << std::endl << std::endl;
	//std::cout << std::endl << "n_field- n_1_mat = " << n_field-n_1_mat << std::endl << std::endl;
	// Step 4: Apply Snell's Law in old coordinates (x,y,z) using the new n_1(u,v)
	std::vector<cv::Mat> Refracted = SnellsLawMat(DirectionCosines, PlaneDefinition, n_0_mat, n_1_mat);
	// Step 5: Transform Refracted Rays to new coordinates u, v, w
	cv::Mat T_x = Refracted.at(0);
	cv::Mat T_y = Refracted.at(1);
	cv::Mat T_z = Refracted.at(2);
	std::cout << std::endl << "T_x = " << std::endl << T_x << std::endl << std::endl;
	//std::cout << std::endl << "T_y = " << std::endl << T_y << std::endl << std::endl;
	std::cout << std::endl << "T_z = " << std::endl << T_z << std::endl << std::endl;
	// Apply Coordinate Transformation to each Ray
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			cv::Vec3d ray(T_x.at<double>(i,j), T_y.at<double>(i,j), T_z.at<double>(i,j));
			cv::Mat ResultRay = Q*cv::Mat(ray, false);
			T_x.at<double>(i,j) = ResultRay.at<double>(0);
			T_y.at<double>(i,j) = ResultRay.at<double>(1);
			T_z.at<double>(i,j) = ResultRay.at<double>(2);
		}
	}
	// Step 6: Runge-Kutta time
	double h = L_t/Number_Of_Steps;
	// Solves ray equation: d/ds (n(r) dr/ds) = grad n
	// T = n dr/ds
	// d/ds T = grad n
	// d/ds r = T/n
	// Uses 4th order Runge Kutta
	
	// Step 7: Obtain grad n = ( dn/du, dn/dv, dn/dw)
	 // Assumption: dn / dw = 0

	//cv::Mat S_z0 = S_z.clone();
	//double Sz0 = calculateMean(S_z);
	double Sz0 = S_z.at<double>(0.0);
	std::cout << "Sz0 = " << Sz0 << std::endl; 
	//std::cout << std::endl << "S_x = " << std::endl << S_x << std::endl << std::endl;
	std::cout << std::endl << "S_z = " << std::endl << S_z << std::endl << std::endl;
	std::cout << std::endl << "T_x = " << std::endl << T_x << std::endl << std::endl;
	//std::cout << std::endl << "T_y = " << std::endl << T_y << std::endl << std::endl;
	std::cout << std::endl << "T_z = " << std::endl << T_z << std::endl << std::endl;
	double hold = h;
	cv::Mat done(S_x.size(), CV_8U, Scalar(0));
	cv::Scalar total = done.total();
	
	std::vector<int> IndexI(S_x.rows);
	std::iota(IndexI.begin(), IndexI.end(), 0);
	std::vector<int> IndexJ(S_x.cols);
	std::iota(IndexJ.begin(), IndexJ.end(), 0);
	std::vector<int> Index(total[0]);
	std::iota(Index.begin(), Index.end(), 0);
	
	std::cout << L_t << std::endl;
	while (cv::sum(done)[0] < total[0])
	{
		cv::Mat k1_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k1_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k1_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l1_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l1_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l1_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k2_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k2_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k2_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l2_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l2_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l2_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k3_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k3_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k3_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l3_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l3_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l3_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k4_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k4_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat k4_w(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l4_u(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l4_v(S_x.size(), CV_64FC1, Scalar(0));
		cv::Mat l4_w(S_x.size(), CV_64FC1, Scalar(0));
		for (auto const& singleindex: Index)
		{
			unsigned int i = singleindex/cols;
			unsigned int j = singleindex%cols;
			//if (S_z.at<double>(i,j) - S_z0.at<double>(i,j) >= L_t)
			if (S_z.at<double>(i,j) - Sz0 >= L_t)
			{
				done.at<uchar>(i,j) = 1;
				Index.erase( std::remove( Index.begin(), Index.end(), singleindex ), Index.end() ); 
			}
			else
			{
				//h = hold;
				//if (S_z.at<double>(i,j) - S_z0.at<double>(i,j) + h > L_t)
				//{
				//	h = h/1.0;
				//}
				// Step 1
				double X = S_x.at<double>(i,j);
				double Y = S_y.at<double>(i,j);
				double I = cols*(X-minSx)/(maxSx-minSx);
				double J = rows*(Y-minSy)/(maxSy-minSy);
				k1_u.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 0);
				k1_v.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 1);
				k1_w.at<double>(i,j) = h * 0;
				
				l1_u.at<double>(i,j) = h * T_x.at<double>(i,j) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l1_v.at<double>(i,j) = h * T_y.at<double>(i,j) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l1_w.at<double>(i,j) = h * T_z.at<double>(i,j) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				
				// Step 2
				X = S_x.at<double>(i,j) + l1_u.at<double>(i,j);
				Y = S_y.at<double>(i,j) + l1_v.at<double>(i,j);
				I = cols*(X-minSx)/(maxSx-minSx);
				J = rows*(Y-minSy)/(maxSy-minSy);
				k2_u.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 0);
				k2_v.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 1);
				k2_w.at<double>(i,j) = h * 0;
				
				l2_u.at<double>(i,j) = h * (T_x.at<double>(i,j) + 0.5*k1_u.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l2_v.at<double>(i,j) = h * (T_y.at<double>(i,j) + 0.5*k1_v.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l2_w.at<double>(i,j) = h * (T_z.at<double>(i,j) + 0.5*k1_w.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				
				// Step 3
				X = S_x.at<double>(i,j) + l2_u.at<double>(i,j);
				Y = S_y.at<double>(i,j) + l2_v.at<double>(i,j);
				I = cols*(X-minSx)/(maxSx-minSx);
				J = rows*(Y-minSy)/(maxSy-minSy);
				k3_u.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 0);
				k3_v.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 1);
				k3_w.at<double>(i,j) = h * 0;
				
				l3_u.at<double>(i,j) = h * (T_x.at<double>(i,j) + 0.5*k2_u.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l3_v.at<double>(i,j) = h * (T_y.at<double>(i,j) + 0.5*k2_v.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l3_w.at<double>(i,j) = h * (T_z.at<double>(i,j) + 0.5*k2_w.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				
				// Step 4
				X = S_x.at<double>(i,j) + l3_u.at<double>(i,j);
				Y = S_y.at<double>(i,j) + l3_v.at<double>(i,j);
				I = cols*(X-minSx)/(maxSx-minSx);
				J = rows*(Y-minSy)/(maxSy-minSy);
				k4_u.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 0);
				k4_v.at<double>(i,j) = h * getDerivativeValue(fptr_img2, cols, rows, I, J, SplineDegree, 1);
				k4_w.at<double>(i,j) = h * 0;
				
				l4_u.at<double>(i,j) = h * (T_x.at<double>(i,j) + 0.5*k3_u.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l4_v.at<double>(i,j) = h * (T_y.at<double>(i,j) + 0.5*k3_v.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
				l4_w.at<double>(i,j) = h * (T_z.at<double>(i,j) + 0.5*k3_w.at<double>(i,j)) / InterpolatedValue(fptr_img2, cols, rows, I, J, SplineDegree);
			}
		}
		
		cv::Mat DTu = 1.0/6.0*(k1_u+2.0*k2_u+2.0*k3_u+k4_u);
		cv::Mat DTv = 1.0/6.0*(k1_v+2.0*k2_v+2.0*k3_v+k4_v);
		cv::Mat DTw = 1.0/6.0*(k1_w+2.0*k2_w+2.0*k3_w+k4_w);
		
		cv::Mat DSu = 1.0/6.0*(l1_u+2.0*l2_u+2.0*l3_u+l4_u);
		cv::Mat DSv = 1.0/6.0*(l1_v+2.0*l2_v+2.0*l3_v+l4_v);
		cv::Mat DSw = 1.0/6.0*(l1_w+2.0*l2_w+2.0*l3_w+l4_w);
		
		S_x = S_x + DSu;
		S_y = S_y + DSv;
		S_z = S_z + DSw;
		
		T_x = T_x + DTu;
		T_y = T_y + DTu;
		T_z = T_z + DTu;
	}
	//std::cout << std::endl << "S_x = " << std::endl << S_x << std::endl << std::endl;
	//std::cout << std::endl << "S_z = " << std::endl << S_z << std::endl << std::endl;
	//std::cout << std::endl << "S_z-S_z0_L_t = " << std::endl << S_z-S_z0-L_t << std::endl << std::endl;
	//std::cout << "L_t = " << L_t << std::endl;
	//std::cout << std::endl << "T_z = " << std::endl << T_z << std::endl << std::endl;
	 
	 // Step 7: Transform incomping rays back to coordinates x y z
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			cv::Vec3d ray(T_x.at<double>(i,j), T_y.at<double>(i,j), T_z.at<double>(i,j));
			cv::Mat ResultRay = Qt*cv::Mat(ray, false);
			T_x.at<double>(i,j) = ResultRay.at<double>(0);
			T_y.at<double>(i,j) = ResultRay.at<double>(1);
			T_z.at<double>(i,j) = ResultRay.at<double>(2);
		}
	}
	//std::cout << std::endl << "T_z = " << std::endl << T_z << std::endl << std::endl;
	
	// Step 8: Apply Snell's Law in old coordinates (x,y,z) using the new n_1(u,v)
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			double X = S_x.at<double>(i,j);
			double Y = S_y.at<double>(i,j);
			double I = cols*(X-minSx)/(maxSx-minSx);
			double J = rows*(Y-minSy)/(maxSy-minSy);
			n_1_mat.at<double>(j,i) = InterpolatedValue(fptr_img2,  cols, rows, I, J, SplineDegree);
		}
	}
	cv::Mat n_2_mat(S_x.size(), CV_64FC1, n_2);
	std::vector<cv::Mat> DirectionCosinesExit;
	DirectionCosinesExit.push_back(T_x);
	DirectionCosinesExit.push_back(T_y);
	DirectionCosinesExit.push_back(T_z);
	std::vector<cv::Mat> RefractedExit = SnellsLawMat(DirectionCosinesExit, PlaneDefinition, n_1_mat, n_2_mat);
	
	// Step 9: Transform Final Positions to new old coordinates x, y, z
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
		{
			//cv::Vec3d position(S_x.at<double>(i,j), S_y.at<double>(i,j), S_z.at<double>(i,j));
			cv::Vec3d position(S_x.at<double>(i,j), S_y.at<double>(i,j), Sz0+L_t);
			double S_u = position.dot(unit_x);
			double S_v = position.dot(unit_y);
			double S_w = position.dot(unit_z);
			S_x.at<double>(i,j) = S_u;
			S_y.at<double>(i,j) = S_v;
			S_z.at<double>(i,j) = S_w; //Should not all be equal
		}
	}
	std::vector<cv::Mat> ReturnPosition;
	ReturnPosition.push_back(S_x);
	ReturnPosition.push_back(S_y);
	ReturnPosition.push_back(S_z);
	PositionDirection Return(ReturnPosition, RefractedExit);
	return Return;
}
/*--------------------------------------------------------------------------*/
extern PositionDirection calculateIntersectionConstantRefraction(const PositionDirection &InitialPositionDirection, const std::vector<double> &PlaneDefinition, const double &n_initial, const double &n_final)
{
	std::vector<cv::Mat> InitialPosition = InitialPositionDirection.Position;
	std::vector<cv::Mat> InitialDirection = InitialPositionDirection.Direction;
	std::vector<cv::Mat> Intersection = calculateIntersectionPlaneLine(InitialPosition, InitialDirection, PlaneDefinition);
	std::vector<cv::Mat> Refracted = SnellsLaw(InitialDirection, PlaneDefinition, n_initial, n_final);
	PositionDirection Return(Intersection, Refracted);
	return Return;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> ForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const double &Lm, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const cv::Mat &n_field, const unsigned int &SplineDegree, const unsigned int &Number_Of_Steps)
{
	double L_c = Lengths[0];
	double L_g = Lengths[1];
	double L_t = Lengths[2];
	double L_s = Lengths[4];
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	double d = PlaneDefinition[3];
	double d5 = d+(a+b+c)*L_s;
	double d4 = d+(a+b+c)*(L_s+L_g);
	double d3 = d+(a+b+c)*(L_s+L_g+L_t);
	double d2 = d+(a+b+c)*(L_s+2.0*L_g+L_t);
	double d1 = d+(a+b+c)*(L_s+2.0*L_g+L_t+L_c);
	//std::vector<double> PlaneDefinition1{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d1}; // Should pass through Origin (0,0,0)
	std::vector<double> PlaneDefinition2{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d2};
	std::vector<double> PlaneDefinition3{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d3};
	std::vector<double> PlaneDefinition4{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d4};
	std::vector<double> PlaneDefinition5{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d5};
	
	std::cout << "Plane 6 Definition: " << PlaneDefinition[0] << "x + " << PlaneDefinition[1] << "y + " << PlaneDefinition[2] << "z + " << PlaneDefinition[3] << " = 0" << ": (x,y)=(0,0) => z_6 = " << -PlaneDefinition[3]/PlaneDefinition[2] << std::endl << std::endl;
	std::cout << "Plane 5 Definition: " << PlaneDefinition[0] << "x + " << PlaneDefinition[1] << "y + " << PlaneDefinition[2] << "z + " << d5 << " = 0 " << ": (x,y)=(0,0) => z_5 = " << -d5/PlaneDefinition[2] << std::endl << std::endl;
	std::cout << "Plane 4 Definition: " << PlaneDefinition[0] << "x + " << PlaneDefinition[1] << "y + " << PlaneDefinition[2] << "z + " << d4 << " = 0 " << ": (x,y)=(0,0) => z_4 = " << -d4/PlaneDefinition[2] << std::endl << std::endl;
	std::cout << "Plane 3 Definition: " << PlaneDefinition[0] << "x + " << PlaneDefinition[1] << "y + " << PlaneDefinition[2] << "z + " << d3 << " = 0 " << ": (x,y)=(0,0) => z_3 = " << -d3/PlaneDefinition[2] << std::endl << std::endl;
	std::cout << "Plane 2 Definition: " << PlaneDefinition[0] << "x + " << PlaneDefinition[1] << "y + " << PlaneDefinition[2] << "z + " << d2 << " = 0 " << ": (x,y)=(0,0) => z_2 = " << -d2/PlaneDefinition[2] << std::endl << std::endl;
	std::cout << "Plane 1 Definition: " << PlaneDefinition[0] << "x + " << PlaneDefinition[1] << "y + " << PlaneDefinition[2] << "z + " << d1 << " = 0 " << ": (x,y)=(0,0) => z_1 = " << -d1/PlaneDefinition[2] << std::endl << std::endl;
	
	/* // Test 1
	std::vector<cv::Mat> InitialPosition;
	cv::Mat S_x(1,1, CV_64FC1, Scalar(0));
	cv::Mat S_y(1,1, CV_64FC1, Scalar(0));
	cv::Mat S_z(1,1, CV_64FC1, Scalar(0));
	InitialPosition.push_back(S_x);
	InitialPosition.push_back(S_y);
	InitialPosition.push_back(S_z);
	cv::Mat Wx(1,1, CV_64FC1, Scalar(.707107));
	cv::Mat Wy(1,1, CV_64FC1, Scalar(0));
	cv::Mat Wz(1,1, CV_64FC1, Scalar(-.707107));
	std::vector<cv::Mat> InitialDirections;
	InitialDirections.push_back(Wx);
	InitialDirections.push_back(Wy);
	InitialDirections.push_back(Wz);
	PositionDirection Plane1(InitialPosition, InitialDirections);
	std::vector<double> PlaneDefinitionW{0, 0, 1};
	double n0W = 0.9; double n1W = 1.0;
	PositionDirection PlaneW = calculateIntersectionConstantRefraction(Plane1,  PlaneDefinitionW, n0W, n1W);
	*/ // Test 2 
	/*
	std::vector<cv::Mat> InitialPosition;	
	cv::Mat S_x(1,1, CV_64FC1, Scalar(0));
	cv::Mat S_y(1,1, CV_64FC1, Scalar(0));
	cv::Mat S_z(1,1, CV_64FC1, Scalar(0));
	InitialPosition.push_back(S_x);
	InitialPosition.push_back(S_y);
	InitialPosition.push_back(S_z);
	cv::Mat Wx(1,1, CV_64FC1, Scalar(4));
	cv::Mat Wy(1,1, CV_64FC1, Scalar(1));
	cv::Mat Wz(1,1, CV_64FC1, Scalar(1));
	std::vector<cv::Mat> InitialDirections;
	InitialDirections.push_back(Wx);
	InitialDirections.push_back(Wy);
	InitialDirections.push_back(Wz);
	PositionDirection Plane1(InitialPosition, InitialDirections);
	std::vector<double> PlaneDefinitionW{0, -2, -1};
	double n0W = 1; double n1W = 1.5;
	PositionDirection PlaneW = calculateIntersectionConstantRefraction(Plane1,  PlaneDefinitionW, n0W, n1W);	
	*/
	// Full Sim
	// Unknown n
	double L_f = calculateLf(focal_length, Lm);
	std::vector<cv::Mat> InitialDirection = calculateDirectionCosines(GridX, GridY, L_f, Distance_From_Pixels_To_Meters);
	std::vector<cv::Mat> InitialPosition;
	cv::Mat S_x(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat S_y(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat S_z(GridX.size(), CV_64FC1, Scalar(0));
	InitialPosition.push_back(S_x);
	InitialPosition.push_back(S_y);
	InitialPosition.push_back(S_z);
	PositionDirection Plane1(InitialPosition, InitialDirection);
	std::cout << std::endl << "--- Plane 2 --- " << std::endl << std::endl;
	PositionDirection Plane2 = calculateIntersectionConstantRefraction(Plane1, PlaneDefinition2, n_0, n_1);
	std::cout << std::endl << "Plane 2 Position X = " << std::endl << Plane2.Position.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 2 Position Y = " << std::endl << Plane2.Position.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 2 Position Z = " << std::endl << Plane2.Position.at(2) << std::endl << std::endl;
	std::cout << std::endl << "Plane 2 Direction X = " << std::endl << Plane2.Direction.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 2 Direction Y = " << std::endl << Plane2.Direction.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 2 Direction Z = " << std::endl << Plane2.Direction.at(2) << std::endl << std::endl;
	std::cout << std::endl << "--- Plane 3 --- " << std::endl << std::endl;
	std::vector<cv::Mat> PositionPlane3 = calculateIntersectionPlaneLine(Plane2.Position, Plane2.Direction, PlaneDefinition3);
	std::cout << std::endl << "Plane 3 Position X = " << std::endl << PositionPlane3.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 3 Position Y = " << std::endl << PositionPlane3.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 3 Position Z = " << std::endl << PositionPlane3.at(2) << std::endl << std::endl;
	std::cout << std::endl << "Plane 3 Direction X = " << std::endl << Plane2.Direction.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 3 Direction Y = " << std::endl << Plane2.Direction.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 3 Direction Z = " << std::endl << Plane2.Direction.at(2) << std::endl << std::endl;
	std::cout << std::endl << "--- Plane 4 --- " << std::endl << std::endl;
	PositionDirection Plane4 = calculateIntersectionRungeKutta(PositionPlane3, Plane2.Direction, PlaneDefinition4, n_field, n_1, n_1, SplineDegree,  L_t, Number_Of_Steps);
	std::cout << std::endl << "Plane 4 Position X = " << std::endl << Plane4.Position.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 4 Position Y = " << std::endl << Plane4.Position.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 4 Position Z = " << std::endl << Plane4.Position.at(2) << std::endl << std::endl;
	std::cout << std::endl << "Plane 4 Direction X = " << std::endl << Plane4.Direction.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 4 Direction Y = " << std::endl << Plane4.Direction.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 4 Direction Z = " << std::endl << Plane4.Direction.at(2) << std::endl << std::endl;
	std::cout << std::endl << "--- Plane 5 --- " << std::endl << std::endl;
	PositionDirection Plane5 = calculateIntersectionConstantRefraction(Plane4, PlaneDefinition5, n_1, n_0);
	std::cout << std::endl << "Plane 5 Position X = " << std::endl << Plane5.Position.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 5 Position Y = " << std::endl << Plane5.Position.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 5 Position Z = " << std::endl << Plane5.Position.at(2) << std::endl << std::endl;
	std::cout << std::endl << "Plane 5 Direction X = " << std::endl << Plane5.Direction.at(0) << std::endl << std::endl;
	//std::cout << std::endl << "Plane 5 Direction Y = " << std::endl << Plane5.Direction.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Plane 5 Direction Z = " << std::endl << Plane5.Direction.at(2) << std::endl << std::endl;
	std::cout << std::endl << "--- Plane 6 --- " << std::endl << std::endl;
	std::vector<cv::Mat> PositionPlane6_unknown = calculateIntersectionPlaneLine(Plane5.Position, Plane5.Direction, PlaneDefinition);
	
	// Known n
	double n_2 = n_1;
	std::cout << std::endl << "Initial Direction X = " << std::endl << InitialDirection.at(0) << std::endl << std::endl;
	std::cout << std::endl << "Initial Direction Y = " << std::endl << InitialDirection.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Initial Direction Z = " << std::endl << InitialDirection.at(2) << std::endl << std::endl;
	InitialDirection = calculateDirectionCosines(GridX+Dx, GridY+Dy, L_f, Distance_From_Pixels_To_Meters);
	std::cout << std::endl << "Initial Direction X = " << std::endl << InitialDirection.at(0) << std::endl << std::endl;
	std::cout << std::endl << "Initial Direction Y = " << std::endl << InitialDirection.at(1) << std::endl << std::endl;
	std::cout << std::endl << "Initial Direction Z = " << std::endl << InitialDirection.at(2) << std::endl << std::endl;
	Plane2 = calculateIntersectionConstantRefraction(Plane1, PlaneDefinition2, n_0, n_1);
	PositionDirection Plane3 = calculateIntersectionConstantRefraction(Plane2, PlaneDefinition3, n_1, n_2);
	Plane4 = calculateIntersectionConstantRefraction(Plane3, PlaneDefinition3, n_2, n_1);
	Plane5 = calculateIntersectionConstantRefraction(Plane4, PlaneDefinition5, n_1, n_0);
	std::vector<cv::Mat> PositionPlane6_known = calculateIntersectionPlaneLine(Plane5.Position, Plane5.Direction, PlaneDefinition);
	
	std::vector<cv::Mat> Displacement;
	cv::Mat DX = PositionPlane6_known[0] - PositionPlane6_unknown[0];
	cv::Mat DY = PositionPlane6_known[1] - PositionPlane6_unknown[1];
	cv::Mat DZ = PositionPlane6_known[2] - PositionPlane6_unknown[2];
	Displacement.push_back(DX);
	Displacement.push_back(DY);
	Displacement.push_back(DZ);
	return Displacement;
	
}
/*--------------------------------------------------------------------------*/
static double getDerivativeValue(float *fptr_img1, const unsigned int &cols, const unsigned int &rows, const double &x, const double &y, const unsigned int &SplineDegree, const unsigned int &direction)
{
	return (double)InterpolatedValueDerivative(fptr_img1, cols, rows, x+0.5*(1.0-(double)direction), y+0.5*(double)direction, SplineDegree-1*(1-direction), SplineDegree-1*direction)
                -(double)InterpolatedValueDerivative(fptr_img1, cols, rows, x-0.5*(1.0-(double)direction), y-0.5*(double)direction, SplineDegree-1*(1-direction), SplineDegree-1*direction);
}
/*--------------------------------------------------------------------------*/
static cv::Mat calculateTransformationMatrix(const std::vector<double> &PlaneDefinition)
{
	double a[] = {PlaneDefinition[0], PlaneDefinition[1], PlaneDefinition[2]};
    double b[] = {0, 0, 1};
	
	double costheta = std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.0);
	double theta = std::acos(costheta);
	
	double normRotationVector = 1;
	if (PlaneDefinition[0] != 0 || PlaneDefinition[1] !=0)
		normRotationVector = std::sqrt(PlaneDefinition[1]*PlaneDefinition[1]+PlaneDefinition[0]*PlaneDefinition[0]);
	
	double rotationVector[] = {-PlaneDefinition[1]/normRotationVector, PlaneDefinition[0]/normRotationVector, 0};
	
	double q0 = std::cos(theta/2.0);
	double q1 = std::sin(theta/2.0)*rotationVector[0];
	double q2 = std::sin(theta/2.0)*rotationVector[1];
	double q3 = std::sin(theta/2.0)*rotationVector[2];
	
	cv::Mat Q(3, 3, CV_64FC1);
	Q.at<double>(0,0) = q0*q0+q1*q1-q2*q2-q3*q3;
	Q.at<double>(0,1) = 2*(q1*q2-q0*q3);
	Q.at<double>(0,2) = 2*(q1*q3+q0*q2);
	
	Q.at<double>(1,0) = 2*(q2*q1+q0*q3);
	Q.at<double>(1,1) = q0*q0-q1*q1+q2*q2-q3*q3;
	Q.at<double>(1,2) = 2*(q2*q3-q0*q1);
	
	Q.at<double>(2,0) = 2*(q3*q1-q0*q2);
	Q.at<double>(2,1) = 2*(q3*q2+q0*q1);
	Q.at<double>(2,2) = q0*q0-q1*q1-q2*q2+q3*q3;
	
	return Q;
	
}
/*--------------------------------------------------------------------------*/ 
namespace {
      // define y(x) = Poly(a, x) in the empty namespace
      template <class Type>
      Type Poly(const std::vector<double> &a, const Type &x)
      {     size_t k  = a.size();
            Type y   = 0.;  // initialize summation
            Type x_i = 1.;  // initialize x^i
            size_t i;
            for(i = 0; i < k; i++)
            {     y   += a[i] * x_i;  // y   = y + a_i * x^i
                  x_i *= x;           // x_i = x_i * x
            }
            return y;
      }
}
/*--------------------------------------------------------------------------*/ 
// main program
extern int poly_test(void)
{   
	std::cout << "test poly " << std::endl;  
	using CppAD::AD;           // use AD as abbreviation for CppAD::AD
      size_t i;                  // a temporary index

      // vector of polynomial coefficients
      size_t k = 5;              // number of polynomial coefficients
      std::vector<double> a(k);       // vector of polynomial coefficients
      for(i = 0; i < k; i++)
            a[i] = 1.;           // value of polynomial coefficients
	std::cout << "test poly 1" << std::endl;  
	
      // domain space vector
      size_t n = 1;              // number of domain space variables
      std::vector< AD<double> > X(1); // vector of domain space variables
      X[0] = 3.0;                 // value corresponding to operation sequence
	  std::cout << "X[0] = " << X[0] << std::endl;
	std::cout << "test poly 3" << std::endl;  

      // declare independent variables and start recording operation sequence
      CppAD::Independent(X);
	std::cout << "test poly 3.1" << std::endl;  

      // range space vector
      size_t m = 1;              // number of ranges space variables
      std::vector< AD<double> > Y(m); // vector of ranges space variables
	std::cout << "test poly 3.5" << std::endl;  
      Y[0] = Poly(a, X[0]);      // value during recording of operations
	std::cout << "test poly 4" << std::endl;  

      // store operation sequence in f: X -> Y and stop recording
      CppAD::ADFun<double> f(X, Y);
	std::cout << "test poly 5" << std::endl;  

      // compute derivative using operation sequence stored in f
      std::vector<double> jac(m * n); // Jacobian of f (m by n matrix)
      std::vector<double> x(n);       // domain space vector
      x[0] = 3.;                 // argument value for derivative
      jac  = f.Jacobian(x);      // Jacobian for operation sequence

      // print the results
      std::cout << "f'(3) computed by CppAD = " << jac[0] << std::endl;

      // check if the derivative is correct
      int error_code;
      if( jac[0] == 142. )
            error_code = 0;      // return code for correct case
      else  error_code = 1;      // return code for incorrect case

      return error_code;
	  return 0;
}
/*--------------------------------------------------------------------------*/ 