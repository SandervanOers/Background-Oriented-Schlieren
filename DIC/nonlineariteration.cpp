# include "nonlineariteration.hpp"
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_rigid_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_affine_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_irregular_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_quadratic_LM_optim(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold);
/*--------------------------------------------------------------------------*/
double Correlation_Coefficient_ZNSSD(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength);
/*--------------------------------------------------------------------------*/
double Correlation_Coefficient_ZNSSD_optim(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &um_g_minus_gm_squared);
/*--------------------------------------------------------------------------*/
double Correlation_Coefficient_ZNSSD_optim_short(const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &sum_g_minus_gm_squared, const std::vector<double> &f_values, const std::vector<double> &g_values);
/*--------------------------------------------------------------------------*/
void calculate_Hessian_Jacobian(cv::Mat &Hessian, cv::Mat &Jacobian, const cv::Mat &img, float *fptr_img1, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &ShapeFunction, const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &sum_g_minus_gm_squared, const double &lambda);
/*--------------------------------------------------------------------------*/
static double getDerivativeValue(float *fptr_img1, const unsigned int &cols, const unsigned int &rows, const double &x, const double &y, const unsigned int &SplineDegree, const unsigned int &direction);
/*--------------------------------------------------------------------------*/
extern std::vector<double> iteration(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction)
{
	
    switch(ShapeFunction)
    {
    case 0:
        return iteration_rigid_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    case 1:
        return iteration_affine_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
		break;
    case 2:
        return iteration_irregular_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    case 3:
        return iteration_quadratic_LM_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    default:
        std::cout << "Shape Function Incorrect" << std::endl;
        break;
    }    
	/*std::vector<double> Solution;
	switch(ShapeFunction)
    {
    case 0:
        Solution = iteration_rigid_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    case 1:
        Solution =  iteration_affine_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
		break;
    case 2:
        Solution =  iteration_irregular_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    case 3:
        Solution =  iteration_quadratic_LM_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    default:
        std::cout << "Shape Function Incorrect" << std::endl;
        break;
    }
	Solution.push_back(i);
	Solution.push_back(j);
	return Solution;*/
}
/*--------------------------------------------------------------------------*/
void iteration_thread(std::vector<std::vector<double>> &SolutionNeighbours, const unsigned int k, const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction)
{
	std::vector<double> Solution;
    switch(ShapeFunction)
    {
    case 0:
        Solution = iteration_rigid_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);	
        break;
    case 1:
         Solution = iteration_affine_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
		break;
    case 2:
         Solution = iteration_irregular_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    case 3:
         Solution = iteration_quadratic_LM_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    default:
        std::cout << "Shape Function Incorrect" << std::endl;
        break;
    }
		Solution.push_back(k);
		SolutionNeighbours.push_back(Solution);
}
/*--------------------------------------------------------------------------*/
void iteration_Thread(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, std::vector<Points_With_Value> &Locations_Best_Correlation)
{
	std::vector<double> Solution;
    switch(ShapeFunction)
    {
    case 0:
        Solution = iteration_rigid_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);	
        break;
    case 1:
         Solution = iteration_affine_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
		break;
    case 2:
         Solution = iteration_irregular_LM(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    case 3:
         Solution = iteration_quadratic_LM_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold);
        break;
    default:
        std::cout << "Shape Function Incorrect" << std::endl;
        break;
    }
	cv::Point NeighbourLocation(i, j);
	DispX.at<double>(NeighbourLocation) = Solution[0];
	DispY.at<double>(NeighbourLocation) = Solution[1];
	Ux.at<double>(NeighbourLocation) = Solution[2];
	Vx.at<double>(NeighbourLocation) = Solution[3];
	Uy.at<double>(NeighbourLocation) = Solution[4];
	Vy.at<double>(NeighbourLocation) = Solution[5];
	Uxy.at<double>(NeighbourLocation) = Solution[6];
	Vxy.at<double>(NeighbourLocation) = Solution[7];
	Uxx.at<double>(NeighbourLocation) = Solution[8];
	Vxx.at<double>(NeighbourLocation) = Solution[9];
	Uyy.at<double>(NeighbourLocation) = Solution[10];
	Vyy.at<double>(NeighbourLocation) = Solution[11];
	CorrelationCoefficient.at<double>(NeighbourLocation) = Solution.back();
	Computed_Points.at<uchar>(NeighbourLocation) = 1;
	Locations_Best_Correlation.push_back(Points_With_Value(Solution.back(), NeighbourLocation));
}
/*--------------------------------------------------------------------------*/
double get_fm(const cv::Mat &img,const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength)
{
	double fm = 0;
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			fm += (double)img.at<float>(jj,ii);
		}
	}
	fm /= (2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1);	
	return fm;
}
/*--------------------------------------------------------------------------*/
double get_sum_f_minus_fm_squared(const cv::Mat &img, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength, const double &fm)
{
	double sum_f_minus_fm_squared = 0;
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			sum_f_minus_fm_squared += ((double)img.at<float>(jj,ii)-fm)*((double)img.at<float>(jj,ii)-fm);
		}
	}	
	return sum_f_minus_fm_squared;
}
/*--------------------------------------------------------------------------*/
std::pair<double, double> get_gm_and_sum_g_minus_gm_squared(float *fptr_img1, const std::vector<double> &P, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength, const unsigned int &SplineDegree, const unsigned int &img_cols, const unsigned int &img_rows)
{
	double gm = 0;
	double sum_g_minus_gm_squared = 0;
    double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }	
	std::vector<double> g_values((2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1));
	unsigned int k = 0;
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		double deltax = ((double)ii-(double)Indexi);
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			double deltay = ((double)jj-(double)Indexj);
			g_values[k] = InterpolatedValue(fptr_img1, img_cols, img_rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree);
			k++;
		}
	}	
	gm = std::accumulate(g_values.begin(), g_values.end(), 0.0)/g_values.size();
	for(auto it = g_values.begin(); it != g_values.end(); ++it)
	{
		sum_g_minus_gm_squared += ((*it)-gm)*((*it)-gm);
	}	
	return std::make_pair(gm, sum_g_minus_gm_squared);   
}
/*--------------------------------------------------------------------------*/
double get_gm(float *fptr_img1, const std::vector<double> &P, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength, const unsigned int &SplineDegree, const unsigned int &img_cols, const unsigned int &img_rows)
{
	double gm = 0;
    double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }	
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		double deltax = ((double)ii-(double)Indexi);
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			double deltay = ((double)jj-(double)Indexj);
			gm += InterpolatedValue(fptr_img1, img_cols, img_rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree);					
		}
	}
	gm /= (2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1);	
	return gm;
}
/*--------------------------------------------------------------------------*/
std::pair<double, double> get_gm_and_sum_g_minus_gm_squared_short(const std::vector<double> g_values)
{
	double gm = std::accumulate(g_values.begin(), g_values.end(), 0.0)/g_values.size();
	double sum_g_minus_gm_squared = 0.0;
	for(auto it = g_values.begin(); it != g_values.end(); ++it)
	{
		sum_g_minus_gm_squared += ((*it)-gm)*((*it)-gm);
	}	
	return std::make_pair(gm, sum_g_minus_gm_squared);   
}
/*--------------------------------------------------------------------------*/
std::vector<double> get_f_values(const cv::Mat &img,const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength)
{
	std::vector<double> f_values((2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1));
	unsigned int k = 0;
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			f_values[k] = (double)img.at<float>(jj,ii);
			k++;
		}
	}
	return f_values;
}
/*--------------------------------------------------------------------------*/
std::vector<double> get_g_values(float *fptr_img1, const std::vector<double> &P, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength, const unsigned int &SplineDegree, const unsigned int &img_cols, const unsigned int &img_rows)
{
	std::vector<double> g_values((2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1));
    double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }	
	unsigned int k = 0;
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		double deltax = ((double)ii-(double)Indexi);
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			double deltay = ((double)jj-(double)Indexj);
			g_values[k] = InterpolatedValue(fptr_img1, img_cols, img_rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree);
			k++;
		}
	}
	return g_values;
}
/*--------------------------------------------------------------------------*/
double get_sum_g_minus_gm_squared(float *fptr_img1, const std::vector<double> &P, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength, const unsigned int &SplineDegree, const double &gm, const unsigned int &img_cols, const unsigned int &img_rows)
{
	double sum_g_minus_gm_squared = 0;
    double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }	
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		double deltax = ((double)ii-(double)Indexi);
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			double deltay = ((double)jj-(double)Indexj);
			double g = InterpolatedValue(fptr_img1, img_cols, img_rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree);
			sum_g_minus_gm_squared += (g-gm)*(g-gm);			
		}
	}
	return sum_g_minus_gm_squared;
}
/*--------------------------------------------------------------------------*/
bool is_data_uniform(const double &sum_g_minus_gm_squared)	
{	
	if (sum_g_minus_gm_squared < 1e-8)
	{
		return 1;
	}	
	else 
	{
		return 0;
	}
}
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_quadratic_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold)
{
		// To Be Removed
        double abs_tolerance = 1;
        double rel_tolerance = 1;
        double max_val_nu = 1e7;
        unsigned int max_iterations = 1e4;
        unsigned int iterations = 0;
        double U0 = P0[0];
        double V0 = P0[1];
        double Ux = P0[2];
        double Vx = P0[3];
        double Uy = P0[4];
        double Vy = P0[5];
        double Uxy = P0[6];
        double Vxy = P0[7];
        double Uxx = P0[8];
        double Vxx = P0[9];
        double Uyy = P0[10];
        double Vyy = P0[11];

        double lambda = 1e-10;
        double nu = 2.0;
        double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength);
        Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;

        while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold))
        {
            iterations++;
            unsigned int Indexi = SubsetLength/2 + i * GridLength;
            unsigned int Indexj = SubsetLength/2 + j * GridLength;
            double sum_g_minus_gm_squared = 0;
            double gm = 0;
            double fm = 0;
            double sum_f_minus_fm_squared = 0;
            for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
            {

                for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
                {
                    double deltax = ((double)ii-(double)Indexi);
                    double deltay = ((double)jj-(double)Indexj);
                    fm += (double)img.at<float>(jj,ii);
                    gm += InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree);
                }
            }
            gm /= (2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1);
            fm /= (2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1);
            for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
            {
                for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
                {
                    double deltax = ((double)ii-(double)Indexi);
                    double deltay = ((double)jj-(double)Indexj);
                    sum_g_minus_gm_squared += (InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm)
                                             *(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm);
                    sum_f_minus_fm_squared += ((double)img.at<float>(jj,ii)-fm)*((double)img.at<float>(jj,ii)-fm);
                }
            }
            cv::Mat Hessian(12, 12,CV_64F, cv::Scalar(0));
            cv::Mat Jacobian(12,1, CV_64F, cv::Scalar(0));
            for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
            {
                for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
                {
                    double deltax = ((double)ii-(double)Indexi);
                    double deltay = ((double)jj-(double)Indexj);
                    double xder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree, 0);
                    double yder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree, 1);

                    /// Hessian Matrix  =
                    /// [dg/dU0  dg/dU0, dg/dU0  dg/dV0, dg/dU0  dg/dUx, dg/dU0  dg/dVx, dg/dU0  dg/Uy, dg/dU0  dg/Vy, dg/dU0  dg/dUxy, dg/dU0  dg/dVxy, dg/dU0  dg/dUxx, dg/dU0  dg/dVxx, dg/dU0  dg/dUyy, dg/dU0  dg/dVyy;
                    ///  dg/dV0  dg/dU0, dg/dV0  dg/dV0, dg/dV0  dg/dUx, dg/dV0  dg/dVx, dg/dV0  dg/Uy, dg/dV0  dg/Vy, dg/dV0  dg/dUxy, dg/dV0  dg/dVxy, dg/dV0  dg/dUxx, dg/dV0  dg/dVxx, dg/dV0  dg/dUyy, dg/dV0  dg/dVyy;
                    ///  dg/dUx  dg/dU0, dg/dUx  dg/dV0, dg/dUx  dg/dUx, dg/dUx  dg/dVx, dg/dUx  dg/Uy, dg/dUx  dg/Vy, dg/dUx  dg/dUxy, dg/dUx  dg/dVxy, dg/dUx  dg/dUxx, dg/dUx  dg/dVxx, dg/dUx  dg/dUyy, dg/dUx  dg/dVyy;
                    ///  dg/dVx  dg/dU0, dg/dVx  dg/dV0, dg/dVx  dg/dUx, dg/dVx  dg/dVx, dg/dVx  dg/Uy, dg/dVx  dg/Vy, dg/dVx  dg/dUxy, dg/dVx  dg/dVxy, dg/dVx  dg/dUxx, dg/dVx  dg/dVxx, dg/dVx  dg/dUyy, dg/dVx  dg/dVyy;
                    ///  dg/dUy  dg/dU0, dg/dUy  dg/dV0, dg/dUy  dg/dUx, dg/dUy  dg/dVx, dg/dUy  dg/Uy, dg/dUy  dg/Vy, dg/dUy  dg/dUxy, dg/dUy  dg/dVxy, dg/dUy  dg/dUxx, dg/dUy  dg/dVxx, dg/dUy  dg/dUyy, dg/dUy  dg/dVyy;
                    ///  dg/dVy  dg/dU0, dg/dVy  dg/dV0, dg/dVy  dg/dUx, dg/dVy  dg/dVx, dg/dVy  dg/Uy, dg/dVy  dg/Vy, dg/dVy  dg/dUxy, dg/dVy  dg/dVxy, dg/dVy  dg/dUxx, dg/dVy  dg/dVxx, dg/dVy  dg/dUyy, dg/dVy  dg/dVyy;
                    ///  dg/dUxy dg/dU0, dg/dUxy dg/dV0, dg/dUxy dg/dUx, dg/dUxy dg/dVx, dg/dUxy dg/Uy, dg/dUxy dg/Vy, dg/dUxy dg/dUxy, dg/dUxy dg/dVxy, dg/dUxy dg/dUxx, dg/dUxy dg/dVxx, dg/dUxy dg/dUyy, dg/dUxy dg/dVyy;
                    ///  dg/dVxy dg/dU0, dg/dVxy dg/dV0, dg/dVxy dg/dUx, dg/dVxy dg/dVx, dg/dVxy dg/Uy, dg/dVxy dg/Vy, dg/dVxy dg/dUxy, dg/dVxy dg/dVxy, dg/dVxy dg/dUxx, dg/dVxy dg/dVxx, dg/dVxy dg/dUyy, dg/dvxy dg/dVyy;
                    ///  dg/dUxx dg/dU0, dg/dUxx dg/dV0, dg/dUxx dg/dUx, dg/dUxx dg/dVx, dg/dUxx dg/Uy, dg/dUxx dg/Vy, dg/dUxx dg/dUxy, dg/dUxx dg/dVxy, dg/dUxx dg/dUxx, dg/dUxx dg/dVxx, dg/dUxx dg/dUyy, dg/dUxx dg/dVyy;
                    ///  dg/dVxx dg/dU0, dg/dVxx dg/dV0, dg/dVxx dg/dUx, dg/dVxx dg/dVx, dg/dVxx dg/Uy, dg/dVxx dg/Vy, dg/dVxx dg/dUxy, dg/dVxx dg/dVxy, dg/dVxx dg/dUxx, dg/dVxx dg/dVxx, dg/dVxx dg/dUyy, dg/dVxx dg/dVyy;
                    ///  dg/dUyy dg/dU0, dg/dUyy dg/dV0, dg/dUyy dg/dUx, dg/dUyy dg/dVx, dg/dUyy dg/Uy, dg/dUyy dg/Vy, dg/dUyy dg/dUxy, dg/dUyy dg/dVxy, dg/dUyy dg/dUxx, dg/dUyy dg/dVxx, dg/dUyy dg/dUyy, dg/dUyy dg/dVyy;
                    ///  dg/dVyy dg/dU0, dg/dVyy dg/dV0, dg/dVyy dg/dUx, dg/dVyy dg/dVx, dg/dVyy dg/Uy, dg/dVyy dg/Vy, dg/dVyy dg/dUxy, dg/dVyy dg/dVxy, dg/dVyy dg/dUxx, dg/dVyy dg/dVxx, dg/dVyy dg/dUyy, dg/dVyy dg/dVyy;

                    Hessian.at<double>(0,0) += xder*xder*(1.0+lambda);
                    Hessian.at<double>(0,1) += xder*yder;
                    Hessian.at<double>(0,2) += xder*deltax*xder;
                    Hessian.at<double>(0,3) += xder*deltax*yder;
                    Hessian.at<double>(0,4) += xder*deltay*xder;
                    Hessian.at<double>(0,5) += xder*deltay*yder;
                    Hessian.at<double>(0,6) += xder*deltax*deltay*xder;
                    Hessian.at<double>(0,7) += xder*deltax*deltay*yder;
                    Hessian.at<double>(0,8) += xder*deltax*deltax*xder;
                    Hessian.at<double>(0,9) += xder*deltax*deltax*yder;
                    Hessian.at<double>(0,10) += xder*deltay*deltay*xder;
                    Hessian.at<double>(0,11) += xder*deltay*deltay*yder;

                    Hessian.at<double>(1,0) += yder*xder;
                    Hessian.at<double>(1,1) += yder*yder*(1.0+lambda);
                    Hessian.at<double>(1,2) += yder*deltax*xder;
                    Hessian.at<double>(1,3) += yder*deltax*yder;
                    Hessian.at<double>(1,4) += yder*deltay*xder;
                    Hessian.at<double>(1,5) += yder*deltay*yder;
                    Hessian.at<double>(1,6) += yder*deltax*deltay*xder;
                    Hessian.at<double>(1,7) += yder*deltax*deltay*yder;
                    Hessian.at<double>(1,8) += yder*deltax*deltax*xder;
                    Hessian.at<double>(1,9) += yder*deltax*deltax*yder;
                    Hessian.at<double>(1,10) += yder*deltay*deltay*xder;
                    Hessian.at<double>(1,11) += yder*deltay*deltay*yder;

                    Hessian.at<double>(2,0) += deltax*xder*xder;
                    Hessian.at<double>(2,1) += deltax*xder*yder;
                    Hessian.at<double>(2,2) += deltax*xder*deltax*xder*(1.0+lambda);
                    Hessian.at<double>(2,3) += deltax*xder*deltax*yder;
                    Hessian.at<double>(2,4) += deltax*xder*deltay*xder;
                    Hessian.at<double>(2,5) += deltax*xder*deltay*yder;
                    Hessian.at<double>(2,6) += deltax*xder*deltax*deltay*xder;
                    Hessian.at<double>(2,7) += deltax*xder*deltax*deltay*yder;
                    Hessian.at<double>(2,8) += deltax*xder*deltax*deltax*xder;
                    Hessian.at<double>(2,9) += deltax*xder*deltax*deltax*yder;
                    Hessian.at<double>(2,10) += deltax*xder*deltay*deltay*xder;
                    Hessian.at<double>(2,11) += deltax*xder*deltay*deltay*yder;

                    Hessian.at<double>(3,0) += deltax*yder*xder;
                    Hessian.at<double>(3,1) += deltax*yder*yder;
                    Hessian.at<double>(3,2) += deltax*yder*deltax*xder;
                    Hessian.at<double>(3,3) += deltax*yder*deltax*yder*(1.0+lambda);
                    Hessian.at<double>(3,4) += deltax*yder*deltay*xder;
                    Hessian.at<double>(3,5) += deltax*yder*deltay*yder;
                    Hessian.at<double>(3,6) += deltax*yder*deltax*deltay*xder;
                    Hessian.at<double>(3,7) += deltax*yder*deltax*deltay*yder;
                    Hessian.at<double>(3,8) += deltax*yder*deltax*deltax*xder;
                    Hessian.at<double>(3,9) += deltax*yder*deltax*deltax*yder;
                    Hessian.at<double>(3,10) += deltax*yder*deltay*deltay*xder;
                    Hessian.at<double>(3,11) += deltax*yder*deltay*deltay*yder;

                    Hessian.at<double>(4,0) += deltay*xder*xder;
                    Hessian.at<double>(4,1) += deltay*xder*yder;
                    Hessian.at<double>(4,2) += deltay*xder*deltax*xder;
                    Hessian.at<double>(4,3) += deltay*xder*deltax*yder;
                    Hessian.at<double>(4,4) += deltay*xder*deltay*xder*(1.0+lambda);
                    Hessian.at<double>(4,5) += deltay*xder*deltay*yder;
                    Hessian.at<double>(4,6) += deltay*xder*deltax*deltay*xder;
                    Hessian.at<double>(4,7) += deltay*xder*deltax*deltay*yder;
                    Hessian.at<double>(4,8) += deltay*xder*deltax*deltax*xder;
                    Hessian.at<double>(4,9) += deltay*xder*deltax*deltax*yder;
                    Hessian.at<double>(4,10) += deltay*xder*deltay*deltay*xder;
                    Hessian.at<double>(4,11) += deltay*xder*deltay*deltay*yder;

                    Hessian.at<double>(5,0) += deltay*yder*xder;
                    Hessian.at<double>(5,1) += deltay*yder*yder;
                    Hessian.at<double>(5,2) += deltay*yder*deltax*xder;
                    Hessian.at<double>(5,3) += deltay*yder*deltax*yder;
                    Hessian.at<double>(5,4) += deltay*yder*deltay*xder;
                    Hessian.at<double>(5,5) += deltay*yder*deltay*yder*(1.0+lambda);
                    Hessian.at<double>(5,6) += deltay*yder*deltax*deltay*xder;
                    Hessian.at<double>(5,7) += deltay*yder*deltax*deltay*yder;
                    Hessian.at<double>(5,8) += deltay*yder*deltax*deltax*xder;
                    Hessian.at<double>(5,9) += deltay*yder*deltax*deltax*yder;
                    Hessian.at<double>(5,10) += deltay*yder*deltay*deltay*xder;
                    Hessian.at<double>(5,11) += deltay*yder*deltay*deltay*yder;

                    Hessian.at<double>(6,0) += deltax*deltay*xder*xder;
                    Hessian.at<double>(6,1) += deltax*deltay*xder*yder;
                    Hessian.at<double>(6,2) += deltax*deltay*xder*deltax*xder;
                    Hessian.at<double>(6,3) += deltax*deltay*xder*deltax*yder;
                    Hessian.at<double>(6,4) += deltax*deltay*xder*deltay*xder;
                    Hessian.at<double>(6,5) += deltax*deltay*xder*deltay*yder;
                    Hessian.at<double>(6,6) += deltax*deltay*xder*deltax*deltay*xder*(1.0+lambda);
                    Hessian.at<double>(6,7) += deltax*deltay*xder*deltax*deltay*yder;
                    Hessian.at<double>(6,8) += deltax*deltay*xder*deltax*deltax*xder;
                    Hessian.at<double>(6,9) += deltax*deltay*xder*deltax*deltax*yder;
                    Hessian.at<double>(6,10) += deltax*deltay*xder*deltay*deltay*xder;
                    Hessian.at<double>(6,11) += deltax*deltay*xder*deltay*deltay*yder;

                    Hessian.at<double>(7,0) += deltax*deltay*yder*xder;
                    Hessian.at<double>(7,1) += deltax*deltay*yder*yder;
                    Hessian.at<double>(7,2) += deltax*deltay*yder*deltax*xder;
                    Hessian.at<double>(7,3) += deltax*deltay*yder*deltax*yder;
                    Hessian.at<double>(7,4) += deltax*deltay*yder*deltay*xder;
                    Hessian.at<double>(7,5) += deltax*deltay*yder*deltay*yder;
                    Hessian.at<double>(7,6) += deltax*deltay*yder*deltax*deltay*xder;
                    Hessian.at<double>(7,7) += deltax*deltay*yder*deltax*deltay*yder*(1.0+lambda);
                    Hessian.at<double>(7,8) += deltax*deltay*yder*deltax*deltax*xder;
                    Hessian.at<double>(7,9) += deltax*deltay*yder*deltax*deltax*yder;
                    Hessian.at<double>(7,10) += deltax*deltay*yder*deltay*deltay*xder;
                    Hessian.at<double>(7,11) += deltax*deltay*yder*deltay*deltay*yder;

                    Hessian.at<double>(8,0) += deltax*deltax*xder*xder;
                    Hessian.at<double>(8,1) += deltax*deltax*xder*yder;
                    Hessian.at<double>(8,2) += deltax*deltax*xder*deltax*xder;
                    Hessian.at<double>(8,3) += deltax*deltax*xder*deltax*yder;
                    Hessian.at<double>(8,4) += deltax*deltax*xder*deltay*xder;
                    Hessian.at<double>(8,5) += deltax*deltax*xder*deltay*yder;
                    Hessian.at<double>(8,6) += deltax*deltax*xder*deltax*deltay*xder;
                    Hessian.at<double>(8,7) += deltax*deltax*xder*deltax*deltay*yder;
                    Hessian.at<double>(8,8) += deltax*deltax*xder*deltax*deltax*xder*(1.0+lambda);
                    Hessian.at<double>(8,9) += deltax*deltax*xder*deltax*deltax*yder;
                    Hessian.at<double>(8,10) += deltax*deltax*xder*deltay*deltay*xder;
                    Hessian.at<double>(8,11) += deltax*deltax*xder*deltay*deltay*yder;

                    Hessian.at<double>(9,0) += deltax*deltax*yder*xder;
                    Hessian.at<double>(9,1) += deltax*deltax*yder*yder;
                    Hessian.at<double>(9,2) += deltax*deltax*yder*deltax*xder;
                    Hessian.at<double>(9,3) += deltax*deltax*yder*deltax*yder;
                    Hessian.at<double>(9,4) += deltax*deltax*yder*deltay*xder;
                    Hessian.at<double>(9,5) += deltax*deltax*yder*deltay*yder;
                    Hessian.at<double>(9,6) += deltax*deltax*yder*deltax*deltay*xder;
                    Hessian.at<double>(9,7) += deltax*deltax*yder*deltax*deltay*yder;
                    Hessian.at<double>(9,8) += deltax*deltax*yder*deltax*deltax*xder;
                    Hessian.at<double>(9,9) += deltax*deltax*yder*deltax*deltax*yder*(1.0+lambda);
                    Hessian.at<double>(9,10) += deltax*deltax*yder*deltay*deltay*xder;
                    Hessian.at<double>(9,11) += deltax*deltax*yder*deltay*deltay*yder;

                    Hessian.at<double>(10,0) += deltay*deltay*xder*xder;
                    Hessian.at<double>(10,1) += deltay*deltay*xder*yder;
                    Hessian.at<double>(10,2) += deltay*deltay*xder*deltax*xder;
                    Hessian.at<double>(10,3) += deltay*deltay*xder*deltax*yder;
                    Hessian.at<double>(10,4) += deltay*deltay*xder*deltay*xder;
                    Hessian.at<double>(10,5) += deltay*deltay*xder*deltay*yder;
                    Hessian.at<double>(10,6) += deltay*deltay*xder*deltax*deltay*xder;
                    Hessian.at<double>(10,7) += deltay*deltay*xder*deltax*deltay*yder;
                    Hessian.at<double>(10,8) += deltay*deltay*xder*deltax*deltax*xder;
                    Hessian.at<double>(10,9) += deltay*deltay*xder*deltax*deltax*yder;
                    Hessian.at<double>(10,10) += deltay*deltay*xder*deltay*deltay*xder*(1.0+lambda);
                    Hessian.at<double>(10,11) += deltay*deltay*xder*deltay*deltay*yder;

                    Hessian.at<double>(11,0) += deltay*deltay*yder*xder;
                    Hessian.at<double>(11,1) += deltay*deltay*yder*yder;
                    Hessian.at<double>(11,2) += deltay*deltay*yder*deltax*xder;
                    Hessian.at<double>(11,3) += deltay*deltay*yder*deltax*yder;
                    Hessian.at<double>(11,4) += deltay*deltay*yder*deltay*xder;
                    Hessian.at<double>(11,5) += deltay*deltay*yder*deltay*yder;
                    Hessian.at<double>(11,6) += deltay*deltay*yder*deltax*deltay*xder;
                    Hessian.at<double>(11,7) += deltay*deltay*yder*deltax*deltay*yder;
                    Hessian.at<double>(11,8) += deltay*deltay*yder*deltax*deltax*xder;
                    Hessian.at<double>(11,9) += deltay*deltay*yder*deltax*deltax*yder;
                    Hessian.at<double>(11,10) += deltay*deltay*yder*deltay*deltay*xder;
                    Hessian.at<double>(11,11) += deltay*deltay*yder*deltay*deltay*yder*(1.0+lambda);

                    /// Jacobian Matrix [dg/dU0; dg/dV0; dg/dUx; dg/dVx; dg/dUy; dg/dVy; dg/dUxy; dg/dVxy; dg/dUxx; dg/dVxx; dg/dUyy; dg/dVyy]
                    double C_components = (img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared);
                    Jacobian.at<double>(0,0) += xder*C_components;
                    Jacobian.at<double>(1,0) += yder*C_components;
                    Jacobian.at<double>(2,0) += deltax*xder*C_components;
                    Jacobian.at<double>(3,0) += deltax*yder*C_components;
                    Jacobian.at<double>(4,0) += deltay*xder*C_components;
                    Jacobian.at<double>(5,0) += deltay*yder*C_components;
                    Jacobian.at<double>(6,0) += deltax*deltay*xder*C_components;
                    Jacobian.at<double>(7,0) += deltax*deltay*yder*C_components;
                    Jacobian.at<double>(8,0) += deltax*deltax*xder*C_components;
                    Jacobian.at<double>(9,0) += deltax*deltax*yder*C_components;
                    Jacobian.at<double>(10,0) += deltay*deltay*xder*C_components;
                    Jacobian.at<double>(11,0) += deltay*deltay*yder*C_components;
                }
            }
            /// Suggestion: Multiply both side with sqrt(sum_g_minus_gm_squared)/2
            Jacobian *= 2/sqrt(sum_g_minus_gm_squared);
            Hessian *= 2/sum_g_minus_gm_squared;

            cv::Mat X(12,1,CV_64F);
            cv::solve(Hessian, Jacobian, X, cv::DECOMP_CHOLESKY);

            std::vector<double> Suggested_Solution = {U0+X.at<double>(0), V0+X.at<double>(1), Ux+X.at<double>(2), Vx+X.at<double>(3), Uy+X.at<double>(4), Vy+X.at<double>(5),  Uxy+X.at<double>(6), Vxy+X.at<double>(7), Uxx+X.at<double>(8), Vxx+X.at<double>(9), Uyy+X.at<double>(10), Vyy+X.at<double>(11)};
            double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength);
            Suggested_Correlation_Coefficient = 1.0-0.5*Suggested_Correlation_Coefficient;

            if (Suggested_Correlation_Coefficient > Correlation_Coefficient_old)
            {
                // Improvement
                Correlation_Coefficient_old = Suggested_Correlation_Coefficient;
                abs_tolerance = cv::sum(cv::abs(X)).val[0];
                rel_tolerance = cv::sum(cv::abs(X)).val[0]/(abs(P0[0])+abs(P0[1])+abs(P0[2])+abs(P0[3])+abs(P0[4])+abs(P0[5])+abs(P0[6])+abs(P0[7])+abs(P0[8])+abs(P0[9])+abs(P0[10])+abs(P0[11]));

                U0 += X.at<double>(0);
                V0 += X.at<double>(1);
                Ux += X.at<double>(2);
                Vx += X.at<double>(3);
                Uy += X.at<double>(4);
                Vy += X.at<double>(5);
                Uxy += X.at<double>(6);
                Vxy += X.at<double>(7);
                Uxx += X.at<double>(8);
                Vxx += X.at<double>(9);
                Uyy += X.at<double>(10);
                Vyy += X.at<double>(11);
                // Decrease lambda => More Gauss-Newton, less Gradient Search
                nu = 2;
                lambda /= 3.0;
            }
            else
            {
                // No Improvement
                // Increase lambda => Less Gauss-Newton, more Gradient Search
                lambda *= nu;
                nu *= 2.0;
            }
        }
        //if (nu>max_val_nu)
        //    std::cout << "Stuck" << std::endl;
        if (iterations>max_iterations)
            std::cout << "Max Iterations Reached" << std::endl;

        std::vector<double> returnvector;
        returnvector.push_back(U0);
        returnvector.push_back(V0);
        returnvector.push_back(Ux);
        returnvector.push_back(Vx);
        returnvector.push_back(Uy);
        returnvector.push_back(Vy);
        returnvector.push_back(Uxy);
        returnvector.push_back(Vxy);
        returnvector.push_back(Uxx);
        returnvector.push_back(Vxx);
        returnvector.push_back(Uyy);
        returnvector.push_back(Vyy);
        double Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, returnvector, SplineDegree, SubsetLength, GridLength);
        Correlation_Coefficient = 1.0-0.5*Correlation_Coefficient;
        returnvector.push_back(Correlation_Coefficient);

        return returnvector;
}
/*--------------------------------------------------------------------------*/
void calculate_Hessian_Jacobian(cv::Mat &Hessian, cv::Mat &Jacobian, const cv::Mat &img, float *fptr_img1, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &ShapeFunction, const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &sum_g_minus_gm_squared, const double &lambda)
{   
	double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		double deltax = ((double)ii-(double)Indexi);
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			double deltay = ((double)jj-(double)Indexj);
			double xder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree, 0);
			double yder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree, 1);

			/// Hessian Matrix  =
			/// [dg/dU0  dg/dU0, dg/dU0  dg/dV0, dg/dU0  dg/dUx, dg/dU0  dg/dVx, dg/dU0  dg/Uy, dg/dU0  dg/Vy, dg/dU0  dg/dUxy, dg/dU0  dg/dVxy, dg/dU0  dg/dUxx, dg/dU0  dg/dVxx, dg/dU0  dg/dUyy, dg/dU0  dg/dVyy;
			///  dg/dV0  dg/dU0, dg/dV0  dg/dV0, dg/dV0  dg/dUx, dg/dV0  dg/dVx, dg/dV0  dg/Uy, dg/dV0  dg/Vy, dg/dV0  dg/dUxy, dg/dV0  dg/dVxy, dg/dV0  dg/dUxx, dg/dV0  dg/dVxx, dg/dV0  dg/dUyy, dg/dV0  dg/dVyy;
			///  dg/dUx  dg/dU0, dg/dUx  dg/dV0, dg/dUx  dg/dUx, dg/dUx  dg/dVx, dg/dUx  dg/Uy, dg/dUx  dg/Vy, dg/dUx  dg/dUxy, dg/dUx  dg/dVxy, dg/dUx  dg/dUxx, dg/dUx  dg/dVxx, dg/dUx  dg/dUyy, dg/dUx  dg/dVyy;
			///  dg/dVx  dg/dU0, dg/dVx  dg/dV0, dg/dVx  dg/dUx, dg/dVx  dg/dVx, dg/dVx  dg/Uy, dg/dVx  dg/Vy, dg/dVx  dg/dUxy, dg/dVx  dg/dVxy, dg/dVx  dg/dUxx, dg/dVx  dg/dVxx, dg/dVx  dg/dUyy, dg/dVx  dg/dVyy;
			///  dg/dUy  dg/dU0, dg/dUy  dg/dV0, dg/dUy  dg/dUx, dg/dUy  dg/dVx, dg/dUy  dg/Uy, dg/dUy  dg/Vy, dg/dUy  dg/dUxy, dg/dUy  dg/dVxy, dg/dUy  dg/dUxx, dg/dUy  dg/dVxx, dg/dUy  dg/dUyy, dg/dUy  dg/dVyy;
			///  dg/dVy  dg/dU0, dg/dVy  dg/dV0, dg/dVy  dg/dUx, dg/dVy  dg/dVx, dg/dVy  dg/Uy, dg/dVy  dg/Vy, dg/dVy  dg/dUxy, dg/dVy  dg/dVxy, dg/dVy  dg/dUxx, dg/dVy  dg/dVxx, dg/dVy  dg/dUyy, dg/dVy  dg/dVyy;
			///  dg/dUxy dg/dU0, dg/dUxy dg/dV0, dg/dUxy dg/dUx, dg/dUxy dg/dVx, dg/dUxy dg/Uy, dg/dUxy dg/Vy, dg/dUxy dg/dUxy, dg/dUxy dg/dVxy, dg/dUxy dg/dUxx, dg/dUxy dg/dVxx, dg/dUxy dg/dUyy, dg/dUxy dg/dVyy;
			///  dg/dVxy dg/dU0, dg/dVxy dg/dV0, dg/dVxy dg/dUx, dg/dVxy dg/dVx, dg/dVxy dg/Uy, dg/dVxy dg/Vy, dg/dVxy dg/dUxy, dg/dVxy dg/dVxy, dg/dVxy dg/dUxx, dg/dVxy dg/dVxx, dg/dVxy dg/dUyy, dg/dvxy dg/dVyy;
			///  dg/dUxx dg/dU0, dg/dUxx dg/dV0, dg/dUxx dg/dUx, dg/dUxx dg/dVx, dg/dUxx dg/Uy, dg/dUxx dg/Vy, dg/dUxx dg/dUxy, dg/dUxx dg/dVxy, dg/dUxx dg/dUxx, dg/dUxx dg/dVxx, dg/dUxx dg/dUyy, dg/dUxx dg/dVyy;
			///  dg/dVxx dg/dU0, dg/dVxx dg/dV0, dg/dVxx dg/dUx, dg/dVxx dg/dVx, dg/dVxx dg/Uy, dg/dVxx dg/Vy, dg/dVxx dg/dUxy, dg/dVxx dg/dVxy, dg/dVxx dg/dUxx, dg/dVxx dg/dVxx, dg/dVxx dg/dUyy, dg/dVxx dg/dVyy;
			///  dg/dUyy dg/dU0, dg/dUyy dg/dV0, dg/dUyy dg/dUx, dg/dUyy dg/dVx, dg/dUyy dg/Uy, dg/dUyy dg/Vy, dg/dUyy dg/dUxy, dg/dUyy dg/dVxy, dg/dUyy dg/dUxx, dg/dUyy dg/dVxx, dg/dUyy dg/dUyy, dg/dUyy dg/dVyy;
			///  dg/dVyy dg/dU0, dg/dVyy dg/dV0, dg/dVyy dg/dUx, dg/dVyy dg/dVx, dg/dVyy dg/Uy, dg/dVyy dg/Vy, dg/dVyy dg/dUxy, dg/dVyy dg/dVxy, dg/dVyy dg/dUxx, dg/dVyy dg/dVxx, dg/dVyy dg/dUyy, dg/dVyy dg/dVyy;

			Hessian.at<double>(0,0) += xder*xder*(1.0+lambda);
			Hessian.at<double>(0,1) += xder*yder;
			Hessian.at<double>(0,2) += xder*deltax*xder;
			Hessian.at<double>(0,3) += xder*deltax*yder;
			Hessian.at<double>(0,4) += xder*deltay*xder;
			Hessian.at<double>(0,5) += xder*deltay*yder;
			Hessian.at<double>(0,6) += xder*deltax*deltay*xder;
			Hessian.at<double>(0,7) += xder*deltax*deltay*yder;
			Hessian.at<double>(0,8) += xder*deltax*deltax*xder;
			Hessian.at<double>(0,9) += xder*deltax*deltax*yder;
			Hessian.at<double>(0,10) += xder*deltay*deltay*xder;
			Hessian.at<double>(0,11) += xder*deltay*deltay*yder;

			Hessian.at<double>(1,0) += yder*xder;
			Hessian.at<double>(1,1) += yder*yder*(1.0+lambda);
			Hessian.at<double>(1,2) += yder*deltax*xder;
			Hessian.at<double>(1,3) += yder*deltax*yder;
			Hessian.at<double>(1,4) += yder*deltay*xder;
			Hessian.at<double>(1,5) += yder*deltay*yder;
			Hessian.at<double>(1,6) += yder*deltax*deltay*xder;
			Hessian.at<double>(1,7) += yder*deltax*deltay*yder;
			Hessian.at<double>(1,8) += yder*deltax*deltax*xder;
			Hessian.at<double>(1,9) += yder*deltax*deltax*yder;
			Hessian.at<double>(1,10) += yder*deltay*deltay*xder;
			Hessian.at<double>(1,11) += yder*deltay*deltay*yder;

			Hessian.at<double>(2,0) += deltax*xder*xder;
			Hessian.at<double>(2,1) += deltax*xder*yder;
			Hessian.at<double>(2,2) += deltax*xder*deltax*xder*(1.0+lambda);
			Hessian.at<double>(2,3) += deltax*xder*deltax*yder;
			Hessian.at<double>(2,4) += deltax*xder*deltay*xder;
			Hessian.at<double>(2,5) += deltax*xder*deltay*yder;
			Hessian.at<double>(2,6) += deltax*xder*deltax*deltay*xder;
			Hessian.at<double>(2,7) += deltax*xder*deltax*deltay*yder;
			Hessian.at<double>(2,8) += deltax*xder*deltax*deltax*xder;
			Hessian.at<double>(2,9) += deltax*xder*deltax*deltax*yder;
			Hessian.at<double>(2,10) += deltax*xder*deltay*deltay*xder;
			Hessian.at<double>(2,11) += deltax*xder*deltay*deltay*yder;

			Hessian.at<double>(3,0) += deltax*yder*xder;
			Hessian.at<double>(3,1) += deltax*yder*yder;
			Hessian.at<double>(3,2) += deltax*yder*deltax*xder;
			Hessian.at<double>(3,3) += deltax*yder*deltax*yder*(1.0+lambda);
			Hessian.at<double>(3,4) += deltax*yder*deltay*xder;
			Hessian.at<double>(3,5) += deltax*yder*deltay*yder;
			Hessian.at<double>(3,6) += deltax*yder*deltax*deltay*xder;
			Hessian.at<double>(3,7) += deltax*yder*deltax*deltay*yder;
			Hessian.at<double>(3,8) += deltax*yder*deltax*deltax*xder;
			Hessian.at<double>(3,9) += deltax*yder*deltax*deltax*yder;
			Hessian.at<double>(3,10) += deltax*yder*deltay*deltay*xder;
			Hessian.at<double>(3,11) += deltax*yder*deltay*deltay*yder;

			Hessian.at<double>(4,0) += deltay*xder*xder;
			Hessian.at<double>(4,1) += deltay*xder*yder;
			Hessian.at<double>(4,2) += deltay*xder*deltax*xder;
			Hessian.at<double>(4,3) += deltay*xder*deltax*yder;
			Hessian.at<double>(4,4) += deltay*xder*deltay*xder*(1.0+lambda);
			Hessian.at<double>(4,5) += deltay*xder*deltay*yder;
			Hessian.at<double>(4,6) += deltay*xder*deltax*deltay*xder;
			Hessian.at<double>(4,7) += deltay*xder*deltax*deltay*yder;
			Hessian.at<double>(4,8) += deltay*xder*deltax*deltax*xder;
			Hessian.at<double>(4,9) += deltay*xder*deltax*deltax*yder;
			Hessian.at<double>(4,10) += deltay*xder*deltay*deltay*xder;
			Hessian.at<double>(4,11) += deltay*xder*deltay*deltay*yder;

			Hessian.at<double>(5,0) += deltay*yder*xder;
			Hessian.at<double>(5,1) += deltay*yder*yder;
			Hessian.at<double>(5,2) += deltay*yder*deltax*xder;
			Hessian.at<double>(5,3) += deltay*yder*deltax*yder;
			Hessian.at<double>(5,4) += deltay*yder*deltay*xder;
			Hessian.at<double>(5,5) += deltay*yder*deltay*yder*(1.0+lambda);
			Hessian.at<double>(5,6) += deltay*yder*deltax*deltay*xder;
			Hessian.at<double>(5,7) += deltay*yder*deltax*deltay*yder;
			Hessian.at<double>(5,8) += deltay*yder*deltax*deltax*xder;
			Hessian.at<double>(5,9) += deltay*yder*deltax*deltax*yder;
			Hessian.at<double>(5,10) += deltay*yder*deltay*deltay*xder;
			Hessian.at<double>(5,11) += deltay*yder*deltay*deltay*yder;

			Hessian.at<double>(6,0) += deltax*deltay*xder*xder;
			Hessian.at<double>(6,1) += deltax*deltay*xder*yder;
			Hessian.at<double>(6,2) += deltax*deltay*xder*deltax*xder;
			Hessian.at<double>(6,3) += deltax*deltay*xder*deltax*yder;
			Hessian.at<double>(6,4) += deltax*deltay*xder*deltay*xder;
			Hessian.at<double>(6,5) += deltax*deltay*xder*deltay*yder;
			Hessian.at<double>(6,6) += deltax*deltay*xder*deltax*deltay*xder*(1.0+lambda);
			Hessian.at<double>(6,7) += deltax*deltay*xder*deltax*deltay*yder;
			Hessian.at<double>(6,8) += deltax*deltay*xder*deltax*deltax*xder;
			Hessian.at<double>(6,9) += deltax*deltay*xder*deltax*deltax*yder;
			Hessian.at<double>(6,10) += deltax*deltay*xder*deltay*deltay*xder;
			Hessian.at<double>(6,11) += deltax*deltay*xder*deltay*deltay*yder;

			Hessian.at<double>(7,0) += deltax*deltay*yder*xder;
			Hessian.at<double>(7,1) += deltax*deltay*yder*yder;
			Hessian.at<double>(7,2) += deltax*deltay*yder*deltax*xder;
			Hessian.at<double>(7,3) += deltax*deltay*yder*deltax*yder;
			Hessian.at<double>(7,4) += deltax*deltay*yder*deltay*xder;
			Hessian.at<double>(7,5) += deltax*deltay*yder*deltay*yder;
			Hessian.at<double>(7,6) += deltax*deltay*yder*deltax*deltay*xder;
			Hessian.at<double>(7,7) += deltax*deltay*yder*deltax*deltay*yder*(1.0+lambda);
			Hessian.at<double>(7,8) += deltax*deltay*yder*deltax*deltax*xder;
			Hessian.at<double>(7,9) += deltax*deltay*yder*deltax*deltax*yder;
			Hessian.at<double>(7,10) += deltax*deltay*yder*deltay*deltay*xder;
			Hessian.at<double>(7,11) += deltax*deltay*yder*deltay*deltay*yder;

			Hessian.at<double>(8,0) += deltax*deltax*xder*xder;
			Hessian.at<double>(8,1) += deltax*deltax*xder*yder;
			Hessian.at<double>(8,2) += deltax*deltax*xder*deltax*xder;
			Hessian.at<double>(8,3) += deltax*deltax*xder*deltax*yder;
			Hessian.at<double>(8,4) += deltax*deltax*xder*deltay*xder;
			Hessian.at<double>(8,5) += deltax*deltax*xder*deltay*yder;
			Hessian.at<double>(8,6) += deltax*deltax*xder*deltax*deltay*xder;
			Hessian.at<double>(8,7) += deltax*deltax*xder*deltax*deltay*yder;
			Hessian.at<double>(8,8) += deltax*deltax*xder*deltax*deltax*xder*(1.0+lambda);
			Hessian.at<double>(8,9) += deltax*deltax*xder*deltax*deltax*yder;
			Hessian.at<double>(8,10) += deltax*deltax*xder*deltay*deltay*xder;
			Hessian.at<double>(8,11) += deltax*deltax*xder*deltay*deltay*yder;

			Hessian.at<double>(9,0) += deltax*deltax*yder*xder;
			Hessian.at<double>(9,1) += deltax*deltax*yder*yder;
			Hessian.at<double>(9,2) += deltax*deltax*yder*deltax*xder;
			Hessian.at<double>(9,3) += deltax*deltax*yder*deltax*yder;
			Hessian.at<double>(9,4) += deltax*deltax*yder*deltay*xder;
			Hessian.at<double>(9,5) += deltax*deltax*yder*deltay*yder;
			Hessian.at<double>(9,6) += deltax*deltax*yder*deltax*deltay*xder;
			Hessian.at<double>(9,7) += deltax*deltax*yder*deltax*deltay*yder;
			Hessian.at<double>(9,8) += deltax*deltax*yder*deltax*deltax*xder;
			Hessian.at<double>(9,9) += deltax*deltax*yder*deltax*deltax*yder*(1.0+lambda);
			Hessian.at<double>(9,10) += deltax*deltax*yder*deltay*deltay*xder;
			Hessian.at<double>(9,11) += deltax*deltax*yder*deltay*deltay*yder;

			Hessian.at<double>(10,0) += deltay*deltay*xder*xder;
			Hessian.at<double>(10,1) += deltay*deltay*xder*yder;
			Hessian.at<double>(10,2) += deltay*deltay*xder*deltax*xder;
			Hessian.at<double>(10,3) += deltay*deltay*xder*deltax*yder;
			Hessian.at<double>(10,4) += deltay*deltay*xder*deltay*xder;
			Hessian.at<double>(10,5) += deltay*deltay*xder*deltay*yder;
			Hessian.at<double>(10,6) += deltay*deltay*xder*deltax*deltay*xder;
			Hessian.at<double>(10,7) += deltay*deltay*xder*deltax*deltay*yder;
			Hessian.at<double>(10,8) += deltay*deltay*xder*deltax*deltax*xder;
			Hessian.at<double>(10,9) += deltay*deltay*xder*deltax*deltax*yder;
			Hessian.at<double>(10,10) += deltay*deltay*xder*deltay*deltay*xder*(1.0+lambda);
			Hessian.at<double>(10,11) += deltay*deltay*xder*deltay*deltay*yder;

			Hessian.at<double>(11,0) += deltay*deltay*yder*xder;
			Hessian.at<double>(11,1) += deltay*deltay*yder*yder;
			Hessian.at<double>(11,2) += deltay*deltay*yder*deltax*xder;
			Hessian.at<double>(11,3) += deltay*deltay*yder*deltax*yder;
			Hessian.at<double>(11,4) += deltay*deltay*yder*deltay*xder;
			Hessian.at<double>(11,5) += deltay*deltay*yder*deltay*yder;
			Hessian.at<double>(11,6) += deltay*deltay*yder*deltax*deltay*xder;
			Hessian.at<double>(11,7) += deltay*deltay*yder*deltax*deltay*yder;
			Hessian.at<double>(11,8) += deltay*deltay*yder*deltax*deltax*xder;
			Hessian.at<double>(11,9) += deltay*deltay*yder*deltax*deltax*yder;
			Hessian.at<double>(11,10) += deltay*deltay*yder*deltay*deltay*xder;
			Hessian.at<double>(11,11) += deltay*deltay*yder*deltay*deltay*yder*(1.0+lambda);

			/// Jacobian Matrix [dg/dU0; dg/dV0; dg/dUx; dg/dVx; dg/dUy; dg/dVy; dg/dUxy; dg/dVxy; dg/dUxx; dg/dVxx; dg/dUyy; dg/dVyy]
			double C_components = (img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared);
			Jacobian.at<double>(0,0) += xder*C_components;
			Jacobian.at<double>(1,0) += yder*C_components;
			Jacobian.at<double>(2,0) += deltax*xder*C_components;
			Jacobian.at<double>(3,0) += deltax*yder*C_components;
			Jacobian.at<double>(4,0) += deltay*xder*C_components;
			Jacobian.at<double>(5,0) += deltay*yder*C_components;
			Jacobian.at<double>(6,0) += deltax*deltay*xder*C_components;
			Jacobian.at<double>(7,0) += deltax*deltay*yder*C_components;
			Jacobian.at<double>(8,0) += deltax*deltax*xder*C_components;
			Jacobian.at<double>(9,0) += deltax*deltax*yder*C_components;
			Jacobian.at<double>(10,0) += deltay*deltay*xder*C_components;
			Jacobian.at<double>(11,0) += deltay*deltay*yder*C_components;
		}
	}
	Hessian *= 1.0/sqrt(sum_g_minus_gm_squared);	
}
/*--------------------------------------------------------------------------*/
void calculate_Hessian_Jacobian_short(cv::Mat &Hessian, cv::Mat &Jacobian, const cv::Mat &img, float *fptr_img1, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &ShapeFunction, const std::vector<double> &f_values, const double &fm, const double &sum_f_minus_fm_squared, const std::vector<double> &g_values, const double &gm, const double &sum_g_minus_gm_squared, const double &lambda)
{   
	double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }
	unsigned int k=0;
	for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
	{
		double deltax = ((double)ii-(double)Indexi);
		for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
		{
			double deltay = ((double)jj-(double)Indexj);
			double xder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree, 0);
			double yder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree, 1);

			/// Hessian Matrix  =
			/// [dg/dU0  dg/dU0, dg/dU0  dg/dV0, dg/dU0  dg/dUx, dg/dU0  dg/dVx, dg/dU0  dg/Uy, dg/dU0  dg/Vy, dg/dU0  dg/dUxy, dg/dU0  dg/dVxy, dg/dU0  dg/dUxx, dg/dU0  dg/dVxx, dg/dU0  dg/dUyy, dg/dU0  dg/dVyy;
			///  dg/dV0  dg/dU0, dg/dV0  dg/dV0, dg/dV0  dg/dUx, dg/dV0  dg/dVx, dg/dV0  dg/Uy, dg/dV0  dg/Vy, dg/dV0  dg/dUxy, dg/dV0  dg/dVxy, dg/dV0  dg/dUxx, dg/dV0  dg/dVxx, dg/dV0  dg/dUyy, dg/dV0  dg/dVyy;
			///  dg/dUx  dg/dU0, dg/dUx  dg/dV0, dg/dUx  dg/dUx, dg/dUx  dg/dVx, dg/dUx  dg/Uy, dg/dUx  dg/Vy, dg/dUx  dg/dUxy, dg/dUx  dg/dVxy, dg/dUx  dg/dUxx, dg/dUx  dg/dVxx, dg/dUx  dg/dUyy, dg/dUx  dg/dVyy;
			///  dg/dVx  dg/dU0, dg/dVx  dg/dV0, dg/dVx  dg/dUx, dg/dVx  dg/dVx, dg/dVx  dg/Uy, dg/dVx  dg/Vy, dg/dVx  dg/dUxy, dg/dVx  dg/dVxy, dg/dVx  dg/dUxx, dg/dVx  dg/dVxx, dg/dVx  dg/dUyy, dg/dVx  dg/dVyy;
			///  dg/dUy  dg/dU0, dg/dUy  dg/dV0, dg/dUy  dg/dUx, dg/dUy  dg/dVx, dg/dUy  dg/Uy, dg/dUy  dg/Vy, dg/dUy  dg/dUxy, dg/dUy  dg/dVxy, dg/dUy  dg/dUxx, dg/dUy  dg/dVxx, dg/dUy  dg/dUyy, dg/dUy  dg/dVyy;
			///  dg/dVy  dg/dU0, dg/dVy  dg/dV0, dg/dVy  dg/dUx, dg/dVy  dg/dVx, dg/dVy  dg/Uy, dg/dVy  dg/Vy, dg/dVy  dg/dUxy, dg/dVy  dg/dVxy, dg/dVy  dg/dUxx, dg/dVy  dg/dVxx, dg/dVy  dg/dUyy, dg/dVy  dg/dVyy;
			///  dg/dUxy dg/dU0, dg/dUxy dg/dV0, dg/dUxy dg/dUx, dg/dUxy dg/dVx, dg/dUxy dg/Uy, dg/dUxy dg/Vy, dg/dUxy dg/dUxy, dg/dUxy dg/dVxy, dg/dUxy dg/dUxx, dg/dUxy dg/dVxx, dg/dUxy dg/dUyy, dg/dUxy dg/dVyy;
			///  dg/dVxy dg/dU0, dg/dVxy dg/dV0, dg/dVxy dg/dUx, dg/dVxy dg/dVx, dg/dVxy dg/Uy, dg/dVxy dg/Vy, dg/dVxy dg/dUxy, dg/dVxy dg/dVxy, dg/dVxy dg/dUxx, dg/dVxy dg/dVxx, dg/dVxy dg/dUyy, dg/dvxy dg/dVyy;
			///  dg/dUxx dg/dU0, dg/dUxx dg/dV0, dg/dUxx dg/dUx, dg/dUxx dg/dVx, dg/dUxx dg/Uy, dg/dUxx dg/Vy, dg/dUxx dg/dUxy, dg/dUxx dg/dVxy, dg/dUxx dg/dUxx, dg/dUxx dg/dVxx, dg/dUxx dg/dUyy, dg/dUxx dg/dVyy;
			///  dg/dVxx dg/dU0, dg/dVxx dg/dV0, dg/dVxx dg/dUx, dg/dVxx dg/dVx, dg/dVxx dg/Uy, dg/dVxx dg/Vy, dg/dVxx dg/dUxy, dg/dVxx dg/dVxy, dg/dVxx dg/dUxx, dg/dVxx dg/dVxx, dg/dVxx dg/dUyy, dg/dVxx dg/dVyy;
			///  dg/dUyy dg/dU0, dg/dUyy dg/dV0, dg/dUyy dg/dUx, dg/dUyy dg/dVx, dg/dUyy dg/Uy, dg/dUyy dg/Vy, dg/dUyy dg/dUxy, dg/dUyy dg/dVxy, dg/dUyy dg/dUxx, dg/dUyy dg/dVxx, dg/dUyy dg/dUyy, dg/dUyy dg/dVyy;
			///  dg/dVyy dg/dU0, dg/dVyy dg/dV0, dg/dVyy dg/dUx, dg/dVyy dg/dVx, dg/dVyy dg/Uy, dg/dVyy dg/Vy, dg/dVyy dg/dUxy, dg/dVyy dg/dVxy, dg/dVyy dg/dUxx, dg/dVyy dg/dVxx, dg/dVyy dg/dUyy, dg/dVyy dg/dVyy;

			Hessian.at<double>(0,0) += xder*xder*(1.0+lambda);
			Hessian.at<double>(0,1) += xder*yder;
			Hessian.at<double>(0,2) += xder*deltax*xder;
			Hessian.at<double>(0,3) += xder*deltax*yder;
			Hessian.at<double>(0,4) += xder*deltay*xder;
			Hessian.at<double>(0,5) += xder*deltay*yder;
			Hessian.at<double>(0,6) += xder*deltax*deltay*xder;
			Hessian.at<double>(0,7) += xder*deltax*deltay*yder;
			Hessian.at<double>(0,8) += xder*deltax*deltax*xder;
			Hessian.at<double>(0,9) += xder*deltax*deltax*yder;
			Hessian.at<double>(0,10) += xder*deltay*deltay*xder;
			Hessian.at<double>(0,11) += xder*deltay*deltay*yder;

			Hessian.at<double>(1,0) += yder*xder;
			Hessian.at<double>(1,1) += yder*yder*(1.0+lambda);
			Hessian.at<double>(1,2) += yder*deltax*xder;
			Hessian.at<double>(1,3) += yder*deltax*yder;
			Hessian.at<double>(1,4) += yder*deltay*xder;
			Hessian.at<double>(1,5) += yder*deltay*yder;
			Hessian.at<double>(1,6) += yder*deltax*deltay*xder;
			Hessian.at<double>(1,7) += yder*deltax*deltay*yder;
			Hessian.at<double>(1,8) += yder*deltax*deltax*xder;
			Hessian.at<double>(1,9) += yder*deltax*deltax*yder;
			Hessian.at<double>(1,10) += yder*deltay*deltay*xder;
			Hessian.at<double>(1,11) += yder*deltay*deltay*yder;

			Hessian.at<double>(2,0) += deltax*xder*xder;
			Hessian.at<double>(2,1) += deltax*xder*yder;
			Hessian.at<double>(2,2) += deltax*xder*deltax*xder*(1.0+lambda);
			Hessian.at<double>(2,3) += deltax*xder*deltax*yder;
			Hessian.at<double>(2,4) += deltax*xder*deltay*xder;
			Hessian.at<double>(2,5) += deltax*xder*deltay*yder;
			Hessian.at<double>(2,6) += deltax*xder*deltax*deltay*xder;
			Hessian.at<double>(2,7) += deltax*xder*deltax*deltay*yder;
			Hessian.at<double>(2,8) += deltax*xder*deltax*deltax*xder;
			Hessian.at<double>(2,9) += deltax*xder*deltax*deltax*yder;
			Hessian.at<double>(2,10) += deltax*xder*deltay*deltay*xder;
			Hessian.at<double>(2,11) += deltax*xder*deltay*deltay*yder;

			Hessian.at<double>(3,0) += deltax*yder*xder;
			Hessian.at<double>(3,1) += deltax*yder*yder;
			Hessian.at<double>(3,2) += deltax*yder*deltax*xder;
			Hessian.at<double>(3,3) += deltax*yder*deltax*yder*(1.0+lambda);
			Hessian.at<double>(3,4) += deltax*yder*deltay*xder;
			Hessian.at<double>(3,5) += deltax*yder*deltay*yder;
			Hessian.at<double>(3,6) += deltax*yder*deltax*deltay*xder;
			Hessian.at<double>(3,7) += deltax*yder*deltax*deltay*yder;
			Hessian.at<double>(3,8) += deltax*yder*deltax*deltax*xder;
			Hessian.at<double>(3,9) += deltax*yder*deltax*deltax*yder;
			Hessian.at<double>(3,10) += deltax*yder*deltay*deltay*xder;
			Hessian.at<double>(3,11) += deltax*yder*deltay*deltay*yder;

			Hessian.at<double>(4,0) += deltay*xder*xder;
			Hessian.at<double>(4,1) += deltay*xder*yder;
			Hessian.at<double>(4,2) += deltay*xder*deltax*xder;
			Hessian.at<double>(4,3) += deltay*xder*deltax*yder;
			Hessian.at<double>(4,4) += deltay*xder*deltay*xder*(1.0+lambda);
			Hessian.at<double>(4,5) += deltay*xder*deltay*yder;
			Hessian.at<double>(4,6) += deltay*xder*deltax*deltay*xder;
			Hessian.at<double>(4,7) += deltay*xder*deltax*deltay*yder;
			Hessian.at<double>(4,8) += deltay*xder*deltax*deltax*xder;
			Hessian.at<double>(4,9) += deltay*xder*deltax*deltax*yder;
			Hessian.at<double>(4,10) += deltay*xder*deltay*deltay*xder;
			Hessian.at<double>(4,11) += deltay*xder*deltay*deltay*yder;

			Hessian.at<double>(5,0) += deltay*yder*xder;
			Hessian.at<double>(5,1) += deltay*yder*yder;
			Hessian.at<double>(5,2) += deltay*yder*deltax*xder;
			Hessian.at<double>(5,3) += deltay*yder*deltax*yder;
			Hessian.at<double>(5,4) += deltay*yder*deltay*xder;
			Hessian.at<double>(5,5) += deltay*yder*deltay*yder*(1.0+lambda);
			Hessian.at<double>(5,6) += deltay*yder*deltax*deltay*xder;
			Hessian.at<double>(5,7) += deltay*yder*deltax*deltay*yder;
			Hessian.at<double>(5,8) += deltay*yder*deltax*deltax*xder;
			Hessian.at<double>(5,9) += deltay*yder*deltax*deltax*yder;
			Hessian.at<double>(5,10) += deltay*yder*deltay*deltay*xder;
			Hessian.at<double>(5,11) += deltay*yder*deltay*deltay*yder;

			Hessian.at<double>(6,0) += deltax*deltay*xder*xder;
			Hessian.at<double>(6,1) += deltax*deltay*xder*yder;
			Hessian.at<double>(6,2) += deltax*deltay*xder*deltax*xder;
			Hessian.at<double>(6,3) += deltax*deltay*xder*deltax*yder;
			Hessian.at<double>(6,4) += deltax*deltay*xder*deltay*xder;
			Hessian.at<double>(6,5) += deltax*deltay*xder*deltay*yder;
			Hessian.at<double>(6,6) += deltax*deltay*xder*deltax*deltay*xder*(1.0+lambda);
			Hessian.at<double>(6,7) += deltax*deltay*xder*deltax*deltay*yder;
			Hessian.at<double>(6,8) += deltax*deltay*xder*deltax*deltax*xder;
			Hessian.at<double>(6,9) += deltax*deltay*xder*deltax*deltax*yder;
			Hessian.at<double>(6,10) += deltax*deltay*xder*deltay*deltay*xder;
			Hessian.at<double>(6,11) += deltax*deltay*xder*deltay*deltay*yder;

			Hessian.at<double>(7,0) += deltax*deltay*yder*xder;
			Hessian.at<double>(7,1) += deltax*deltay*yder*yder;
			Hessian.at<double>(7,2) += deltax*deltay*yder*deltax*xder;
			Hessian.at<double>(7,3) += deltax*deltay*yder*deltax*yder;
			Hessian.at<double>(7,4) += deltax*deltay*yder*deltay*xder;
			Hessian.at<double>(7,5) += deltax*deltay*yder*deltay*yder;
			Hessian.at<double>(7,6) += deltax*deltay*yder*deltax*deltay*xder;
			Hessian.at<double>(7,7) += deltax*deltay*yder*deltax*deltay*yder*(1.0+lambda);
			Hessian.at<double>(7,8) += deltax*deltay*yder*deltax*deltax*xder;
			Hessian.at<double>(7,9) += deltax*deltay*yder*deltax*deltax*yder;
			Hessian.at<double>(7,10) += deltax*deltay*yder*deltay*deltay*xder;
			Hessian.at<double>(7,11) += deltax*deltay*yder*deltay*deltay*yder;

			Hessian.at<double>(8,0) += deltax*deltax*xder*xder;
			Hessian.at<double>(8,1) += deltax*deltax*xder*yder;
			Hessian.at<double>(8,2) += deltax*deltax*xder*deltax*xder;
			Hessian.at<double>(8,3) += deltax*deltax*xder*deltax*yder;
			Hessian.at<double>(8,4) += deltax*deltax*xder*deltay*xder;
			Hessian.at<double>(8,5) += deltax*deltax*xder*deltay*yder;
			Hessian.at<double>(8,6) += deltax*deltax*xder*deltax*deltay*xder;
			Hessian.at<double>(8,7) += deltax*deltax*xder*deltax*deltay*yder;
			Hessian.at<double>(8,8) += deltax*deltax*xder*deltax*deltax*xder*(1.0+lambda);
			Hessian.at<double>(8,9) += deltax*deltax*xder*deltax*deltax*yder;
			Hessian.at<double>(8,10) += deltax*deltax*xder*deltay*deltay*xder;
			Hessian.at<double>(8,11) += deltax*deltax*xder*deltay*deltay*yder;

			Hessian.at<double>(9,0) += deltax*deltax*yder*xder;
			Hessian.at<double>(9,1) += deltax*deltax*yder*yder;
			Hessian.at<double>(9,2) += deltax*deltax*yder*deltax*xder;
			Hessian.at<double>(9,3) += deltax*deltax*yder*deltax*yder;
			Hessian.at<double>(9,4) += deltax*deltax*yder*deltay*xder;
			Hessian.at<double>(9,5) += deltax*deltax*yder*deltay*yder;
			Hessian.at<double>(9,6) += deltax*deltax*yder*deltax*deltay*xder;
			Hessian.at<double>(9,7) += deltax*deltax*yder*deltax*deltay*yder;
			Hessian.at<double>(9,8) += deltax*deltax*yder*deltax*deltax*xder;
			Hessian.at<double>(9,9) += deltax*deltax*yder*deltax*deltax*yder*(1.0+lambda);
			Hessian.at<double>(9,10) += deltax*deltax*yder*deltay*deltay*xder;
			Hessian.at<double>(9,11) += deltax*deltax*yder*deltay*deltay*yder;

			Hessian.at<double>(10,0) += deltay*deltay*xder*xder;
			Hessian.at<double>(10,1) += deltay*deltay*xder*yder;
			Hessian.at<double>(10,2) += deltay*deltay*xder*deltax*xder;
			Hessian.at<double>(10,3) += deltay*deltay*xder*deltax*yder;
			Hessian.at<double>(10,4) += deltay*deltay*xder*deltay*xder;
			Hessian.at<double>(10,5) += deltay*deltay*xder*deltay*yder;
			Hessian.at<double>(10,6) += deltay*deltay*xder*deltax*deltay*xder;
			Hessian.at<double>(10,7) += deltay*deltay*xder*deltax*deltay*yder;
			Hessian.at<double>(10,8) += deltay*deltay*xder*deltax*deltax*xder;
			Hessian.at<double>(10,9) += deltay*deltay*xder*deltax*deltax*yder;
			Hessian.at<double>(10,10) += deltay*deltay*xder*deltay*deltay*xder*(1.0+lambda);
			Hessian.at<double>(10,11) += deltay*deltay*xder*deltay*deltay*yder;

			Hessian.at<double>(11,0) += deltay*deltay*yder*xder;
			Hessian.at<double>(11,1) += deltay*deltay*yder*yder;
			Hessian.at<double>(11,2) += deltay*deltay*yder*deltax*xder;
			Hessian.at<double>(11,3) += deltay*deltay*yder*deltax*yder;
			Hessian.at<double>(11,4) += deltay*deltay*yder*deltay*xder;
			Hessian.at<double>(11,5) += deltay*deltay*yder*deltay*yder;
			Hessian.at<double>(11,6) += deltay*deltay*yder*deltax*deltay*xder;
			Hessian.at<double>(11,7) += deltay*deltay*yder*deltax*deltay*yder;
			Hessian.at<double>(11,8) += deltay*deltay*yder*deltax*deltax*xder;
			Hessian.at<double>(11,9) += deltay*deltay*yder*deltax*deltax*yder;
			Hessian.at<double>(11,10) += deltay*deltay*yder*deltay*deltay*xder;
			Hessian.at<double>(11,11) += deltay*deltay*yder*deltay*deltay*yder*(1.0+lambda);

			/// Jacobian Matrix [dg/dU0; dg/dV0; dg/dUx; dg/dVx; dg/dUy; dg/dVy; dg/dUxy; dg/dVxy; dg/dUxx; dg/dVxx; dg/dUyy; dg/dVyy]
			double C_components = (f_values[k]-fm)/sqrt(sum_f_minus_fm_squared)-(g_values[k]-gm)/sqrt(sum_g_minus_gm_squared);
			k++;
			Jacobian.at<double>(0,0) += xder*C_components;
			Jacobian.at<double>(1,0) += yder*C_components;
			Jacobian.at<double>(2,0) += deltax*xder*C_components;
			Jacobian.at<double>(3,0) += deltax*yder*C_components;
			Jacobian.at<double>(4,0) += deltay*xder*C_components;
			Jacobian.at<double>(5,0) += deltay*yder*C_components;
			Jacobian.at<double>(6,0) += deltax*deltay*xder*C_components;
			Jacobian.at<double>(7,0) += deltax*deltay*yder*C_components;
			Jacobian.at<double>(8,0) += deltax*deltax*xder*C_components;
			Jacobian.at<double>(9,0) += deltax*deltax*yder*C_components;
			Jacobian.at<double>(10,0) += deltay*deltay*xder*C_components;
			Jacobian.at<double>(11,0) += deltay*deltay*yder*C_components;
		}
	}
	Hessian *= 1.0/sqrt(sum_g_minus_gm_squared);	
}
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_quadratic_LM_optim(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold)
{
        double abs_tolerance = 1;
        double rel_tolerance = 1;
        double max_val_nu = 1e7;
		double lambda = 1e-10;
        double nu = 2.0;
        unsigned int max_iterations = 1e4;
        unsigned int iterations = 0;
        double U0 = P0[0];
        double V0 = P0[1];
        double Ux = P0[2];
        double Vx = P0[3];
        double Uy = P0[4];
        double Vy = P0[5];
        double Uxy = P0[6];
        double Vxy = P0[7];
        double Uxx = P0[8];
        double Vxx = P0[9];
        double Uyy = P0[10];
        double Vyy = P0[11];

        unsigned int Indexi = SubsetLength/2 + i * GridLength;
        unsigned int Indexj = SubsetLength/2 + j * GridLength;

		//auto tr1 = std::chrono::high_resolution_clock::now();
		double fm, sum_f_minus_fm_squared;
        //fm = get_fm(img, Indexi,  Indexj,  SubsetLength);
        //sum_f_minus_fm_squared = get_sum_f_minus_fm_squared(img, Indexi,  Indexj,  SubsetLength, fm);
		std::vector<double> f_values = get_f_values(img, Indexi, Indexj, SubsetLength);
		std::tie(fm, sum_f_minus_fm_squared) = get_gm_and_sum_g_minus_gm_squared_short(f_values);
		
		double gm, sum_g_minus_gm_squared;
		std::vector<double> g_values = get_g_values(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
		std::tie(gm, sum_g_minus_gm_squared) = get_gm_and_sum_g_minus_gm_squared_short(g_values);
		//gm = get_gm(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
		//sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
		bool data_uniform = (is_data_uniform(sum_f_minus_fm_squared)||is_data_uniform(sum_g_minus_gm_squared));
		//auto tr2= std::chrono::high_resolution_clock::now();
		
        //double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
        double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD_optim_short(fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared, f_values, g_values);
		
        Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;
		
		cv::Mat Hessian(12, 12,CV_64F, cv::Scalar(0));
		cv::Mat Jacobian(12,1, CV_64F, cv::Scalar(0));
		//calculate_Hessian_Jacobian(Hessian, Jacobian, img, fptr_img1, P0, SplineDegree, SubsetLength, Indexi, Indexj, 3, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared, lambda);
		//auto tr3 = std::chrono::high_resolution_clock::now();
		calculate_Hessian_Jacobian_short(Hessian, Jacobian, img, fptr_img1, P0, SplineDegree, SubsetLength, Indexi, Indexj, 3, f_values, fm, sum_f_minus_fm_squared, g_values, gm, sum_g_minus_gm_squared, lambda);
		//auto tr4= std::chrono::high_resolution_clock::now();
		//std::cout << "fg took" << std::chrono::duration_cast<std::chrono::milliseconds>(tr2-tr1).count() << " milliseconds, HJ took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr4-tr3).count() << " milliseconds" << std::endl;
		//std::cout << "fg took " << std::chrono::duration_cast<std::chrono::microseconds>(tr2-tr1).count() << " microseconds, HJ took " << std::chrono::duration_cast<std::chrono::microseconds>(tr4-tr3).count() << " microseconds" << std::endl;
        while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold) && data_uniform == 0)
        {
            iterations++;
            cv::Mat X(12,1,CV_64F);
            cv::solve(Hessian, Jacobian, X, cv::DECOMP_CHOLESKY);

            std::vector<double> Suggested_Solution = {U0+X.at<double>(0), V0+X.at<double>(1), Ux+X.at<double>(2), Vx+X.at<double>(3), Uy+X.at<double>(4), Vy+X.at<double>(5),  Uxy+X.at<double>(6), Vxy+X.at<double>(7), Uxx+X.at<double>(8), Vxx+X.at<double>(9), Uyy+X.at<double>(10), Vyy+X.at<double>(11)};
			double suggested_gm, suggested_sum_g_minus_gm_squared;
			//suggested_gm = get_gm(fptr_img1, Suggested_Solution, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
			//suggested_sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, suggested_gm, img.cols, img.rows);
			//double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, suggested_gm, suggested_sum_g_minus_gm_squared);
			std::vector<double> suggested_g_values = get_g_values(fptr_img1,  Suggested_Solution, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
			std::tie(suggested_gm, suggested_sum_g_minus_gm_squared) = get_gm_and_sum_g_minus_gm_squared_short(suggested_g_values);
			double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD_optim_short(fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared, f_values, suggested_g_values); 
            Suggested_Correlation_Coefficient = 1.0-0.5*Suggested_Correlation_Coefficient;

            if (Suggested_Correlation_Coefficient > Correlation_Coefficient_old)
            {
                // Improvement
                Correlation_Coefficient_old = Suggested_Correlation_Coefficient;
                abs_tolerance = cv::sum(cv::abs(X)).val[0];
                rel_tolerance = cv::sum(cv::abs(X)).val[0]/(abs(P0[0])+abs(P0[1])+abs(P0[2])+abs(P0[3])+abs(P0[4])+abs(P0[5])+abs(P0[6])+abs(P0[7])+abs(P0[8])+abs(P0[9])+abs(P0[10])+abs(P0[11]));

                U0 += X.at<double>(0);
                V0 += X.at<double>(1);
                Ux += X.at<double>(2);
                Vx += X.at<double>(3);
                Uy += X.at<double>(4);
                Vy += X.at<double>(5);
                Uxy += X.at<double>(6);
                Vxy += X.at<double>(7);
                Uxx += X.at<double>(8);
                Vxx += X.at<double>(9);
                Uyy += X.at<double>(10);
                Vyy += X.at<double>(11);
				g_values = suggested_g_values;
				gm = suggested_gm;
				sum_g_minus_gm_squared = suggested_sum_g_minus_gm_squared;
				data_uniform = is_data_uniform(sum_g_minus_gm_squared);
				
                // Decrease lambda => More Gauss-Newton, less Gradient Search
                nu = 2;
                lambda /= 3.0;
				
				// Recompute Hessian and Jacobian
				//calculate_Hessian_Jacobian(Hessian, Jacobian, img, fptr_img1, Suggested_Solution, SplineDegree, SubsetLength, Indexi, Indexj, 3, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared, lambda);
				calculate_Hessian_Jacobian_short(Hessian, Jacobian, img, fptr_img1, Suggested_Solution, SplineDegree, SubsetLength, Indexi, Indexj, 3, f_values, fm, sum_f_minus_fm_squared, g_values, gm, sum_g_minus_gm_squared, lambda);
				
				/// calculate_Hessian_Jacobian amd get_gm contain the same call for InterpolatedValue 
				/// Should do this once -> get InterpolatedValue (in vector) and use this in get_gm and calculate_Hessian_Jacobian
            }
            else
            {
                // No Improvement
                // Increase lambda => Less Gauss-Newton, more Gradient Search
                double lambda_new = lambda*nu;
                nu *= 2.0;
				
				// Scale Diagonal of Hessian
				// Jacobian stays the same
				for (unsigned int i = 0; i < 12; i++)
				{
					Hessian.at<double>(i,i) *= (1.0+lambda_new)/(1.0+lambda);
				}
				lambda = lambda_new;
            }
        }

        std::vector<double> returnvector;
        returnvector.push_back(U0);
        returnvector.push_back(V0);
        returnvector.push_back(Ux);
        returnvector.push_back(Vx);
        returnvector.push_back(Uy);
        returnvector.push_back(Vy);
        returnvector.push_back(Uxy);
        returnvector.push_back(Vxy);
        returnvector.push_back(Uxx);
        returnvector.push_back(Vxx);
        returnvector.push_back(Uyy);
        returnvector.push_back(Vyy);
		if (data_uniform)
			returnvector.push_back(-1);
		else
		{
			returnvector.push_back(Correlation_Coefficient_old);
		}
        

        return returnvector;
}
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_irregular_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold)
{
        double abs_tolerance = 1;
        double rel_tolerance = 1;
        double max_val_nu = 1e7;
        unsigned int max_iterations = 1e4;
        unsigned int iterations = 0;
        double U0 = P0[0];
        double V0 = P0[1];
        double Ux = P0[2];
        double Vx = P0[3];
        double Uy = P0[4];
        double Vy = P0[5];
        double Uxy = P0[6];
        double Vxy = P0[7];

        double lambda = 1e-10;
        double nu = 2.0;
        //double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength);
        //Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;
		
		unsigned int Indexi = SubsetLength/2 + i * GridLength;
		unsigned int Indexj = SubsetLength/2 + j * GridLength;        
		double fm = get_fm(img, Indexi,  Indexj,  SubsetLength);
        double sum_f_minus_fm_squared = get_sum_f_minus_fm_squared(img, Indexi,  Indexj,  SubsetLength, fm);
		double gm = get_gm(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
		double sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
		bool data_uniform = (is_data_uniform(sum_f_minus_fm_squared)||is_data_uniform(sum_g_minus_gm_squared));
        double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
        Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;

        while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold) && data_uniform == 0)
        {
            iterations++;			
            cv::Mat Hessian(8,8,CV_64F, cv::Scalar(0));
            cv::Mat Jacobian(8,1, CV_64F, cv::Scalar(0));
            for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
            {
				double deltax = ((double)ii-(double)Indexi);
                for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
                {
                    double deltay = ((double)jj-(double)Indexj);
                    double xder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay, SplineDegree, 0);
                    double yder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay, SplineDegree, 1);

                    /// Hessian Matrix [dg/dU0 dg/dU0, dg/dU0 dg/dV0, dg/dU0 dg/dUx, dg/dU0 dg/dVx, dg/dU0 dg/Uy, dg/dU0 dg/Vy, dg/dU0 dg/dUxy, dg/dU0 dg/dVxy;
                    ///                 dg/dV0 dg/dU0, dg/dV0 dg/dV0, dg/dV0 dg/dUx, dg/dV0 dg/dVx, dg/dV0 dg/Uy, dg/dV0 dg/Vy, dg/dV0 dg/dUxy, dg/dV0 dg/dVxy;
                    ///                 dg/dUx dg/dU0, dg/dUx dg/dV0, dg/dUx dg/dUx, dg/dUx dg/dVx, dg/dUx dg/Uy, dg/dUx dg/Vy, dg/dUx dg/dUxy, dg/dUx dg/dVxy;
                    ///                 dg/dVx dg/dU0, dg/dVx dg/dV0, dg/dVx dg/dUx, dg/dVx dg/dVx, dg/dVx dg/Uy, dg/dVx dg/Vy, dg/dVx dg/dUxy, dg/dVx dg/dVxy;
                    ///                 dg/dUy dg/dU0, dg/dUy dg/dV0, dg/dUy dg/dUx, dg/dUy dg/dVx, dg/dUy dg/Uy, dg/dUy dg/Vy, dg/dUy dg/dUxy, dg/dUy dg/dVxy;
                    ///                 dg/dVy dg/dU0, dg/dVy dg/dV0, dg/dVy dg/dUx, dg/dVy dg/dVx, dg/dVy dg/Uy, dg/dVy dg/Vy, dg/dVy dg/dUxy, dg/dVy dg/dVxy;
                    ///                 dg/dUxy dg/dU0, dg/dUxy dg/dV0, dg/dUxy dg/dUx, dg/dUxy dg/dVx, dg/dUxy dg/Uy, dg/dUxy dg/Vy, dg/dUxy dg/dUxy, dg/dUxy dg/dVxy;
                    ///                 dg/dVxy dg/dU0, dg/dVxy dg/dV0, dg/dVxy dg/dUx, dg/dVxy dg/dVx, dg/dVxy dg/Uy, dg/dVxy dg/Vy, dg/dVxy dg/dUxy, dg/dVxy dg/dVxy]
                    Hessian.at<double>(0,0) += xder*xder*(1.0+lambda);
                    Hessian.at<double>(0,1) += xder*yder;
                    Hessian.at<double>(0,2) += xder*deltax*xder;
                    Hessian.at<double>(0,3) += xder*deltax*yder;
                    Hessian.at<double>(0,4) += xder*deltay*xder;
                    Hessian.at<double>(0,5) += xder*deltay*yder;
                    Hessian.at<double>(0,6) += xder*deltax*deltay*xder;
                    Hessian.at<double>(0,7) += xder*deltax*deltay*yder;

                    Hessian.at<double>(1,0) += yder*xder;
                    Hessian.at<double>(1,1) += yder*yder*(1.0+lambda);
                    Hessian.at<double>(1,2) += yder*deltax*xder;
                    Hessian.at<double>(1,3) += yder*deltax*yder;
                    Hessian.at<double>(1,4) += yder*deltay*xder;
                    Hessian.at<double>(1,5) += yder*deltay*yder;
                    Hessian.at<double>(1,6) += yder*deltax*deltay*xder;
                    Hessian.at<double>(1,7) += yder*deltax*deltay*yder;

                    Hessian.at<double>(2,0) += deltax*xder*xder;
                    Hessian.at<double>(2,1) += deltax*xder*yder;
                    Hessian.at<double>(2,2) += deltax*xder*deltax*xder*(1.0+lambda);
                    Hessian.at<double>(2,3) += deltax*xder*deltax*yder;
                    Hessian.at<double>(2,4) += deltax*xder*deltay*xder;
                    Hessian.at<double>(2,5) += deltax*xder*deltay*yder;
                    Hessian.at<double>(2,6) += deltax*xder*deltax*deltay*xder;
                    Hessian.at<double>(2,7) += deltax*xder*deltax*deltay*yder;

                    Hessian.at<double>(3,0) += deltax*yder*xder;
                    Hessian.at<double>(3,1) += deltax*yder*yder;
                    Hessian.at<double>(3,2) += deltax*yder*deltax*xder;
                    Hessian.at<double>(3,3) += deltax*yder*deltax*yder*(1.0+lambda);
                    Hessian.at<double>(3,4) += deltax*yder*deltay*xder;
                    Hessian.at<double>(3,5) += deltax*yder*deltay*yder;
                    Hessian.at<double>(3,6) += deltax*yder*deltax*deltay*xder;
                    Hessian.at<double>(3,7) += deltax*yder*deltax*deltay*yder;

                    Hessian.at<double>(4,0) += deltay*xder*xder;
                    Hessian.at<double>(4,1) += deltay*xder*yder;
                    Hessian.at<double>(4,2) += deltay*xder*deltax*xder;
                    Hessian.at<double>(4,3) += deltay*xder*deltax*yder;
                    Hessian.at<double>(4,4) += deltay*xder*deltay*xder*(1.0+lambda);
                    Hessian.at<double>(4,5) += deltay*xder*deltay*yder;
                    Hessian.at<double>(4,6) += deltay*xder*deltax*deltay*xder;
                    Hessian.at<double>(4,7) += deltay*xder*deltax*deltay*yder;

                    Hessian.at<double>(5,0) += deltay*yder*xder;
                    Hessian.at<double>(5,1) += deltay*yder*yder;
                    Hessian.at<double>(5,2) += deltay*yder*deltax*xder;
                    Hessian.at<double>(5,3) += deltay*yder*deltax*yder;
                    Hessian.at<double>(5,4) += deltay*yder*deltay*xder;
                    Hessian.at<double>(5,5) += deltay*yder*deltay*yder*(1.0+lambda);
                    Hessian.at<double>(5,6) += deltay*yder*deltax*deltay*xder;
                    Hessian.at<double>(5,7) += deltay*yder*deltax*deltay*yder;

                    Hessian.at<double>(6,0) += deltax*deltay*xder*xder;
                    Hessian.at<double>(6,1) += deltax*deltay*xder*yder;
                    Hessian.at<double>(6,2) += deltax*deltay*xder*deltax*xder;
                    Hessian.at<double>(6,3) += deltax*deltay*xder*deltax*yder;
                    Hessian.at<double>(6,4) += deltax*deltay*xder*deltay*xder;
                    Hessian.at<double>(6,5) += deltax*deltay*xder*deltay*yder;
                    Hessian.at<double>(6,6) += deltax*deltay*xder*deltax*deltay*xder*(1.0+lambda);
                    Hessian.at<double>(6,7) += deltax*deltay*xder*deltax*deltay*yder;

                    Hessian.at<double>(7,0) += deltax*deltay*yder*xder;
                    Hessian.at<double>(7,1) += deltax*deltay*yder*yder;
                    Hessian.at<double>(7,2) += deltax*deltay*yder*deltax*xder;
                    Hessian.at<double>(7,3) += deltax*deltay*yder*deltax*yder;
                    Hessian.at<double>(7,4) += deltax*deltay*yder*deltay*xder;
                    Hessian.at<double>(7,5) += deltax*deltay*yder*deltay*yder;
                    Hessian.at<double>(7,6) += deltax*deltay*yder*deltax*deltay*xder;
                    Hessian.at<double>(7,7) += deltax*deltay*yder*deltax*deltay*yder*(1.0+lambda);
                    /// Jacobian Matrix [dg/dU0; dg/dV0; dg/dUx; dg/dVx; dg/dUy; dg/dVy; dg/dUxy; dg/dVxy]
                    double C_components = (img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared);
                    Jacobian.at<double>(0,0) += xder*C_components;
                    Jacobian.at<double>(1,0) += yder*C_components;
                    Jacobian.at<double>(2,0) += deltax*xder*C_components;
                    Jacobian.at<double>(3,0) += deltax*yder*C_components;
                    Jacobian.at<double>(4,0) += deltay*xder*C_components;
                    Jacobian.at<double>(5,0) += deltay*yder*C_components;
                    Jacobian.at<double>(6,0) += deltax*deltay*xder*C_components;
                    Jacobian.at<double>(7,0) += deltax*deltay*yder*C_components;
                }
            }
            /// Suggestion: Multiply both side with sqrt(sum_g_minus_gm_squared)/2
            Jacobian *= 2/sqrt(sum_g_minus_gm_squared);
            Hessian *= 2/sum_g_minus_gm_squared;
            cv::Mat X(8,1,CV_64F);
            cv::solve(Hessian, Jacobian, X, cv::DECOMP_CHOLESKY);

            std::vector<double> Suggested_Solution = {U0+X.at<double>(0), V0+X.at<double>(1), Ux+X.at<double>(2), Vx+X.at<double>(3), Uy+X.at<double>(4), Vy+X.at<double>(5),  Uxy+X.at<double>(6), Vxy+X.at<double>(7)};
            //double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength);
			double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
            Suggested_Correlation_Coefficient = 1.0-0.5*Suggested_Correlation_Coefficient;

            if (Suggested_Correlation_Coefficient > Correlation_Coefficient_old)
            {
                // Improvement
                Correlation_Coefficient_old = Suggested_Correlation_Coefficient;
                abs_tolerance = cv::sum(cv::abs(X)).val[0];
                rel_tolerance = cv::sum(cv::abs(X)).val[0]/(abs(P0[0])+abs(P0[1])+abs(P0[2])+abs(P0[3])+abs(P0[4])+abs(P0[5])+abs(P0[6])+abs(P0[7]));

                U0 += X.at<double>(0);
                V0 += X.at<double>(1);
                Ux += X.at<double>(2);
                Vx += X.at<double>(3);
                Uy += X.at<double>(4);
                Vy += X.at<double>(5);
                Uxy += X.at<double>(6);
                Vxy += X.at<double>(7);
				
				gm = get_gm(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
				sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
				data_uniform = is_data_uniform(sum_g_minus_gm_squared);
				if (data_uniform) break;
							
                // Decrease lambda => More Gauss-Newton, less Gradient Search
                nu = 2;
                lambda /= 3.0;
            }
            else
            {
                // No Improvement
                // Increase lambda => Less Gauss-Newton, more Gradient Search
                lambda *= nu;
                nu *= 2.0;
            }
        }
        if (iterations>max_iterations)
            std::cout << "Max Iterations Reached" << std::endl;

        std::vector<double> returnvector;
        returnvector.push_back(U0);
        returnvector.push_back(V0);
        returnvector.push_back(Ux);
        returnvector.push_back(Vx);
        returnvector.push_back(Uy);
        returnvector.push_back(Vy);
        returnvector.push_back(Uxy);
        returnvector.push_back(Vxy);
        // Return empty second order terms
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
		double Correlation_Coefficient;
		if (data_uniform)
			Correlation_Coefficient = -1;
		else
		{
			Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, returnvector, SplineDegree, SubsetLength, GridLength);
			Correlation_Coefficient = 1-0.5*Correlation_Coefficient;
		}
        returnvector.push_back(Correlation_Coefficient);

        return returnvector;
}
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_affine_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold)
{
        double abs_tolerance = 1;
        double rel_tolerance = 1;
        double max_val_nu = 1e7;
        unsigned int max_iterations = 1e4;
        unsigned int iterations = 0;
        double U0 = P0[0];
        double V0 = P0[1];
        double Ux = P0[2];
        double Vx = P0[3];
        double Uy = P0[4];
        double Vy = P0[5];

        double lambda = 1e-10;
        double nu = 2.0;
        //double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength);
        //Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;
		
		unsigned int Indexi = SubsetLength/2 + i * GridLength;
		unsigned int Indexj = SubsetLength/2 + j * GridLength;        
		double fm = get_fm(img, Indexi,  Indexj,  SubsetLength);
        double sum_f_minus_fm_squared = get_sum_f_minus_fm_squared(img, Indexi,  Indexj,  SubsetLength, fm);
		double gm = get_gm(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
		double sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
		bool data_uniform = (is_data_uniform(sum_f_minus_fm_squared)||is_data_uniform(sum_g_minus_gm_squared));
        double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
        Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;

        while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold) && data_uniform == 0)
        {
            iterations++;
			data_uniform = is_data_uniform(sum_g_minus_gm_squared);
			if (data_uniform) break;
			
            cv::Mat Hessian(6,6,CV_64F, cv::Scalar(0));
            cv::Mat Jacobian(6,1, CV_64F, cv::Scalar(0));
            for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
            {
				double deltax = ((double)ii-(double)Indexi);
                for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
                {
                    double deltay = ((double)jj-(double)Indexj);
                    double xder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay, jj+V0+Vx*deltax+Vy*deltay, SplineDegree, 0);
                    double yder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay, jj+V0+Vx*deltax+Vy*deltay, SplineDegree, 1);

                    /// Hessian Matrix [dg/dU0 dg/dU0, dg/dU0 dg/dV0, dg/dU0 dg/dUx, dg/dU0 dg/dVx, dg/dU0 dg/Uy, dg/dU0 dg/Vy;
                    ///                 dg/dV0 dg/dU0, dg/dV0 dg/dV0, dg/dV0 dg/dUx, dg/dV0 dg/dVx, dg/dV0 dg/Uy, dg/dV0 dg/Vy;
                    ///                 dg/dUx dg/dU0, dg/dUx dg/dV0, dg/dUx dg/dUx, dg/dUx dg/dVx, dg/dUx dg/Uy, dg/dUx dg/Vy;
                    ///                 dg/dVx dg/dU0, dg/dVx dg/dV0, dg/dVx dg/dUx, dg/dVx dg/dVx, dg/dVx dg/Uy, dg/dVx dg/Vy;
                    ///                 dg/dUy dg/dU0, dg/dUy dg/dV0, dg/dUy dg/dUx, dg/dUy dg/dVx, dg/dUy dg/Uy, dg/dUy dg/Vy;
                    ///                 dg/dVy dg/dU0, dg/dVy dg/dV0, dg/dVy dg/dUx, dg/dVy dg/dVx, dg/dVy dg/Uy, dg/dVy dg/Vy;
                    Hessian.at<double>(0,0) += xder*xder*(1.0+lambda);
                    Hessian.at<double>(0,1) += xder*yder;
                    Hessian.at<double>(0,2) += xder*deltax*xder;
                    Hessian.at<double>(0,3) += xder*deltax*yder;
                    Hessian.at<double>(0,4) += xder*deltay*xder;
                    Hessian.at<double>(0,5) += xder*deltay*yder;

                    Hessian.at<double>(1,0) += yder*xder;
                    Hessian.at<double>(1,1) += yder*yder*(1.0+lambda);
                    Hessian.at<double>(1,2) += yder*deltax*xder;
                    Hessian.at<double>(1,3) += yder*deltax*yder;
                    Hessian.at<double>(1,4) += yder*deltay*xder;
                    Hessian.at<double>(1,5) += yder*deltay*yder;

                    Hessian.at<double>(2,0) += deltax*xder*xder;
                    Hessian.at<double>(2,1) += deltax*xder*yder;
                    Hessian.at<double>(2,2) += deltax*xder*deltax*xder*(1.0+lambda);
                    Hessian.at<double>(2,3) += deltax*xder*deltax*yder;
                    Hessian.at<double>(2,4) += deltax*xder*deltay*xder;
                    Hessian.at<double>(2,5) += deltax*xder*deltay*yder;

                    Hessian.at<double>(3,0) += deltax*yder*xder;
                    Hessian.at<double>(3,1) += deltax*yder*yder;
                    Hessian.at<double>(3,2) += deltax*yder*deltax*xder;
                    Hessian.at<double>(3,3) += deltax*yder*deltax*yder*(1.0+lambda);
                    Hessian.at<double>(3,4) += deltax*yder*deltay*xder;
                    Hessian.at<double>(3,5) += deltax*yder*deltay*yder;

                    Hessian.at<double>(4,0) += deltay*xder*xder;
                    Hessian.at<double>(4,1) += deltay*xder*yder;
                    Hessian.at<double>(4,2) += deltay*xder*deltax*xder;
                    Hessian.at<double>(4,3) += deltay*xder*deltax*yder;
                    Hessian.at<double>(4,4) += deltay*xder*deltay*xder*(1.0+lambda);
                    Hessian.at<double>(4,5) += deltay*xder*deltay*yder;

                    Hessian.at<double>(5,0) += deltay*yder*xder;
                    Hessian.at<double>(5,1) += deltay*yder*yder;
                    Hessian.at<double>(5,2) += deltay*yder*deltax*xder;
                    Hessian.at<double>(5,3) += deltay*yder*deltax*yder;
                    Hessian.at<double>(5,4) += deltay*yder*deltay*xder;
                    Hessian.at<double>(5,5) += deltay*yder*deltay*yder*(1.0+lambda);
                    /// Jacobian Matrix [dg/dU0; dg/dV0; dg/dUx; dg/dVx; dg/dUy; dg/dVy]
                    double C_components = (img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay, jj+V0+Vx*deltax+Vy*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared);
                    Jacobian.at<double>(0,0) += xder*C_components;
                    Jacobian.at<double>(1,0) += yder*C_components;
                    Jacobian.at<double>(2,0) += deltax*xder*C_components;
                    Jacobian.at<double>(3,0) += deltax*yder*C_components;
                    Jacobian.at<double>(4,0) += deltay*xder*C_components;
                    Jacobian.at<double>(5,0) += deltay*yder*C_components;
                }
            }
            /// Suggestion: Multiply both side with sqrt(sum_g_minus_gm_squared)/2
            Jacobian *= 2/sqrt(sum_g_minus_gm_squared);
            Hessian *= 2/sum_g_minus_gm_squared;
            cv::Mat X(6,1,CV_64F);
            cv::solve(Hessian, Jacobian, X, cv::DECOMP_CHOLESKY);

            std::vector<double> Suggested_Solution = {U0+X.at<double>(0), V0+X.at<double>(1), Ux+X.at<double>(2), Vx+X.at<double>(3), Uy+X.at<double>(4), Vy+X.at<double>(5)};
            //double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength);
			double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
            Suggested_Correlation_Coefficient = 1.0-0.5*Suggested_Correlation_Coefficient;

            if (Suggested_Correlation_Coefficient > Correlation_Coefficient_old)
            {
                // Improvement
                Correlation_Coefficient_old = Suggested_Correlation_Coefficient;
                abs_tolerance = cv::sum(cv::abs(X)).val[0];
                rel_tolerance = cv::sum(cv::abs(X)).val[0]/(abs(P0[0])+abs(P0[1])+abs(P0[2])+abs(P0[3])+abs(P0[4])+abs(P0[5]));
                U0 += X.at<double>(0);
                V0 += X.at<double>(1);
                Ux += X.at<double>(2);
                Vx += X.at<double>(3);
                Uy += X.at<double>(4);
                Vy += X.at<double>(5);
				
				gm = get_gm(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
				sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
				data_uniform = is_data_uniform(sum_g_minus_gm_squared);
				if (data_uniform) break;
			
                // Decrease lambda => More Gauss-Newton, less Gradient Search
                nu = 2;
                lambda /= 3.0;
            }
            else
            {
                // No Improvement
                // Increase lambda => Less Gauss-Newton, more Gradient Search
                lambda *= nu;
                nu *= 2.0;
            }
        }
        if (iterations>max_iterations)
            std::cout << "Max Iterations Reached" << std::endl;

        std::vector<double> returnvector;
        returnvector.push_back(U0);
        returnvector.push_back(V0);
        returnvector.push_back(Ux);
        returnvector.push_back(Vx);
        returnvector.push_back(Uy);
        returnvector.push_back(Vy);
        // return empty second order terms
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
		double Correlation_Coefficient;
		if (data_uniform)
			Correlation_Coefficient = -1;
		else
		{
			Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, returnvector, SplineDegree, SubsetLength, GridLength);
			Correlation_Coefficient = 1-0.5*Correlation_Coefficient;
		}
        returnvector.push_back(Correlation_Coefficient);

        return returnvector;
}
/*--------------------------------------------------------------------------*/
std::vector<double> iteration_rigid_LM(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold)
{
        double abs_tolerance = 1;
        double rel_tolerance = 1;
        double max_val_nu = 1e7;
        unsigned int max_iterations = 1e4;
        unsigned int iterations = 0;
        double U0 = P0[0];
        double V0 = P0[1];

        double lambda = 1e-10;
        double nu = 2.0;
        //double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength);

		unsigned int Indexi = SubsetLength/2 + i * GridLength;
		unsigned int Indexj = SubsetLength/2 + j * GridLength;
		
        double fm = get_fm(img, Indexi,  Indexj,  SubsetLength);
        double sum_f_minus_fm_squared = get_sum_f_minus_fm_squared(img, Indexi,  Indexj,  SubsetLength, fm);
		double gm = get_gm(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
		double sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, P0, Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
		bool data_uniform = (is_data_uniform(sum_f_minus_fm_squared)||is_data_uniform(sum_g_minus_gm_squared));

        double Correlation_Coefficient_old = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, P0, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
        Correlation_Coefficient_old = 1.0-0.5*Correlation_Coefficient_old;		
		
        while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold) && data_uniform == 0)
        {
            iterations++;			
            cv::Mat Hessian(2,2,CV_64F, cv::Scalar(0));
            cv::Mat Jacobian(2,1, CV_64F, cv::Scalar(0));
			for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
			{
				for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
				{
					double xder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0, jj+V0, SplineDegree, 0);
					double yder = getDerivativeValue(fptr_img1, img.cols, img.rows, ii+U0, jj+V0, SplineDegree, 1);
					/// Hessian Matrix [dg/dU0 dg/dU0, dg/dU0 dg/dV0; dg/dV0 dg/dU0, dg/dV0 dg/dV0]
					Hessian.at<double>(0,0) += xder*xder*(1.0+lambda);
					Hessian.at<double>(0,1) += xder*yder;
					Hessian.at<double>(1,0) += yder*xder;
					Hessian.at<double>(1,1) += yder*yder*(1.0+lambda);
					/// Jacobian Matrix [dg/dU0; dg/dV0]
					double C_components = (img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0, jj+V0, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared);
					Jacobian.at<double>(0,0) += xder*C_components;
					Jacobian.at<double>(1,0) += yder*C_components;
				}
			}
			/// Suggestion: Multiply both side with sqrt(sum_g_minus_gm_squared)/2
			Jacobian *= 2/sqrt(sum_g_minus_gm_squared);
			Hessian *= 2/sum_g_minus_gm_squared;
			cv::Mat X(2,1,CV_64F);
			cv::solve(Hessian, Jacobian, X, cv::DECOMP_CHOLESKY);
			
			std::vector<double> Suggested_Solution = {U0+X.at<double>(0), V0+X.at<double>(1)};
			//double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength);
			double Suggested_Correlation_Coefficient = Correlation_Coefficient_ZNSSD_optim(img, fptr_img1, i, j, Suggested_Solution, SplineDegree, SubsetLength, GridLength, fm, sum_f_minus_fm_squared, gm, sum_g_minus_gm_squared);
			Suggested_Correlation_Coefficient = 1.0-0.5*Suggested_Correlation_Coefficient;
            if (Suggested_Correlation_Coefficient > Correlation_Coefficient_old)
            {
                // Improvement
                Correlation_Coefficient_old = Suggested_Correlation_Coefficient;
                abs_tolerance = cv::sum(cv::abs(X)).val[0];
                rel_tolerance = cv::sum(cv::abs(X)).val[0]/(abs(P0[0])+abs(P0[1]));
                U0 += X.at<double>(0);
                V0 += X.at<double>(1);

				gm = get_gm(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
				sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, Suggested_Solution , Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);
				data_uniform = is_data_uniform(sum_g_minus_gm_squared);
				if (data_uniform) break;				
				
                // Decrease lambda => More Gauss-Newton, less Gradient Search
                nu = 2;
                lambda /= 3.0;
            }
            else
            {
                // No Improvement
                // Increase lambda => Less Gauss-Newton, more Gradient Search
                lambda *= nu;
                nu *= 2.0;
            }
        }
        std::vector<double> returnvector;
        returnvector.push_back(U0);
        returnvector.push_back(V0);
        // return empty first and second order terms
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
        returnvector.push_back(0);
		double Correlation_Coefficient;
		if (data_uniform)
		{
			Correlation_Coefficient = -1;
		}
		else
		{
			Correlation_Coefficient = Correlation_Coefficient_ZNSSD(img, fptr_img1, i, j, returnvector, SplineDegree, SubsetLength, GridLength);
			Correlation_Coefficient = 1-0.5*Correlation_Coefficient;
		}
        returnvector.push_back(Correlation_Coefficient);

        return returnvector;
}
/*--------------------------------------------------------------------------*/
double Correlation_Coefficient_ZNSSD(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength)
{
    double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }

    double Correlation_Coefficient = 0;
    unsigned int Indexi = SubsetLength/2 + i * GridLength;
    unsigned int Indexj = SubsetLength/2 + j * GridLength;/*
    double sum_g_minus_gm_squared = 0;
    double gm = 0;
    double fm = 0;
    double sum_f_minus_fm_squared = 0;
    for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
    {
        double deltax = ((double)ii-(double)Indexi);
        for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
        {
            double deltay = ((double)jj-(double)Indexj);
            fm += (double)img.at<float>(jj,ii);
            gm += InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree);
        }
    }
    gm /= (2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1);
    fm /= (2*(SubsetLength/2)+1)*(2*(SubsetLength/2)+1);
    for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
    {
        double deltax = ((double)ii-(double)Indexi);
        for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
        {
            double deltay = ((double)jj-(double)Indexj);
            double g_minus_gm = InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm;
            sum_g_minus_gm_squared += g_minus_gm*g_minus_gm;
            sum_f_minus_fm_squared += ((double)img.at<float>(jj,ii)-fm)*((double)img.at<float>(jj,ii)-fm);
        }
    }
	*/
	double fm = get_fm(img, Indexi,  Indexj,  SubsetLength);
	double sum_f_minus_fm_squared = get_sum_f_minus_fm_squared(img, Indexi,  Indexj,  SubsetLength, fm);
	double gm = get_gm(fptr_img1, P, Indexi, Indexj, SubsetLength, SplineDegree, img.cols, img.rows);
	double sum_g_minus_gm_squared = get_sum_g_minus_gm_squared(fptr_img1, P, Indexi, Indexj, SubsetLength, SplineDegree, gm, img.cols, img.rows);

    for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
    {
        double deltax = ((double)ii-(double)Indexi);
        for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
        {
            double deltay = ((double)jj-(double)Indexj);
            double C_value = ((img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared));
            Correlation_Coefficient += C_value*C_value;
        }
    }
    if (isnan(Correlation_Coefficient))
    {
        Correlation_Coefficient = 4;
	}
    return Correlation_Coefficient;
}
double Correlation_Coefficient_ZNSSD_optim(const cv::Mat &img, float *fptr_img1, const unsigned int &i, const unsigned int &j, const std::vector<double> &P, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &sum_g_minus_gm_squared)
{
    double U0 = P[0];
    double V0 = P[1];
    double Ux = 0;
    double Vx = 0;
    double Uy = 0;
    double Vy = 0;
    double Uxy = 0;
    double Vxy = 0;
    double Uxx = 0;
    double Vxx = 0;
    double Uyy = 0;
    double Vyy = 0;
    if (P.size() > 2)
    {
        Ux = P[2];
        Vx = P[3];
        Uy = P[4];
        Vy = P[5];
    }
    if (P.size() > 6)
    {
        Uxy = P[6];
        Vxy = P[7];
    }
    if (P.size() > 8)
    {
        Uxx = P[8];
        Vxx = P[9];
        Uyy = P[10];
        Vyy = P[11];
    }

    double Correlation_Coefficient = 0;
    unsigned int Indexi = SubsetLength/2 + i * GridLength;
    unsigned int Indexj = SubsetLength/2 + j * GridLength;
    for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
    {
        double deltax = ((double)ii-(double)Indexi);
        for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
        {
            double deltay = ((double)jj-(double)Indexj);
            double C_value = ((img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared));
            Correlation_Coefficient += C_value*C_value;
        }
    }
    if (isnan(Correlation_Coefficient))
    {
        Correlation_Coefficient = 4;
	}
    return Correlation_Coefficient;
}
/*--------------------------------------------------------------------------*/
double Correlation_Coefficient_ZNSSD_optim_short(const double &fm, const double &sum_f_minus_fm_squared, const double &gm, const double &sum_g_minus_gm_squared, const std::vector<double> &f_values, const std::vector<double> &g_values)
{

    double Correlation_Coefficient = 0;
	//for(auto it = g_values.begin(); it != g_values.end(); ++it)
	//{
		
	//}
	
	for (unsigned int i = 0; i <f_values.size(); i++)
	{
		double C_value = (f_values[i]-fm)/sqrt(sum_f_minus_fm_squared)-(g_values[i]-gm)/sqrt(sum_g_minus_gm_squared);
		Correlation_Coefficient += C_value*C_value;
	}
    if (isnan(Correlation_Coefficient))
    {
        Correlation_Coefficient = 4;
	}
    return Correlation_Coefficient;
}
/*--------------------------------------------------------------------------*/
/*double Correlation_Coefficient_ZNSSD_optim_short(const cv::Mat &img, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &fm, const double &sum_f_minus_fm_squared, const std::vector<double> &g_values, const double &gm, const double &sum_g_minus_gm_squared)
{
    double Correlation_Coefficient = 0;
    for (unsigned int ii = Indexi-SubsetLength/2; ii < Indexi+SubsetLength/2+1; ii++)
    {
        double deltax = ((double)ii-(double)Indexi);
        for (unsigned int jj = Indexj-SubsetLength/2; jj < Indexj+SubsetLength/2+1; jj++)
        {
            double deltay = ((double)jj-(double)Indexj);
            double C_value = ((img.at<float>(jj,ii)-fm)/sqrt(sum_f_minus_fm_squared)-(InterpolatedValue(fptr_img1, img.cols, img.rows, ii+U0+Ux*deltax+Uy*deltay+Uxy*deltax*deltay+Uxx*deltax*deltax+Uyy*deltay*deltay, jj+V0+Vx*deltax+Vy*deltay+Vxy*deltax*deltay+Vxx*deltax*deltax+Vyy*deltay*deltay, SplineDegree)-gm)/sqrt(sum_g_minus_gm_squared));
            Correlation_Coefficient += C_value*C_value;
        }
    }
    if (isnan(Correlation_Coefficient))
    {
        Correlation_Coefficient = 4;
	}
    return Correlation_Coefficient;
}*/
/*--------------------------------------------------------------------------*/
static double getDerivativeValue(float *fptr_img1, const unsigned int &cols, const unsigned int &rows, const double &x, const double &y, const unsigned int &SplineDegree, const unsigned int &direction)
{
    return (double)InterpolatedValueDerivative(fptr_img1, cols, rows, x+0.5*(1.0-(double)direction), y+0.5*(double)direction, SplineDegree-1*(1-direction), SplineDegree-1*direction)
                -(double)InterpolatedValueDerivative(fptr_img1, cols, rows, x-0.5*(1.0-(double)direction), y-0.5*(double)direction, SplineDegree-1*(1-direction), SplineDegree-1*direction);
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Point> get_valid_Neighbours(const cv::Mat &M_valid_points, const cv::Mat &Computed_Points, const unsigned int &x, const unsigned int &y, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned &offset)
{
    unsigned int xtotal = x*GridLength + SubsetLength/2+offset;
    unsigned int ytotal = y*GridLength + SubsetLength/2+offset;

    std::vector<cv::Point> valid_neighbours;
    if ((int)M_valid_points.at<uchar>(ytotal, xtotal-GridLength))
    {
        // left neighbour
        if ((int)Computed_Points.at<uchar>(y,x-1)==0)
        {
            // if not computed
            cv::Point Neighbour;
            Neighbour.x = ((xtotal-GridLength)-SubsetLength/2-offset)/GridLength;
            Neighbour.y = (ytotal-SubsetLength/2-offset)/GridLength;
            valid_neighbours.push_back(Neighbour);
        }
    }
    if ((int)M_valid_points.at<uchar>(ytotal, xtotal+GridLength))
    {
        // right neighbour
        if ((int)Computed_Points.at<uchar>(y,x+1)==0)
        {
            // if not computed
            cv::Point Neighbour;
            Neighbour.x = ((xtotal+GridLength)-SubsetLength/2-offset)/GridLength;
            Neighbour.y = (ytotal-SubsetLength/2-offset)/GridLength;
            valid_neighbours.push_back(Neighbour);
        }
    }
    if ((int)M_valid_points.at<uchar>(ytotal-GridLength, xtotal))
    {
        // upper neighbour
        if ((int)Computed_Points.at<uchar>(y-1,x)==0)
        {
            // if not computed
            cv::Point Neighbour;
            Neighbour.x = (xtotal-SubsetLength/2-offset)/GridLength;
            Neighbour.y = ((ytotal-GridLength)-SubsetLength/2-offset)/GridLength;
            valid_neighbours.push_back(Neighbour);
        }
    }
    if ((int)M_valid_points.at<uchar>(ytotal+GridLength, xtotal))
    {
        // lower neighbour
        if ((int)Computed_Points.at<uchar>(y+1,x)==0)
        {
            // if not computed
            cv::Point Neighbour;
            Neighbour.x = (xtotal-SubsetLength/2-offset)/GridLength;
            Neighbour.y = ((ytotal+GridLength)-SubsetLength/2-offset)/GridLength;
            valid_neighbours.push_back(Neighbour);
        }
    }
    return valid_neighbours;
}
/*--------------------------------------------------------------------------*/
