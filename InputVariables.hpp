#ifndef H_InputVariables
#define H_InputVariables
struct InputVariables {
	
	std::string path;
	unsigned int SplineDegree;
	unsigned int SubsetLength;
	unsigned int GridLength;
	unsigned int ShapeFunction;
	unsigned int propagationfunction;
	unsigned int ordering;
	unsigned int xStart;
	unsigned int xEnd;
	unsigned int yStart;
	unsigned int yEnd;
	unsigned int offset;
	unsigned int Number_Of_Threads;
	unsigned int MaxPixelYVertical;
	double abs_tolerance_threshold;
	double rel_tolerance_threshold;
	double minimum_corrcoeff_IG;
	
	unsigned int xStart_ROI;
	unsigned int yStart_ROI;
	unsigned int horx_ROI;
	unsigned int very_ROI;
	
	unsigned int DICNeeded;
	unsigned int CalibrationNeeded;
	unsigned int CalculateRefractionIndex;
};
#endif