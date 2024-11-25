extern fragReturn generate_FragPairs(var** binStats, char* fragsfileName, numbers nums);

extern double** calculateProbabilities(double** probTuple, var* mainDic, var** binStats, numbers nums);

extern double** fit_Spline(var* mainDic, double* x, double* y, var passNo, numbers nums);
