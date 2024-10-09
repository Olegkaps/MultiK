#include "alglib-cpp/src/interpolation.h"
#include "alglib-cpp/src/stdafx.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

using namespace alglib;

void SmoothingSpline(double* x, double* y, unsigned int N, double w)
{
    try
    {
        real_1d_array x_new;
        real_1d_array y_new;
        spline1dinterpolant s;
        spline1dfitreport rep;
	x_new.setcontent(N, x);
	y_new.setcontent(N, y);


        spline1dfit(x_new, y_new, 50, w, s, rep);

        for(int i = 0; i < N; i++) {
                y[i] = spline1dcalc(s, x[i]);
        }


    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
    }
}



#ifdef __cplusplus
}
#endif
