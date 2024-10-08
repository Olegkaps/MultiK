// PriorC MultiK


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "../Utils/MyTypes.h"
#include "../Utils/MyUtils.h"
#include "PriorCTypes.h"
#include "PriorCUtils.h"
#include "PriorCInteractions.h"
#include "PriorCProbs.h"
#include "PriorCReading.h"
#include "PriorCArgumentParcing.h"


FILE* logfile;


float* splineX;
float* splineXinit;
float* splineY;
float* splineYinit;
var residual;
var* outliers;
var* outliersdist;
size_t outliersdistLength;
float* FDRX;
float* FDRXinit;
float* FDRY; 
float* FDRYinit; 


int main(int argc, char *argv[])
{
    numbers nums;
    nums.binLength = 0;

    priorCArgs args;    
    
    args = parse(args, argc, argv);
    nums.noOfBins = args.noOfBins;
    nums.binCount = args.binCount;
    nums.binLength = args.binLength;
    nums.distUpThres = args.distUpThres;
    nums.distLowThres = args.distLowThres;
    nums.noOfFrags = args.noOfFrags;

    float distScaling = 1000000.0 / nums.binLength;
    if(nums.binCount <= (nums.distUpThres - nums.distLowThres))
    {
        nums.mainDicSize = nums.binCount - nums.distLowThres;
    } else
    {
        nums.mainDicSize = nums.distUpThres - nums.distLowThres + 1;
    }

    var* mainDic = xcalloc(nums.mainDicSize, sizeof(var)); 
    interactions inf = read_Interactions(mainDic, args.contactCountsFileName, outliers, nums);
    nums = inf.nums;
    mainDic = inf.mainDic;
    var** binStats = calloc2d(nums.noOfBins, 5, sizeof(var));
    // binStats
    // 0: lower boundary of distances in this bin
    // 1: upper boundary of distances in this bin
    // 2: no. of possible pairs w/in this range of distances
    // 3: sumoverallContactCounts
    // 4: Sumoveralldistances in this bin in distScaling vals
    binStats = makeBinsFromInteractions(binStats, mainDic, outliersdist, outliersdistLength, nums);
    fragReturn res = generate_FragPairs(binStats, args.fragsfileName, nums);
    binStats = res.binStats;
    nums = res.nums;


    float* x = xcalloc(nums.noOfBins, sizeof(float));
    float* y = xcalloc(nums.noOfBins, sizeof(float));


    float** probTuple = xmalloc(2 * sizeof(float*));
    probTuple[0] = x;
    probTuple[1] = y;
    //probTuple[2] = yerr;
    probTuple = calculateProbabilities(probTuple, mainDic, binStats, nums);
    
    //outliers = SortedList();
    //outliersdist = SortedList();
    probTuple = fit_Spline(mainDic, probTuple[0], probTuple[1], 1, nums);

    for(var i = 1; i < args.noOfPasses; i++)
    {
        /*memset(mainDic, 0, nums.mainDicSize * sizeof(var));
        memset2d_f(probTuple, 2, 0, nums.noOfBins);
        inf = read_Interactions(mainDic, args.contactCountsFileName, outliers, nums);
        nums = inf.nums;
        mainDic = inf.mainDic;
        memset2d_v(binStats, nums.noOfBins, 0, 5);  
        
        binStats = makeBinsFromInteractions(binStats, mainDic, outliersdist, outliersdistLength, nums);
        res = generate_FragPairs(binStats, args.fragsfileName, nums);
        binStats = res.binStats;
        nums = res.nums;
        
        probTuple = calculateProbabilities(probTuple, mainDic, binStats, nums);*/
        probTuple = fit_Spline(mainDic, probTuple[0], probTuple[1],  i, nums);
    }


    splineX = probTuple[0];
    splineY = probTuple[1];

    write_spline(args, splineX, splineY, nums);

    free(mainDic);
    free2d(nums.noOfBins, binStats);
    free2d(2, probTuple);
    fclose(logfile);
    return 0;
}