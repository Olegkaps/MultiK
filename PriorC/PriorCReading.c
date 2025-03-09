#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "MyTypes.h"
#include "MyUtils.h"
#include "PriorCTypes.h"
#include "PriorCReading.h"
#include "hashTable.h"

interactions* read_Interactions(var* mainDic, char* contactCountsFileName, var* outliers, numbers* nums, priorCArgs* args)
{
    
    FILE* contactCountsFile = xopen(contactCountsFileName, "r");
    // THE VERY FIRST LINE OF THE contactCountsFile INCLUDE TECHNICAL INFORMATION
    // !IMPORTANT
    //xscanf(2, contactCountsFile, "Number of bins: %u\nLength of a single bin: %u", &nums.binCount, &nums.binLength);
    
    
    fprintf(stdout, "Reading the contact counts file to generate bins...\n");
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();


    nums->observedInterAllSum = 0;
    nums->observedIntraAllSum = 0;
    nums->observedIntraInRangeSum = 0; 
    
    
    char* ch1 = xmalloc(30*sizeof(char));
    char* ch2 = xmalloc(30*sizeof(char));
    var bin1, bin2, contactCount, distance;

    while(!feof(contactCountsFile)){
        // !IMPORTANT
        xscanf(6, contactCountsFile,"%s\t%u\t%s\t%u\t%u\t%u\n", ch1, &bin1, ch2, &bin2, &distance, &contactCount);

        if ((strcmp(args->dnaDensityFileName, ""))) {
            
            DNA_bins* dnas = find_bins_of_chr(args->density_table, ch2);
            dnas->Count[bin2] += contactCount;
        }
        // TO DO: OUTLIERS //
        // !IMPORTANT
        // TO DO: OUTLIERS //

        if(strcmp(ch1, ch2)) // INTERACTION TYPE - INTER 
        {
            nums->observedInterAllSum += contactCount;
        } 
        else // INTERACTION TYPE - INTRA
        {
            nums->observedIntraAllSum += contactCount;
            if((nums->distLowThres <= distance) && (distance <= nums->distUpThres)) // DISTANCE IN RANGE
            {
                mainDic[distance - nums->distLowThres] += contactCount;
                nums->observedIntraInRangeSum += contactCount;
            }
        }
    }
    fclose(contactCountsFile);


    end = clock();
    cpu_time_used += ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, "Interactions file read. Time took %f sec\n\n", cpu_time_used);
    

    // OUTPUT
    interactions* inf = xmalloc(sizeof(interactions));
    inf->nums = nums;
    inf->mainDic = mainDic;
    xfree(ch1);
    xfree(ch2);

    return inf;
}


void write_spline(priorCArgs* args, double* x, double* y, numbers* nums)
{
    FILE* splineOut = xopen(args->outputPath, "w");

    fprintf(splineOut, "Bin Length: %u\n", nums->binLength);
    fprintf(splineOut, "Num of Bins: %u\n", nums->noOfBins);
    double y_sum = 0.0;
    
    
    for(var i = 0; i < nums->noOfBins; i++)
    {
        if(y[i] > 0.0000000000000099999)
        {
            fprintf(splineOut, "%.1f\t%.9f\n", x[i], y[i]);
        } 
    }
    fclose(splineOut);
    
    fprintf(stdout, "=========================\n");
    fprintf(stdout, "PriorC completed successfully\n\n\n");

}
