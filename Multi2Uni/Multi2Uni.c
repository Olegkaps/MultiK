#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <omp.h>

#include "../Utils/MyTypes.h"
#include "../Utils/MyUtils.h"
#include "Multi2UniTypes.h"
#include "EMUtils.h"
#include "Multi2UniArgumentParcing.h"
#include "Multi2UniReader.h"



int main(int argc, char *argv[])
{
    DATA data;
    
    data = parseArgs(data, argc, argv);
    data = read_spline_prior(data);

    CSR matrix;

    FILE* multiFile = xopen(data.pathes.multiPath, "r");
    //var readsNum, binPairsNum, valuesNum, TotalReadsNum;
    xscanf(5, multiFile, "Bin Length: %u\nTotal reads: %u\nNum of reads: %u\nNum of bin pairs: %u\nNum of items: %u\n#bins\n", &data.spline.binLength2, &data.nums.TotalReadsNum, &data.nums.readsNum, &data.nums.binPairsNum, &data.nums.valuesNum);
    
    matrix = init_matrix(matrix, data);
    matrix = readBinsStrings(multiFile, matrix, data);

    xscanf(1, multiFile, "#reads %u\n", &data.nums.readsPosition);

    matrix = readReadsStrings(multiFile, matrix, data);
    
    fclose(multiFile);
    
    if(data.spline.binLength1 != data.spline.binLength2)
    {
        printf("WARNING: you use files with different bin size!\n");
    }
    // do not use probs and distances anymore
    xfree(data.spline.probs);
    xfree(data.spline.distances);
        

#ifndef NOOMP
    omp_set_dynamic(0);
    var thread_num;        
    var num_threads = omp_get_num_procs();
    float *thread_diffs = xmalloc(num_threads*sizeof(float)); 
    float valuesForThread;
#endif
    var maxIter = 500;
    float diffThre = 0.0001; // 0.01*0.01
    //var selectChangeThre = 1; // DO NOT USE
    float diff = 0.0;
    float pi_new, pi, d;
    clock_t start, end;
    double cpu_time_used;
    var ind;

#ifndef NOOMP
#pragma omp target
#endif
for(var k = 0; k < maxIter; k++)
{
    diff = 0.0;
    zeroMem(matrix.Z.sumBinPair, data.nums.binPairsNum);
#ifndef NOOMP
    zeroMem(thread_diffs, num_threads);
#endif
    cpu_time_used = 0.0;
    start = clock();

    // Z sum
#ifndef NOOMP
    #pragma omp parallel for
#endif
    for(var i = 0; i < data.nums.binPairsNum; i++)
    {
        for(var j = matrix.Z.binPairs_ptr[i]; j < matrix.Z.binPairs_ptr[i+1]; j++)
        {
            matrix.Z.sumBinPair[i] += matrix.Pi.SumID[matrix.Z.pairID[j]];
        }
        //matrix.Z.sumBinPair[i] = 1 / matrix.Z.sumBinPair[i]; // RIGHT CALCULATION
    }

    // Update Pi
#ifndef NOOMP    
    #pragma omp parallel for private(ind, pi, pi_new, d, thread_num)
#endif
    for(var i = 0; i < data.nums.readsNum; i++)
    {
#ifndef NOOMP
        thread_num = omp_get_thread_num();
#endif
        matrix.Pi.SumID[i] = 0.0;
        for(var j = matrix.Pi.pairID_ptr[i]; j < matrix.Pi.pairID_ptr[i+1]; j++)
        {
            ind = matrix.Pi.binPairs[j];
            pi = matrix.Pi.pi_values[j];
            pi_new = pi * matrix.Z.sumBinPair[ind] + matrix.Z.const_to_pi[ind];
            
            d = pi - pi_new;
#ifndef NOOMP
            thread_diffs[thread_num] += d*d;
#else
            diff += d*d;
#endif
            matrix.Pi.pi_values[j] = pi_new;
            matrix.Pi.SumID[i] += matrix.Pi.pi_values[j];
        }
        matrix.Pi.SumID[i] = 1/matrix.Pi.SumID[i];
    }
    
#ifndef NOOMP
    for(int i = 0; i < num_threads; i++)
    {
        diff += thread_diffs[i];
    }
#endif


    end = clock();
    cpu_time_used += ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Iteration Number %u, time used: %f sec, difference: %f.\n", k, cpu_time_used, diff);

    // Pi difference
    if(diff < diffThre)
    {
        break;
    }
}

    // calculate probabilities and write out the result
    FILE* multiReadsFile = xopen(data.pathes.multiPath, "r");
    FILE* outFile = xopen(data.pathes.outputPath, "w");
    fseek(multiReadsFile, data.nums.readsPosition, SEEK_SET);
    
    readAndWriteReads(multiReadsFile, outFile, matrix, data);

    fclose(outFile);
    fclose(multiReadsFile);

    free_matrix(matrix);
   

#ifdef NOOMP
    printf("Multi2Uni NOOMP Done (^_^)\n");
#else
    printf("Multi2Uni Done (^_^)\n");
#endif


    return 0;
}