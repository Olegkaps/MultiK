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


typedef struct {
    var* pairID_ptr;
    double* SumID;
    var* binPairs;
    double* pi_values;
} CSR_Pi;

typedef struct {
    var* binPairs_ptr;
    double* sumBinPair;
    double* const_to_pi;
    var* pairID;
    //double* z_values;
} CSR_Z;

typedef struct {
    CSR_Pi Pi;
    CSR_Z Z;
} CSR;

typedef struct {
    var* distances;
    double* probs;
    double interProb;
    var splineLength;
    var binLength1;
    var binLength2;
} splineData;

typedef struct {
    char *priorPath;
    char *multiPath;
    char *outputPath;
    char *outdir;
} pathesData;

typedef struct {
    var binPairsNum;
    var readsNum;
    var valuesNum;
    var TotalReadsNum;
    var readsPosition;
} numbersData;


typedef struct {
    splineData spline;
    pathesData pathes;
    numbersData nums;
} DATA;

#include "../Multi2Uni/Multi2UniArgumentParcing.h"



CSR init_matrix(CSR matrix, DATA data)
{
    matrix.Pi.pairID_ptr = xmalloc((data.nums.readsNum+1)*sizeof(var));
    matrix.Pi.pairID_ptr[0] = 0;
    matrix.Pi.SumID = xmalloc(data.nums.readsNum*sizeof(double));
    matrix.Pi.binPairs = xmalloc(data.nums.valuesNum*sizeof(var));
    matrix.Pi.pi_values = xmalloc(data.nums.valuesNum*sizeof(double));
    matrix.Z.binPairs_ptr = xmalloc((data.nums.binPairsNum+1)*sizeof(var));
    matrix.Z.binPairs_ptr[0] = 0;
    matrix.Z.sumBinPair = xmalloc(data.nums.binPairsNum*sizeof(double));
    matrix.Z.const_to_pi = xmalloc(data.nums.binPairsNum*sizeof(double));
    matrix.Z.pairID = xmalloc(data.nums.valuesNum*sizeof(var));
    //Z.z_values = malloc(valuesNum*sizeof(double));)

    return matrix;
}


int main(int argc, char *argv[])
{
    DATA data;
    data.nums.binPairsNum = 4;
    data.nums.readsNum = 5;
    data.nums.valuesNum = 15;


    CSR matrix;
    matrix =  init_matrix(matrix, data);

    var ptr[] = {0, 3, 5, 8, 11, 15};
    memcpy(matrix.Pi.pairID_ptr, ptr, 6*sizeof(var));
    double ptr1[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    memcpy(matrix.Pi.SumID, ptr1, 5*sizeof(double));
    var ptr2[] = {0, 2, 3, 1, 3, 0, 1, 2, 1, 2, 3, 0, 1, 2, 3};
    memcpy(matrix.Pi.binPairs, ptr2, 15*sizeof(var));

    var ptr3[] = {0, 3, 7, 11, 15};
    memcpy(matrix.Z.binPairs_ptr, ptr3, 5*sizeof(var));
    double ptr10[] = {0.0, 0.0, 0.0, 0.0};
    memcpy(matrix.Z.sumBinPair, ptr10, 4*sizeof(double));
    double ptr11[] = {100.0, 200.0, 333.0, 1.0};
    memcpy(matrix.Z.const_to_pi, ptr11, 4*sizeof(double));
    var ptr9[] = {0, 2, 4, 1, 2, 3, 4, 0, 2, 3, 4, 0, 1, 3, 4};
    memcpy(matrix.Z.pairID, ptr9, 15*sizeof(var));

    var tests_num = 5;
     
    double t1[] = {0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.9, 0.1, 0.1, 0.9, 0.9};
    double t2[] = {0.9, 0.9, 0.1, 0.9, 0.1, 0.9, 0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.9, 0.1, 0.1};
    double t3[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double t4[] = {0.3, 0.3, 0.4, 0.6, 0.4, 0.3, 0.3, 0.4, 0.2, 0.4, 0.4, 0.2, 0.2, 0.3, 0.3};
    double t5[] = {0.57, 0.23, 0.94, 0.36, 0.894, 0.23, 0.53, 0.94, 0.52, 0.34, 0.84, 0.22, 0.62, 0.23, 0.73};
    
    double* tests[] = {t1, t2, t3, t4, t5};    
    
    FILE* test_res = xopen("test_em.txt", "w");
    for(var ggplot=0; ggplot < tests_num; ggplot++){
    memcpy(matrix.Pi.pi_values, (double*)tests[ggplot], 15*sizeof(double));
    //matrix.Pi.pi_values = tests[ggplot];


    var ind;
    for(var i = 0; i < data.nums.readsNum; i++)
    {
        matrix.Pi.SumID[i] = 0.0;
        for(var j = matrix.Pi.pairID_ptr[i]; j < matrix.Pi.pairID_ptr[i+1]; j++)
        {
            matrix.Pi.SumID[i] += matrix.Pi.pi_values[j];
        }
        matrix.Pi.SumID[i] = 1/matrix.Pi.SumID[i];
    }

#ifndef NOOMP
    omp_set_dynamic(0);
    var thread_num;        
    var num_threads = omp_get_num_procs();
    double *thread_diffs = xmalloc(num_threads*sizeof(double)); 
    double valuesForThread;
#endif
    var maxIter = 500;
    double diffThre = 0.0001; // 0.01*0.01
    //var selectChangeThre = 1; // DO NOT USE
    double diff = 0.0;
    double pi_new, pi, d;
    clock_t start, end;
    double cpu_time_used;

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
        //printf("%f\n", matrix.Z.sumBinPair[i]);
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
    double resultProb;
    //write down
    for(var i = 0; i < data.nums.readsNum; i++)
    {
        for(var j = matrix.Pi.pairID_ptr[i]; j < matrix.Pi.pairID_ptr[i+1]; j++)
        {
            resultProb = matrix.Pi.pi_values[j] * matrix.Pi.SumID[i];
            fprintf(test_res, "%.9f ", resultProb);
        }
    }
    fprintf(test_res, "\n");
}
    fclose(test_res);
    //free_matrix(matrix);
   

#ifdef NOOMP
    printf("Multi2Uni NOOMP Done (^_^)\n");
#else
    printf("Multi2Uni Done (^_^)\n");
#endif


    return 0;
}