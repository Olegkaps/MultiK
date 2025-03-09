void expectationMaximization(DATA* data, CSR* matrix) {
    #ifndef NOOMP
    omp_set_dynamic(0);
    var thread_num;        
    var num_threads = omp_get_num_procs();
    float *thread_diffs = xmalloc(num_threads*sizeof(float)); 
    float valuesForThread;
#endif
    var maxIter = 500;
    float diffThre = 0.000001; // 0.001*0.001
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
    zeroMem(matrix->Z.sumBinPair, data->nums.binPairsNum);
#ifndef NOOMP
    zeroMem(thread_diffs, num_threads);
#endif
    cpu_time_used = 0.0;
    start = clock();

    // Z sum
#ifndef NOOMP
    #pragma omp parallel for
#endif
    for(var i = 0; i < data->nums.binPairsNum; i++)
    {
        matrix->Z.sumBinPair[i] = 0.0;
        for(var j = matrix->Z.binPairs_ptr[i]; j < matrix->Z.binPairs_ptr[i+1]; j++)
        {
            matrix->Z.sumBinPair[i] += matrix->Pi.SumID[matrix->Z.pairID[j]];
        }
    }

    // Update Pi
#ifndef NOOMP    
    #pragma omp parallel for private(ind, pi, pi_new, d, thread_num)
#endif
    for(var i = 0; i < data->nums.readsNum; i++)
    {
#ifndef NOOMP
        thread_num = omp_get_thread_num();
#endif
        matrix->Pi.SumID[i] = 0.0;
        for(var j = matrix->Pi.pairID_ptr[i]; j < matrix->Pi.pairID_ptr[i+1]; j++)
        {
            ind = matrix->Pi.binPairs[j];
            pi = matrix->Pi.pi_values[j];
            pi_new = pi * matrix->Z.sumBinPair[ind] + matrix->Z.const_to_pi[ind];
            
            d = (pi - pi_new) / (pi + pi_new);
#ifndef NOOMP
            thread_diffs[thread_num] += d*d;
#else
            diff += d*d;
#endif
            matrix->Pi.pi_values[j] = pi_new;
            matrix->Pi.SumID[i] += matrix->Pi.pi_values[j];
        }
        matrix->Pi.SumID[i] = 1/matrix->Pi.SumID[i];
    }
    
#ifndef NOOMP
    for(int i = 0; i < num_threads; i++)
    {
        diff += thread_diffs[i];
    }
#endif


    end = clock();
    cpu_time_used += ((double) (end - start)) / CLOCKS_PER_SEC;
    diff /= data->nums.valuesNum;
    printf("Iteration Number %u, time used: %f sec, average difference: %f.\n", k, cpu_time_used, diff);

    // Pi difference
    if(diff < diffThre)
    {
        break;
    }
}

}