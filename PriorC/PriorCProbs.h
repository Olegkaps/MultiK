fragReturn generate_FragPairs(var** binStats, char* fragsfileName, numbers nums)
{
    fprintf(stdout, "Reading fragments file...\n");
    FILE* fragsfile = xopen(fragsfileName, "r");
    
    
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();

    var possibleIntraInRangeCountPerChr;    
    var possibleIntraAllCount = 0;
    nums.possibleInterAllCount = 0;
    nums.possibleIntraInRangeCount = 0;
    nums.interChrProb = 0.0;
    nums.baselineIntraChrProb = 0.0;


    char chrName[10];
    var chrLength;
    
    var up, down, npairs;
    float intxnDistance;
    int c;
    //###########################
    //mapsThres not accounted for
    //###########################
    if(nums.binLength != 0)
    {
        //xscanf(1, fragsfile, "Number of fragments: %u\n", &nums.noOfFrags);                                                                            
        while(!feof(fragsfile))                                                 
        {                                                                       
            xscanf(2, fragsfile, "%s\t%u\n", chrName, &chrLength);                    
            chrLength /= nums.binLength;                                                
            nums.possibleInterAllCount += chrLength * (nums.noOfFrags - chrLength);        
            possibleIntraAllCount += (chrLength * (chrLength + 1)) / 2;          
            for(var i = 0; i < nums.noOfBins; i++)                                   
            {
                if(binStats[i][3] == 0)                                          
                {                                                                
                    nums.noOfBins = i;                                            
                    break;                                                       
                }                                                                
                for(var intxnDistance = binStats[i][0]; intxnDistance <= binStats[i][1]; intxnDistance++)
                {
                    if (intxnDistance < nums.distLowThres)
                    {
                        continue;
                    }
                    if (intxnDistance > nums.distUpThres || intxnDistance > chrLength)
                    {
                        break;
                    }
                    npairs = (var)chrLength - intxnDistance;
                    possibleIntraInRangeCountPerChr += npairs;
                    binStats[i][2] += npairs;
                    binStats[i][4] += intxnDistance * npairs;
                }

            }
            nums.possibleIntraInRangeCount += possibleIntraInRangeCountPerChr;
            possibleIntraInRangeCountPerChr = 0;
        }
        nums.interChrProb = ((nums.possibleInterAllCount != 0) ? (1.0 / nums.possibleInterAllCount) : 0);
        nums.baselineIntraChrProb = 1.0 / possibleIntraAllCount;
    } //TO DO: else
    
    end = clock();
    cpu_time_used += ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, "Fragments file read. Time took %f sec\n\n", cpu_time_used);
    
    fclose(fragsfile);

    fragReturn res;
    res.binStats = binStats;
    res.nums = nums;

    return res;
} 


float** calculateProbabilities(float** probTuple, var* mainDic, var** binStats, numbers nums)
{
    float* x = probTuple[0];
    float* y = probTuple[1];
  
    var sumCC, sumDistB4Scaling, possPairsInRange;
    float avgCC, avgDist;  
    for(var i = 0; i < nums.noOfBins; i++)
    {
        possPairsInRange = binStats[i][2];
        sumCC = binStats[i][3];
        sumDistB4Scaling = binStats[i][4];
        avgCC = ((float)sumCC/possPairsInRange)/nums.observedIntraInRangeSum;
        avgDist = (((float)sumDistB4Scaling)/possPairsInRange);
        //avgDist *= distScaling
        y[i] = avgCC;
        x[i] = avgDist;
    }

    return probTuple;
}



float** fit_Spline(var* mainDic, float* x_, float* y_, var passNo, numbers nums)
{
    fprintf(stdout, "Spline fit Pass %u starting...\n", passNo);
    clock_t start, end;
    double cpu_time_used = 0.0;
    start = clock();
    
    var k = 3;
    var n = nums.noOfBins - k;
    var i, j;

    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //fprintf(stdout, "Number of outliers is %lu\n", nums.outliersLength);
    fprintf(stdout, "Spline fit Pass %u completed. Time took %f\n\n", passNo, cpu_time_used);    
    float** probTuple = xmalloc(2 * sizeof(float*));
    probTuple[0] = x_;
    probTuple[1] = y_;

    return probTuple;
}
