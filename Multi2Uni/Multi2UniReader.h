#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>



DATA read_spline_prior(DATA data)
{
    FILE* splineFile = xopen(data.pathes.priorPath, "r");
    xscanf(2, splineFile, "Bin Length: %u\nNum of Bins: %u\n", &data.spline.binLength1, &data.spline.splineLength);

    data.spline.interProb = INFINITY;

    data.spline.distances = xmalloc(data.spline.splineLength * sizeof(var));
    data.spline.probs = xmalloc(data.spline.splineLength * sizeof(data.spline.probs));
    float prob; 
    var i = 0;
    while(!feof(splineFile))
    {
        xscanf(2, splineFile, "%f\t%f\n", &data.spline.distances[i], &prob);
        data.spline.probs[i] = prob;
        data.spline.interProb = MIN(data.spline.interProb, prob);
        i++;
        if(i >= data.spline.splineLength) {
           break;
        }
    }

    data.spline.splineLength = i;
    fclose(splineFile);
    printf("inter prob %f\n", data.spline.interProb);

    return data;
}


CSR readBinsStrings(FILE* multiFile, CSR matrix, DATA data)
{
    for(var i = 0; i < data.nums.binPairsNum; i++)
    {
        var currBinPairLen, bin1, bin2, intraContactsNum, interContactsNum;

        xscanf(5, multiFile, "%u %u %u %u %u\n", &currBinPairLen, &bin1, &bin2, &intraContactsNum, &interContactsNum);
        //bin2 >= bin1
        
        matrix.Z.binPairs_ptr[i+1] = matrix.Z.binPairs_ptr[i] + currBinPairLen;
        matrix.Z.const_to_pi[i] = intraContactsNum * get_spline_prob("0", "0", abs(bin2 - bin1), data) + interContactsNum * data.spline.interProb;
    }

    return matrix;
}


CSR readReadsStrings(FILE* multiFile, CSR matrix, DATA data)
{
    var currReadLen, binPairIndex;
    char ID[100], chr1[100], chr2[100];
    var* binPairsRead = xcalloc(data.nums.binPairsNum, sizeof(var));
    var bin1, bin2;
    var ind;
    var idx;

    for(var i = 0; i < data.nums.readsNum; i++)
    {
        xscanf(1, multiFile, "%u\n", &currReadLen);
        matrix.Pi.pairID_ptr[i+1] = matrix.Pi.pairID_ptr[i] + currReadLen;
        matrix.Pi.SumID[i] = 0.0;
        for(var j = 0; j < currReadLen; j++)
        {
            xscanf(6, multiFile, "%s %s %u %s %u %u\n", ID, chr1, &bin1, chr2, &bin2, &binPairIndex);
            
            //Pi setter
            ind = matrix.Pi.pairID_ptr[i] + j;
            matrix.Pi.binPairs[ind] = binPairIndex;
            matrix.Pi.pi_values[ind] = get_spline_prob(chr1, chr2, abs(bin2-bin1), data);
            matrix.Pi.SumID[i] += matrix.Pi.pi_values[ind];
            //Z setter
            idx = binPairsRead[binPairIndex] + matrix.Z.binPairs_ptr[binPairIndex];
            matrix.Z.pairID[idx] = i;
            binPairsRead[binPairIndex] += 1;
        }
        matrix.Pi.SumID[i] = 1/matrix.Pi.SumID[i];
    }
    xfree(binPairsRead);

    return matrix;
}


void readAndWriteReads(FILE* multiReadsFile, FILE* outFile, CSR matrix, DATA data)
{
    var currReadLen, binPairIndex;
    char ID[100], chr1[100], chr2[100];
    var* binPairsRead = xcalloc(data.nums.binPairsNum, sizeof(var));
    var bin1, bin2;
    var ind;
    var idx;
    float resultProb;
    var readsPosition;

    xscanf(1, multiReadsFile, "#reads %u\n", &readsPosition);
    for(var i = 0; i < data.nums.readsNum; i++)
    {
        xscanf(1, multiReadsFile, "%u\n", &currReadLen);
        for(var j = matrix.Pi.pairID_ptr[i]; j < matrix.Pi.pairID_ptr[i+1]; j++)
        {
            resultProb = matrix.Pi.pi_values[j] * matrix.Pi.SumID[i];
            xscanf(6, multiReadsFile, "%s %s %u %s %u %u\n", ID, chr1, &bin1, chr2, &bin2, &binPairIndex);
            fprintf(outFile, "%s %s %u %s %u %.9f\n", ID, chr1, bin1, chr2, bin2, resultProb);
        }
    }
}
