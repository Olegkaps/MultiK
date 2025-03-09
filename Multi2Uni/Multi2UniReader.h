#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>



DATA* read_spline_prior(DATA* data) {
    FILE* splineFile = xopen(data->pathes.priorPath, "r");
    xscanf(2, splineFile, "Bin Length: %u\nNum of Bins: %u\n", &data->spline.binLength1, &data->spline.splineLength);

    data->spline.interProb = INFINITY;

    data->spline.distances = xmalloc(data->spline.splineLength * sizeof(var));
    data->spline.probs = xmalloc(data->spline.splineLength * sizeof(data->spline.probs));
    float prob, dist; 
    var i = 0;
    while(!feof(splineFile))
    {
        xscanf(2, splineFile, "%f\t%f\n", &dist, &prob);
        data->spline.distances[i] = (var)dist;
        data->spline.probs[i] = prob;
        data->spline.interProb = MIN(data->spline.interProb, prob);
        i++;
        if(i >= data->spline.splineLength) {
           break;
        }
    }

    data->spline.splineLength = i;
    fclose(splineFile);
    printf("inter prob %f\n", data->spline.interProb);

    return data;
}

void readDNADensity(DATA* data) {
    FILE* density = xopen(data->pathes.density, "r");
    var curr_len;
    var len;
    char* chr = xmalloc(30*sizeof(char));
    xscanf(1, density, "len is %u\n", &len);
    data->dnas = init_dna(data->dnas, len);
        
    data->nums.averageDensity = 0.0;
    data->nums.dnaBinsCount = 0;

    float d;
    for(var i = 0; i < len; i++) {
        xscanf(2, density, "%s %u\n", chr, &curr_len);
        insert_chr_to_DNA_htable(data->dnas, chr, curr_len);
        DNA_bins* dna_bins = find_bins_of_chr(data->dnas, chr);
        for (var j = 0; j < curr_len; j++) {
            xscanf(1, density, "%f\n", &d);
            dna_bins->Density[j] = d;
            data->nums.averageDensity += d;
            data->nums.dnaBinsCount++;
        }
    }
    fclose(density);
    data->nums.averageDensity /= data->nums.dnaBinsCount;

    xfree(chr);
}


CSR* readBinsStrings(FILE* multiFile, CSR* matrix, DATA* data)
{
    fscanf(multiFile, "#bins\n");
    var* binPairsRead = xcalloc(data->nums.binPairsNum, sizeof(var));

    var dna_bin, currBinPairLen, dist, uniContactsNum;
    int isIntraContact;
    char* chr = xmalloc(30*sizeof(char));

    for(var i = 0; i < data->nums.binPairsNum; i++)
    {
        xscanf(6, multiFile, "%s %u %u %d %u %u\n", chr, &dna_bin, &currBinPairLen, &isIntraContact, &dist, &uniContactsNum);   
        
        matrix->Z.binPairs_ptr[i+1] = matrix->Z.binPairs_ptr[i] + currBinPairLen;
        matrix->Z.const_to_pi[i] = uniContactsNum + data->nums.TotalReadsNum * get_spline_prob(isIntraContact, dist, data, dna_bin, chr);
    }
    var binPairIndex, idx;
    for (var i = 0; i < data->nums.readsNum; i++) {
        for(var j = matrix->Pi.pairID_ptr[i]; j < matrix->Pi.pairID_ptr[i+1]; j++) {
            binPairIndex = matrix->Pi.binPairs[j];
            //Z setter
            idx = binPairsRead[binPairIndex] + matrix->Z.binPairs_ptr[binPairIndex];
            matrix->Z.pairID[idx] = i;
            binPairsRead[binPairIndex] += 1;
        }
    }


    xfree(binPairsRead);
    xfree(chr);

    return matrix;
}


CSR* readReadsStrings(FILE* multiFile, CSR* matrix, DATA* data)
{
    var currReadLen, binPairIndex;
    char* ID = xmalloc(30*sizeof(char));
    char* chr1 = xmalloc(30*sizeof(char)); 
    char* chr2 = xmalloc(30*sizeof(char));
    var bin1, bin2;
    var ind;
    var idx;
    var dist;


    for(var i = 0; i < data->nums.readsNum; i++)
    {
        xscanf(2, multiFile, "%u %s\n", &currReadLen, ID);
        matrix->Pi.pairID_ptr[i+1] = matrix->Pi.pairID_ptr[i] + currReadLen;
        matrix->Pi.SumID[i] = 0.0;
        for(var j = 0; j < currReadLen; j++)
        {
            xscanf(6, multiFile, "%s %u %s %u %u %u\n", chr1, &bin1, chr2, &bin2, &dist, &binPairIndex);
            
            //Pi setter
            ind = matrix->Pi.pairID_ptr[i] + j;
            matrix->Pi.binPairs[ind] = binPairIndex;
            matrix->Pi.pi_values[ind] = get_spline_prob(strcmp(chr1, chr2), dist, data, bin2, chr2);
            matrix->Pi.SumID[i] += matrix->Pi.pi_values[ind];
        }
        matrix->Pi.SumID[i] = 1/matrix->Pi.SumID[i];
    }

    xfree(ID);
    xfree(chr1);
    xfree(chr2);

    return matrix;
}


void readAndWriteReads(FILE* multiReadsFile, FILE* outFile, CSR* matrix, DATA* data)
{
    var currReadLen, binPairIndex;
    char* ID = xmalloc(30*sizeof(char));
    char* chr1 = xmalloc(30*sizeof(char));
    char* chr2 = xmalloc(30*sizeof(char));
    var* binPairsRead = xcalloc(data->nums.binPairsNum, sizeof(var));
    var bin1, bin2;
    var ind;
    var idx;
    var dist;
    float resultProb;
    var readsPosition;

    for(var i = 0; i < data->nums.readsNum; i++)
    {
        xscanf(2, multiReadsFile, "%u %s\n", &currReadLen, ID);
        for(var j = matrix->Pi.pairID_ptr[i]; j < matrix->Pi.pairID_ptr[i+1]; j++)
        {
            resultProb = matrix->Pi.pi_values[j] * matrix->Pi.SumID[i];
            xscanf(6, multiReadsFile, "%s %u %s %u %u %u\n", chr1, &bin1, chr2, &bin2, &dist, &binPairIndex);
            fprintf(outFile, "%s\t%s\t%u\t%s\t%u\t%.9f\n", ID, chr1, bin1, chr2, bin2, resultProb);
        }
    }
    xfree(ID);
    xfree(chr1);
    xfree(chr2);
}
