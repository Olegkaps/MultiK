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
#include "EM.h"



int main(int argc, char *argv[])
{
    DATA* data = xmalloc(sizeof(DATA));
    data->dnas = xmalloc(sizeof(DNA_chr_bin_hash_table));
    
    data = parseArgs(data, argc, argv);
    data = read_spline_prior(data);
    if(strcmp(data->pathes.density, "")) {
        readDNADensity(data);
    }


    CSR* matrix = xmalloc(sizeof(CSR));
    

    FILE* multiFile = xopen(data->pathes.multiPath, "r");
    //var readsNum, binPairsNum, valuesNum, TotalReadsNum;
    xscanf(5, multiFile, "Bin Length: %u\nTotal reads: %u\nNum of reads: %u\nNum of bin pairs: %u\nNum of items: %u\n#reads\n", &data->spline.binLength2, &data->nums.TotalReadsNum, &data->nums.readsNum, &data->nums.binPairsNum, &data->nums.valuesNum);
    
    matrix = init_matrix(matrix, data);
    matrix = readReadsStrings(multiFile, matrix, data);
    matrix = readBinsStrings(multiFile, matrix, data);


    
    fclose(multiFile);
    
    if(data->spline.binLength1 != data->spline.binLength2)
    {
        printf("WARNING: you use files with different bin size!\n");
    }
    // do not use probs and distances anymore
    xfree(data->spline.probs);
    xfree(data->spline.distances);
    if(strcmp(data->pathes.density, "")) {
        free_table(data->dnas);
    }


    expectationMaximization(data, matrix);
    // calculate probabilities and write out the result
    FILE* multiReadsFile = xopen(data->pathes.multiPath, "r");
    FILE* outFile = xopen(data->pathes.outputPath, "w");

    xscanf(5, multiFile, "Bin Length: %u\nTotal reads: %u\nNum of reads: %u\nNum of bin pairs: %u\nNum of items: %u\n#reads\n", &data->spline.binLength2, &data->nums.TotalReadsNum, &data->nums.readsNum, &data->nums.binPairsNum, &data->nums.valuesNum);
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
