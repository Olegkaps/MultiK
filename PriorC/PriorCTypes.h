#include "hashTable.h"

typedef struct {
    size_t outliersLength;
    var distLowThres;
    var distUpThres;
    var binCount;
    var binLength;
    var mainDicSize;
    var observedInterAllSum; 
    var observedIntraAllSum;   
    var observedIntraInRangeSum; 
    var noOfBins; 
    var noOfFrags;
    var possibleIntraInRangeCount;
    var possibleInterAllCount;
    float interChrProb;
    float baselineIntraChrProb;
} numbers;

typedef struct {
    var distLowThres;
    var distUpThres;
    var binCount;
    var binLength;
    const char* distLowThresNucl;
    const char* distUpThresNucl;
    char* logfileName;
    char* contactCountsFileName;
    char* fragsfileName;
    char* outputPath;
    var noOfPasses;
    var noOfBins;
    var noOfFrags;
    char* libName;
    char* dnaDensityFileName;
    DNA_chr_bin_hash_table* density_table;
} priorCArgs;


typedef struct {
    var* mainDic;
    numbers* nums;
} interactions;

typedef struct {
    var** binStats;
    numbers* nums;
} fragReturn;
