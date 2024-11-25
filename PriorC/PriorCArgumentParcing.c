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
#include "PriorCArgumentParcing.h"

void showHelp(void)
{
                fprintf(stdout, "PriorC helper:\n\n");
                fprintf(stdout, "-i\t--interactions=FileName\tREQUIRED: interactions between fragment pairs are read from INTERSFILE\n\n");
                fprintf(stdout, "-f\t--fragments=FileName\tREQUIRED: lengthes of chromosomes\n\n");
                fprintf(stdout, "-r\t--resolution=number\tREQUIRED: length of single bin\n\n");
                fprintf(stdout, "-o\t--outfile=FileName\tREQUIRED: output file name\n\n");
                fprintf(stdout, "-p\t--passes=Number\tOPTIONAL: number of spline passes to run. Default is 1\n\n");
                fprintf(stdout, "-b\t--noOfBins\tOPTIONAL: number of equal-occupancy (count) bins. Default is 100\n\n");
                fprintf(stdout, "-l\t--lib\tOPTIONAL: Name of the library that is analyzed to be used for name of file prefixes. DEFAULT is PriorC\n\n");
                fprintf(stdout, "-U\t--upperbound\tOPTIONAL: upper bound on the intra-chromosomal distance range (unit: base pairs). DEFAULT no limit. STRONGLY suggested to have a limit for large genomes, such as human/mouse. ex. '1000000, 5000000, etc.'\n\n");
                fprintf(stdout, "-L\t--lowerbound\tOPTIONAL: lower bound on the intra-chromosomal distance range (unit: base pairs). DEFAULT no limit. Suggested limit is 2x the resolution of the input files\n\n");
                fprintf(stdout, "-x\t--contactType\tOPTIONAL: use this flag to determine which chromosomal regions to study (intraOnly, interOnly, All) DEFAULT is intraOnly\n\n");
                fprintf(stdout, "-h\t--help\tshow this message\n");
}

priorCArgs ProcessArgs(priorCArgs args) {
    isNullStr(args.contactCountsFileName);
    isNullStr(args.fragsfileName);
    isNullStr(args.outputPath);
    if(args.binLength == 0) {
        showHelp();
        exit(0);
    }


    char chrName[10];
    var chrLength;
    FILE* fragsFile = xopen(args.fragsfileName, "r");    
    while(!feof(fragsFile)) {
        xscanf(2, fragsFile, "%s\t%u\n", chrName, &chrLength);
        args.binCount = MAX(args.binCount, chrLength / args.binLength);
        args.noOfFrags += chrLength / args.binLength;
    }
    //noOfFrags, args.binCount, <args.binLength>
    fclose(fragsFile);
    
    if(args.distUpThresNucl != NULL)
    {
        args.distUpThres = (var)(atol(args.distUpThresNucl) / args.binLength);
    } else
    {
        args.distUpThres = (var)INFINITY;
    }
    if(args.distLowThresNucl != NULL)
    {
        args.distLowThres = (var)(atol(args.distLowThresNucl) / args.binLength);
    } else
    {
        args.distLowThres = 0;
    }
    if(args.distLowThres >= args.distUpThres)
    {
        fprintf(stderr, "ERROR: Bad distance boundaries\n");
        exit(0);
    }

    return args;
}

priorCArgs init_args(priorCArgs args) {
    args.noOfPasses = 1;
    args.noOfBins = 100;
    args.binCount = 0;
    args.binLength = 0;
    args.noOfFrags = 0;
    args.libName = "PriorC";
    args.distLowThresNucl = NULL;
    args.distUpThresNucl = NULL;
    args.logfileName = NULL;
    args.contactCountsFileName = NULL;
    args.fragsfileName = NULL;
    args.outputPath;


    return args;
}


priorCArgs parse(priorCArgs args, int argc, char *argv[]) {
    args = init_args(args);
    
    fprintf(stdout, "\nGIVEN PriorC ARGUMENTS\n");
    fprintf(stdout, "=========================\n\n");
        //argparser
    extern char *optarg;
    extern int optind, opterr, optopt;
    int c;
    static struct option long_options[] =
    {
        {"interactions", required_argument, NULL, 'i'},
        {"fragments", required_argument, NULL, 'f'},
        {"resolution", required_argument, NULL, 'r'},
        {"outfile", required_argument, NULL, 'o'},
        {"passes", required_argument, NULL, 'p'},
        {"noOfBins", required_argument, NULL, 'b'},
        {"lib", required_argument, NULL, 'l'},
        {"upperbound", required_argument, NULL, 'U'},
        {"lowerbound", required_argument, NULL, 'L'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
    
    while ((c = getopt_long(argc, argv, ":i:f:r:o:p:b:l:U:L:h", long_options, NULL)) != -1)
    {
        switch (c)
        {
            case 'i':
                fprintf(stdout, "Reading interactions file from: %s\n", optarg);
                args.contactCountsFileName = optarg; 
                break;
            case 'f':
                fprintf(stdout, "Reading fragments file from: %s\n", optarg);
                args.fragsfileName = optarg;
                break;
            case 'r':
                fprintf(stdout, "Bin lenght: %s\n", optarg);
                args.binLength = (var)atoi(optarg);
                break;
            case 'o':
                fprintf(stdout, "Output file: %s\n", optarg);
                args.outputPath = optarg;
                break;
            case 'p':
                fprintf(stdout, "The number of spline passes is %s\n", optarg);
                args.noOfPasses = (var)atoi(optarg);
                break;
            case 'b':
                fprintf(stdout, "The number of bins is %s\n", optarg);
                args.noOfBins = (var)atoi(optarg);
                break;
            case 'l':
                fprintf(stdout, "The name of the library for outputted files will be %s\n", optarg);
                args.libName = optarg;
                break;
            case 'U':
                fprintf(stdout, "Upper Distance threshold is %s\n", optarg);
                args.distUpThresNucl = optarg;
                break;
            case 'L':
                fprintf(stdout, "Lower Distance threshold is %s\n", optarg);
                args.distLowThresNucl = optarg;
                break;
            case 'h':
                showHelp();
                exit(0);
            default:
                showHelp();
                exit(0);
        }
    }

    args = ProcessArgs(args);
    
    fprintf(stdout, "\nAll arguments processed.\n");
    fprintf(stdout, "=========================\n\n");

    return args;
}
