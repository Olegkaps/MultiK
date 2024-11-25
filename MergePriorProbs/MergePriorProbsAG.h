void showHelp(void)
{
                fprintf(stdout, "MergePriorProb helper:\n\n");
                fprintf(stdout, "-1\t--first\tREQUIRED: first spline file.\n\n");
                fprintf(stdout, "-2\t--second\tREQUIRED: second spline file.\n\n");
                fprintf(stdout, "-o\t--output\tREQUIRED: output spline file.\n\n");
                fprintf(stdout, "-h\t--help\tshow this message\n");
}

spline2merge ProcessArgs(spline2merge args) {
    isNullStr(args.first);
    isNullStr(args.second);
    isNullStr(args.output);


    FILE* firstFile = xopen(args.first, "r");
    FILE* secondFile = xopen(args.second, "r");
    var firstBinLength, firstBinCount;
    var secondBinLength, secondBinCount;


    xscanf(2, firstFile, "Bin Length: %d\nNum of Bins: %d\n", &firstBinLength, &firstBinCount);
    xscanf(2, secondFile, "Bin Length: %d\nNum of Bins: %d\n", &secondBinLength, &secondBinCount);
    fclose(firstFile);
    fclose(secondFile);
    
    if(firstBinLength != secondBinLength) {
        printf("ERROR: files with different bin size.\n");
        exit(0);
    }
    if(firstBinCount != secondBinCount) {
        printf("WARNING: files with different bin count, using minimum size.\nCheck are it appropriate files and try different bin size.\n");
    }
    args.binLength = firstBinLength;
    args.binCount = MIN(firstBinCount, secondBinCount);

    return args;
}


spline2merge parse(spline2merge args, int argc, char *argv[]) {
    args.first = NULL;
    args.second = NULL;

    extern char *optarg;
    extern int optind, opterr, optopt;
    int c;
    static struct option long_options[] =
    {
        {"first", required_argument, NULL, '1'},
        {"second", required_argument, NULL, '2'},
        {"output", required_argument, NULL, 'o'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, ":1:2:o:h", long_options, NULL)) != -1)
    {
        switch (c)
        {
            case '1':
                fprintf(stdout, "Reading first file from: %s\n", optarg);
                args.first = optarg; 
                break;
            case '2':
                fprintf(stdout, "Reading second file from: %s\n", optarg);
                args.second = optarg; 
                break;
            case 'o':
                fprintf(stdout, "Output file will be: %s\n", optarg);
                args.output = optarg; 
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

    return args;
}
