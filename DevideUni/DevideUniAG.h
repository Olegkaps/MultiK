void showHelp(void)
{
                fprintf(stdout, "MergePriorProb helper:\n\n");
                fprintf(stdout, "-1\t--down\tREQUIRED: output down Uni file.\n\n");
                fprintf(stdout, "-2\t--up\tREQUIRED: output up Uni file.\n\n");
                fprintf(stdout, "-i\t--input\tREQUIRED: input Uni file.\n\n");
                fprintf(stdout, "-h\t--help\tshow this message\n");
}

devideUniArgs ProcessArgs(devideUniArgs args) {
    isNullStr(args.down);
    isNullStr(args.up);
    isNullStr(args.input);

    return args;
}


devideUniArgs parse(devideUniArgs args, int argc, char *argv[]) {
    args.up = NULL;
    args.down = NULL;

    extern char *optarg;
    extern int optind, opterr, optopt;
    int c;
    static struct option long_options[] =
    {
        {"down", required_argument, NULL, '1'},
        {"second", required_argument, NULL, '2'},
        {"input", required_argument, NULL, 'i'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, ":1:2:i:h", long_options, NULL)) != -1)
    {
        switch (c)
        {
            case '1':
                fprintf(stdout, "Writting down file to: %s\n", optarg);
                args.down = optarg; 
                break;
            case '2':
                fprintf(stdout, "Writting second file to: %s\n", optarg);
                args.up = optarg; 
                break;
            case 'i':
                fprintf(stdout, "Input file is: %s\n", optarg);
                args.input = optarg; 
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