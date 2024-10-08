#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>



extern char *optarg;
extern int optind, opterr, optopt;
int c;
static struct option long_options[] =
{
    {"prior", required_argument, NULL, 'p'},
    {"multi", required_argument, NULL, 'm'},
    {"filename", required_argument, NULL, 'f'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}
};


void showHelp(void)
{
                fprintf(stdout, "Multi2Uni helper:\n\n");
                fprintf(stdout, "-p\t--prior\tREQUIRED: Prior built by PriorC.\n\n");
                fprintf(stdout, "-m\t--multi\tREQUIRED: Multi-mapping contact bin pair file.\n\n");
                fprintf(stdout, "-f\t--filename\tREQUIRED: Multi-mapping contact read pair posterior probability assignment output file name.\n\n");
                fprintf(stdout, "-h\t--help\tshow this message.\n");
}

DATA parseArgs(DATA data, int argc, char *argv[])
{
    data.pathes.priorPath = NULL; 
    data.pathes.multiPath = NULL;
    data.pathes.outputPath = NULL;
    data.pathes.outdir = ".";

    while ((c = getopt_long(argc, argv, ":p:m:f:h", long_options, NULL)) != -1)
    {
        switch (c)
        {
            case 'p':
                fprintf(stdout, "Reading Prior file from: %s\n", optarg);
                data.pathes.priorPath = optarg; 
                break;
            case 'm':
                fprintf(stdout, "Reading Multi file from: %s\n", optarg);
                data.pathes.multiPath = optarg;
                break;
            case 'f':
                fprintf(stdout, "Result file will be: %s\n", optarg);
                data.pathes.outputPath = optarg;
                break;
            case 'h':
                showHelp();
                exit(0);
            default:
                showHelp();
                exit(0);
        }
    }
    
    isNullStr(data.pathes.priorPath);
    isNullStr(data.pathes.multiPath);
    isNullStr(data.pathes.outputPath);


    return data;
}