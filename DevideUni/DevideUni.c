#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "../Utils/MyUtils.h"
#include "../Utils/MyTypes.h"
#include "DevideUniTypes.h"
#include "DevideUniAG.h"


int main(int argc, char *argv[]) {
    devideUniArgs args;
    args = parse(args, argc, argv);

    FILE* input = xopen(args.input, "r");
    FILE* down = xopen(args.down, "w");
    FILE* up = xopen(args.up, "w");

    char ch1[100];
    char ch2[100];
    var bin1, bin2, dist, contactCount;
    while(!feof(input)) {
        xscanf(6, input,"%s\t%u\t%s\t%u\t%u\t%u\n", ch1, &bin1, ch2, &bin2, &dist, &contactCount);
        if (dist == 0) {
            fprintf(up, "%s\t%u\t%s\t%u\t%u\t%u\n", ch1, bin1, ch2, bin2, dist, contactCount);
            fprintf(down, "%s\t%u\t%s\t%u\t%u\t%u\n", ch1, bin1, ch2, bin2, dist, contactCount);
        } else if(bin1 > bin2) {
            fprintf(up, "%s\t%u\t%s\t%u\t%u\t%u\n", ch1, bin1, ch2, bin2, dist, contactCount);
        } else if(bin1 <= bin2) {
            fprintf(down, "%s\t%u\t%s\t%u\t%u\t%u\n", ch1, bin1, ch2, bin2, dist, contactCount);
        } 
    }

    fclose(input);
    fclose(down);
    fclose(up);

    printf("Done.\n");

    return 0;
}
