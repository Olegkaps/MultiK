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
    var bin1, bin2, contactCount;
    var i = 0;
    while(!feof(input)) {
        xscanf(5, input,"%s\t%u\t%s\t%u\t%u\n", ch1, &bin1, ch2, &bin2, &contactCount);
        if(bin1 > bin2) {
            fprintf(up, "%s\t%u\t%s\t%u\t%u\n", ch1, bin1, ch2, bin2, contactCount);
        } else if(bin1 < bin2) {
            fprintf(down, "%s\t%u\t%s\t%u\t%u\n", ch1, bin1, ch2, bin2, contactCount);
        } else {
            if(i % 2) {
                fprintf(up, "%s\t%u\t%s\t%u\t%u\n", ch1, bin1, ch2, bin2, contactCount);
                i++;
            } else {
                fprintf(down, "%s\t%u\t%s\t%u\t%u\n", ch1, bin1, ch2, bin2, contactCount);
                i++;
            }
        }
    }

    xfree(input);
    xfree(down);
    xfree(up);

    printf("Done.\n");

    return 0;
}
