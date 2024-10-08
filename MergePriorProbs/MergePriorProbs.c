#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "../Utils/MyTypes.h"
#include "../Utils/MyUtils.h"
#include "MergePriorTypes.h"
#include "MergePriorProbsAG.h"
#include "MergePiorProbsReading.h"


int main(int argc, char *argv[]) {
    spline2merge args;
    args = parse(args, argc, argv);

    probTuple data;
    data.x = xmalloc(args.binCount * sizeof(float));
    data.y = xmalloc(args.binCount * sizeof(float));

    data = readSplines(data, args);

    writeSpline(data, args);

    return 0;
}
