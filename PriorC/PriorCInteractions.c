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
#include "PriorCInteractions.h"


var** makeBinsFromInteractions(var** binStats, var* mainDic, var* outliersdist, size_t  outliersdistSize, numbers* nums)
{
    var desiredPerBin = (var)((float)nums->observedIntraInRangeSum / (float)nums->noOfBins);

    var interactionTotalForBinTermination = 0;
    var n = 0; // BIN COUNTER SO FAR
    var totalInteractionCountSoFar = 0;
    var currentBinContactCount = 0;

   
    binStats[0][0] = nums->distLowThres;

    for(var i = 0; i < nums->mainDicSize; i++)
    {
        currentBinContactCount += mainDic[i];
        totalInteractionCountSoFar += mainDic[i];
        interactionTotalForBinTermination += mainDic[i];
        if(interactionTotalForBinTermination >= desiredPerBin)
        {
            interactionTotalForBinTermination = 0;
            binStats[n][1] = nums->distLowThres + i;
            binStats[n][3] = currentBinContactCount;
            n += 1;
            currentBinContactCount = 0;
            if(n < nums->noOfBins)
            {
                desiredPerBin = (float)(nums->observedIntraInRangeSum-totalInteractionCountSoFar) / (float)(nums->noOfBins - n);
                if(desiredPerBin == 0)
                {
                    break;
                }
                binStats[n][0] = binStats[n - 1][1] + 1;

            }
            
            
            
        }
    }

    // TO DO: OUTLIERS
    
    return binStats;
}
