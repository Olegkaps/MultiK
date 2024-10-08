probTuple readSplines(probTuple data, spline2merge args) {
    FILE* firstFile = xopen(args.first, "r");
    FILE* secondFile = xopen(args.second, "r");
    float x1, x2, y1, y2;

    xscanf(2, firstFile, "Bin Length: %d\nNum of Bins: %d\n", &args.binLength, &args.binCount);
    xscanf(2, secondFile, "Bin Length: %d\nNum of Bins: %d\n", &args.binLength, &args.binCount);
    for(var i = 0; i < args.binCount; i++) {
        xscanf(2, firstFile, "%f\t%f\n", &x1, &y1);
        xscanf(2, secondFile, "%f\t%f\n", &x2, &y2);
        data.x[i] = (x1 + x2) / 2;
        data.y[i] = (y1 + y2) / 2;
    }

    fclose(firstFile);
    fclose(secondFile);
    
    return data;
}

void writeSpline(probTuple data, spline2merge args) {

    FILE* splineOut = xopen(args.output, "w");
    
    fprintf(splineOut, "Bin Length: %u\n", args.binLength);
    fprintf(splineOut, "Num of Bins: %u\n", args.binCount);
    for(var i = 0; i < args.binCount; i++)
    {
        if(data.y[i] > 0.0000000000000099999)
        {
            fprintf(splineOut, "%.0f\t%.9f\n", data.x[i], data.y[i]);
        }
    }

    fclose(splineOut);
    printf("Done.\n");
}