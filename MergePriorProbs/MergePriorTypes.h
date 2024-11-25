typedef struct {
    char* first;
    char* second;
    char* output;
    var binLength;
    var binCount;
} spline2merge;

typedef struct {
    float* x;
    float* y;
    spline2merge args;
} probTuple;
