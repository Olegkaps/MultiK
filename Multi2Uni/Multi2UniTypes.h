typedef struct {
    var* pairID_ptr;
    float* SumID;
    var* binPairs;
    float* pi_values;
} CSR_Pi;

typedef struct {
    var* binPairs_ptr;
    float* sumBinPair;
    float* const_to_pi;
    var* pairID;
    //float* z_values;
} CSR_Z;

typedef struct {
    CSR_Pi Pi;
    CSR_Z Z;
} CSR;

typedef struct {
    var* distances;
    float* probs;
    float interProb;
    var splineLength;
    var binLength1;
    var binLength2;
} splineData;

typedef struct {
    char *priorPath;
    char *multiPath;
    char *outputPath;
    char *outdir;
} pathesData;

typedef struct {
    var binPairsNum;
    var readsNum;
    var valuesNum;
    var TotalReadsNum;
    var readsPosition;
} numbersData;


typedef struct {
    splineData spline;
    pathesData pathes;
    numbersData nums;
} DATA;

