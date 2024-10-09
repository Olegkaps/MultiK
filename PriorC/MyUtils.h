#define xscanf(argsCount, file, string, ...) {\
    var c = fscanf(file, string, __VA_ARGS__);\
    if (c != argsCount)\
    {\
        fprintf(stderr, "ERROR: Can't read file.\nExpected %u arguments, got %u.\nExpected string: \"%s\".\n", argsCount, c, (string));\
        exit(0);\
    }\
}


extern FILE* xopen(char* filePath, const char* flag);

extern void* xmalloc(size_t dataSize);

extern void* xcalloc(size_t dataSize, size_t itemSize);

extern void zfree(void* ptr);

extern void isNullStr(char* filePath);

#define xfree(ptr) {zfree(ptr); ptr = 0;}
#define zeroMem(a,n) 	 {memset(a,0,n*sizeof(*a));}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
