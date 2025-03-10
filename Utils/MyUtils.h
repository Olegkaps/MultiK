#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifndef my_utils_h
#define my_utils_h

extern void showHelp();

#define xscanf(argsCount, file, string, ...) {\
    var c = fscanf(file, string, __VA_ARGS__);\
    if (c != argsCount)\
    {\
        fprintf(stderr, "ERROR: Can't read file.\nExpected %u arguments, got %u.\nExpected string: \"%s\".\n", argsCount, c, (string));\
        exit(0);\
    }\
}


FILE* xopen(char* filePath, const char* flag){
    FILE* openedFile = fopen(filePath, flag);
    if (openedFile == NULL)
    {
        fprintf(stderr, "ERROR: can't open file \"%s\".\n", filePath);
        exit(0);
    }
    return openedFile;
}

void* xmalloc(size_t dataSize){
    void* ptr = malloc(dataSize);
    if (ptr == NULL)
    {
        fprintf(stderr, "ERROR: Can`t allocate memory.\n");
        exit(0);
    }
    return ptr;
}

void* xcalloc(size_t dataSize, size_t itemSize){
    void* ptr = calloc(dataSize, itemSize);
    if (ptr == NULL)
    {
        fprintf(stderr, "ERROR: Can`t allocate memory.\n");
        exit(0);
    }
    return ptr;
}

void zfree(void* ptr){
    if(ptr == NULL)
    {
        fprintf(stderr, "ERROR: Can`t free memory.\n");
        exit(0);
    }
    free(ptr);
}

void isNullStr(char* filePath)
{
    if(filePath == NULL)
    {
        fprintf(stderr, "ERROR: one of the files is NULL\n\n");
        showHelp();
        exit(0);
    }
}

#define xfree(ptr) {zfree(ptr); ptr = 0;}
#define zeroMem(a,n) 	 {memset(a,0,n*sizeof(*a));}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


#endif