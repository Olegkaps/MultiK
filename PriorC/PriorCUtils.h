
extern void free2d(size_t x, void *matrix);

extern var** calloc2d(size_t x, size_t y, size_t elem_size);

extern void memset2d_v(var** array, size_t x, int number, var end);

extern void memset2d_f(float** array, size_t x, int number, var end);

#define xfree2d(X, M) {free2d(X, M); M=NULL;}
