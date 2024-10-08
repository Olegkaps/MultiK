//from StackOverflow
void free2d(size_t x, void *matrix) {
    if(!matrix) return;

    void **m = matrix;
    for(int i = 0; i < x; i++) {
        if(m[i]) xfree(m[i]);
    }
    xfree(m);

    return;
}

var** calloc2d(size_t x, size_t y, size_t elem_size) {
    var **p = xcalloc(x, sizeof(var*));

    for(int i = 0; i < x; i++) {
        p[i] = xcalloc(y, elem_size);
    }

    return (var**)p;
}

void memset2d_v(var** array, size_t x, int number, var end)
{
    for(var i = 0; i < x; i++)
    {
        memset(array[i], number, end * sizeof(var));
    }
}

void memset2d_f(float** array, size_t x, int number, var end)
{
    for(var i = 0; i < x; i++)
    {
        memset(array[i], number, end * sizeof(float));
    }
}
