float get_spline_prob(char* chr1, char* chr2, var dist, DATA data)
{
    if (strcmp(chr1, chr2) != 0)
    {
        return data.spline.interProb;
    }
    //binary search
    int start = 0;
    int end = data.spline.splineLength - 1;
    int position;
    while(start != end)
    {
        position  = (end + start) / 2;
        if (data.spline.distances[position] == dist)
        {
            return data.spline.probs[position];
        } else if (data.spline.distances[position] < dist)
        {
            start = position + 1;
        } else
        {
            end = position; 
        }
    }

    return data.spline.probs[start];
}

CSR init_matrix(CSR matrix, DATA data)
{
    matrix.Pi.pairID_ptr = xmalloc((data.nums.readsNum+1)*sizeof(var));
    matrix.Pi.pairID_ptr[0] = 0;
    matrix.Pi.SumID = xmalloc(data.nums.readsNum*sizeof(float));
    matrix.Pi.binPairs = xmalloc(data.nums.valuesNum*sizeof(var));
    matrix.Pi.pi_values = xmalloc(data.nums.valuesNum*sizeof(float));
    matrix.Z.binPairs_ptr = xmalloc((data.nums.binPairsNum+1)*sizeof(var));
    matrix.Z.binPairs_ptr[0] = 0;
    matrix.Z.sumBinPair = xmalloc(data.nums.binPairsNum*sizeof(float));
    matrix.Z.const_to_pi = xmalloc(data.nums.binPairsNum*sizeof(float));
    matrix.Z.pairID = xmalloc(data.nums.valuesNum*sizeof(var));
    //Z.z_values = malloc(valuesNum*sizeof(float));)

    return matrix;
}

void free_matrix(CSR matrix)
{
    xfree(matrix.Z.binPairs_ptr);
    xfree(matrix.Z.const_to_pi);
    xfree(matrix.Z.pairID);
    //free(Z.z_values);
    xfree(matrix.Pi.pairID_ptr);
    xfree(matrix.Pi.SumID);
    xfree(matrix.Pi.binPairs);
    xfree(matrix.Pi.pi_values);
}