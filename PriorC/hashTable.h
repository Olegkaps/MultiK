#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "MyUtils.h"
#include "MyTypes.h"

#ifndef hash_table_h
#define hash_table_h

extern uint32_t murmur_32_scramble(uint32_t k);

extern uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed);

typedef struct {
    var *Count;
    float *Density;
    var len;
} DNA_bins;

typedef struct {
    char *chr;
    DNA_bins *chr_bin;
} chrs;

typedef struct {
    var *hashes;
    var len;
    chrs *c;
} DNA_chr_bin_hash_table;



extern void insert_chr_to_DNA_htable(DNA_chr_bin_hash_table* DNAs, char* chr, var numOfBins);

extern DNA_bins* find_bins_of_chr(DNA_chr_bin_hash_table* DNAs, char* chr);

extern DNA_chr_bin_hash_table* init_hash_table(char* chromosomesFileName, var resolution);

extern DNA_bins* compute_density(DNA_bins* dna, double* dists, double* probs);

extern void free_table(DNA_chr_bin_hash_table* dnas);

#endif