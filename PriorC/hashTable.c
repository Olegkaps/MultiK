#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "MyUtils.h"
#include "MyTypes.h"

#ifndef hash_table_h
#define hash_table_h

static inline uint32_t murmur_32_scramble(uint32_t k) {
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    return k;
}

uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed)
{
	uint32_t h = seed;
    uint32_t k;
    /* Read in groups of 4. */
    for (size_t i = len >> 2; i; i--) {
        // Here is a source of differing results across endiannesses.
        // A swap here has no effects on hash properties though.
        memcpy(&k, key, sizeof(uint32_t));
        key += sizeof(uint32_t);
        h ^= murmur_32_scramble(k);
        h = (h << 13) | (h >> 19);
        h = h * 5 + 0xe6546b64;
    }
    /* Read the rest. */
    k = 0;
    for (size_t i = len & 3; i; i--) {
        k <<= 8;
        k |= key[i - 1];
    }
    // A swap is *not* necessary here because the preceding loop already
    // places the low bytes in the low places according to whatever endianness
    // we use. Swaps only apply when the memory is copied in a chunk.
    h ^= murmur_32_scramble(k);
    /* Finalize. */
	h ^= len;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

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


void insert_chr_to_DNA_htable(DNA_chr_bin_hash_table* DNAs, char* chr, var numOfBins) {
    var hash = (var)murmur3_32(chr, strlen(chr), 0x9747b28c);
    var pos = hash % DNAs->len;
    int i = 0;
    while(DNAs->hashes[pos] != 0) {
        pos = (pos + 2*i + 1) % DNAs->len;
        i++;
    }
    DNAs->hashes[pos] = hash;
    DNAs->c[pos].chr = xcalloc(30, sizeof(char));
    strcpy(DNAs->c[pos].chr, chr);
    DNAs->c[pos].chr_bin->Count = xcalloc(numOfBins, sizeof(var));
    DNAs->c[pos].chr_bin->Density = xcalloc(numOfBins, sizeof(float));
    DNAs->c[pos].chr_bin->len = numOfBins;
}

DNA_bins* find_bins_of_chr(DNA_chr_bin_hash_table* DNAs, char* chr) {
    var hash = (var)murmur3_32(chr, strlen(chr), 0x9747b28c);
    var pos = hash % DNAs->len;
    int i = 0;
    while(DNAs->hashes[pos] != hash && strcmp(DNAs->c[pos].chr, chr)) {
        if(DNAs->hashes[pos] == 0) {
            fprintf(stderr, "ERROR: Can`t find density to chromosome: %s\n", chr);
            exit(0);
        }
        pos = (pos + 2*i + 1) % DNAs->len;
        i++;
    }
    return DNAs->c[pos].chr_bin;
}

DNA_chr_bin_hash_table* init_dna(DNA_chr_bin_hash_table* DNAs, var N) {
    DNAs->hashes = xcalloc(2*N, sizeof(var));
    DNAs->c = xmalloc(2*N*sizeof(chrs));
    for(var i = 0; i < 2*N; i++) {
        DNAs->c[i].chr_bin = (DNA_bins*)xmalloc(sizeof(DNA_bins));
        DNAs->c[i].chr = "";
    }
    DNAs->len = 2*N;

    return DNAs;
}

DNA_chr_bin_hash_table* init_hash_table(char* chromosomesFileName, var resolution) {
    FILE *chroms = xopen(chromosomesFileName, "r");
    var N = 0;
    char *chr = xmalloc(30*sizeof(char));
    var chr_len;
    while(!feof(chroms)) {
        xscanf(2, chroms, "%s\t%u\n", chr, &chr_len);
        N++;
    }
    fclose(chroms);
    DNA_chr_bin_hash_table* DNAs = xmalloc(sizeof(DNA_chr_bin_hash_table));
    DNAs = init_dna(DNAs, N);

    chroms = xopen(chromosomesFileName, "r");
    while(!feof(chroms)) {
        xscanf(2, chroms, "%s\t%u\n", chr, &chr_len);
        insert_chr_to_DNA_htable(DNAs, chr, chr_len/resolution);
    }
    fclose(chroms);
    xfree(chr);
    
    return DNAs;
}

DNA_bins* compute_density(DNA_bins* dna, double* dists, double* probs) {
    var lenght = dna->len;
    var diff, pos;
    for(int k = 0; k < lenght; k++) {
        var count = dna->Count[k];
        if(count == 0) {
            continue;
        }
        diff = 0;
        pos = 0;
        dna->Density[k] += (float)count;
        for (int l = k + 1; l < lenght; l++) {
            if (diff > (var)dists[pos]) {
                pos++;
            }
            if(probs[pos] < 0.01) {
                break;
            }
            dna->Density[l] += (float)count*probs[pos];
            diff++;
        }
        diff = 0;
        pos = 0;
        for (int l = k-1; l >= 0; l--) {
            if (diff > (var)dists[pos]) {
                pos++;
            }
            if(probs[pos] < 0.01) {
                break;
            }
            dna->Density[l] += (float)count*probs[pos];
            diff++;
        }
    }
    return dna;
}

void free_table(DNA_chr_bin_hash_table* dnas) {
    for(var i = 0; i < dnas->len; i++) {
        if(dnas->hashes[i] != 0) {
            xfree(dnas->c[i].chr_bin->Count);
            xfree(dnas->c[i].chr_bin->Density);
            xfree(dnas->c[i].chr_bin);
        }
    }
    xfree(dnas->c);
    xfree(dnas->hashes);
}

#endif