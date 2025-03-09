#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <tuple>


struct PairKey {
    std::string chr1;
    int bin1;
    std::string chr2;
    int bin2;
    int dist;
    bool operator<(const PairKey &other) const {
        if (chr1 < other.chr1)
            return true;
        if (chr1 > other.chr1)
            return false;
        if (bin1 < other.bin1)
            return true;
        if (bin1 > other.bin1)
            return false;
        if (chr2 < other.chr2)
            return true;
        if (chr2 > other.chr2)
            return false;
        return bin2 < other.bin2;
    }
};


struct Bin {
    std::string chr;
    int bin;
    int mid;
    bool operator<(const Bin &other) const {
        if (chr < other.chr)
            return true;
        if (chr > other.chr)
            return false;
        if (bin < other.bin)
            return true;
        if (bin > other.bin)
            return false;
        return mid < other.mid;
    }
};



struct UniRecord {
    std::string chr1;
    int bin1;
    std::string chr2;
    int bin2;
    int contactCount;
};


struct EMData {
    int resolution;
    int total;
    int IDsNum;
    int binPairsNum;
    int itemsNum;
    std::map<PairKey, std::array<int, 3>> pairsData; // [multiPairs, uniCountIntra, uniCountInter]
    std::map<std::string, std::vector<PairKey>> multiIDs;
};


struct ProbCalcData {
    int maxBin;
    int resolution;
    std::map<PairKey, int> uniIntraContacts;
};


struct RNA_Annotation
{
    std::string a_name;
    int a_start;
    int a_end;
};

bool compareRNA_Annotation(const RNA_Annotation &a, const RNA_Annotation &b)
{
    return a.a_start < b.a_start;
}
