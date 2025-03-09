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


#include "structs.h"

// Abstract class DataCollecter
class DataCollecter {
    public:
        std::map<PairKey, int> uniPairsCountIntra; // default int = 0
        std::map<PairKey, int> uniPairsCountInter;
        std::map<std::string, std::vector<PairKey>> multiIDs;
        std::map<PairKey, int> multiPairs;
        int resolution;
        int maxBin;
        int total;
        
        DataCollecter(int resolution) : resolution(resolution), maxBin(-1), total(0) {}
        virtual ~DataCollecter() {}
    
        // Pure virtual operator() to process a record.
        virtual void operator()(int bin1, int bin2, const std::string &chr1, const std::string &chr2, 
                                const std::string &readID = "", int contactCount = 0, int dist = -1) = 0;
        
        // Pure virtual functions to give data to probability calculation and EM.
        virtual ProbCalcData give_data_to_prob_calc() = 0;
        virtual EMData give_data_to_em() = 0;
    };
    
    // StandardRNACollecter inheriting from DataCollecter
    class StandardRNACollecter : public DataCollecter {
    public:
        int IDsNum = 0;
        int itemsNum = 0;
        std::string lastID;
    
        StandardRNACollecter(int resolution) : DataCollecter(resolution), IDsNum(0), itemsNum(0), lastID("") {}
    
        // Overloaded operator() to process data records.
        virtual void operator()(int bin1, int bin2, const std::string &chr1, const std::string &chr2, 
                                const std::string &readID = "", int contactCount = 0, int dist = -1) override {
            if (!readID.empty()) { // multi file branch
                itemsNum += 1;
                // Append the bin pair data to multiIDs under the given readID.
                PairKey key = {chr1, bin1, chr2, bin2};
                multiIDs[readID].push_back(key);
                if (lastID != readID) {
                    //total += 1;
                    lastID = readID;
                    IDsNum += 1;
                }
                multiPairs[key] += 1;
            } else if (contactCount) { // uni file branch
                total += contactCount;
                
                if (chr1 == chr2) {
                    PairKey key = {chr1, bin1, chr2, bin2};
                    uniPairsCountIntra[key] += contactCount;
                } else {
                    PairKey key = {chr1, bin1, chr2, bin2};
                    uniPairsCountInter[key] += contactCount;
                }
            }
        }
    
        // Returns data for EM algorithm in the same order as Python yields.
        virtual EMData give_data_to_em() override {
            EMData data;
            data.resolution = resolution;
            data.total = total;
            data.IDsNum = IDsNum;
            data.itemsNum = itemsNum;
            // For each pair in multiPairs, collect counts from uni data maps.
            for (const auto &entry : multiPairs) {
                PairKey key = entry.first;
                int multiCount = entry.second;
                int uniCountIntra = 0;
                int uniCountInter = 0;
                if (uniPairsCountIntra.find(key) != uniPairsCountIntra.end()) {
                    uniCountIntra = uniPairsCountIntra[key];
                }
                if (uniPairsCountInter.find(key) != uniPairsCountInter.end()) {
                    uniCountInter = uniPairsCountInter[key];
                }
                data.pairsData[key] = {multiCount, uniCountIntra, uniCountInter};
            }
            data.binPairsNum = static_cast<int>(multiPairs.size());
            data.multiIDs = multiIDs;
            return data;
        }
    
        // Returns data for probability calculation in the same order as Python yields.
        virtual ProbCalcData give_data_to_prob_calc() override {
            ProbCalcData data;
            data.maxBin = maxBin;
            data.resolution = resolution;
            return data;
        }
    };
    


    class contactsRNACollecter : public DataCollecter {
        public:
            int IDsNum = 0;
            int itemsNum = 0;
            int UniIgnored = 0;
            int MultiIgnored = 0;
            std::string lastID;
            std::vector<UniRecord> uniContacts;
        
            contactsRNACollecter(int resolution) : DataCollecter(resolution), IDsNum(0), itemsNum(0), lastID("") {}
        
            // Overloaded operator() to process data records.
            virtual void operator()(int bin1, int bin2, const std::string &chr1, const std::string &chr2,
                                    const std::string &readID = "", int contactCount = 0, int dist = -1) override {
                if(bin1 == -1) {
                    if (!readID.empty()) {
                        MultiIgnored += 1;
                    } else if (contactCount) {
                        UniIgnored += contactCount;
                    }
                } else if (!readID.empty()) { // multi file branch
                    itemsNum += 1;
                    // Append the bin pair data to multiIDs under the given readID.
                    multiIDs[readID].push_back({chr1, bin1, chr2, bin2, dist});
                    if (lastID != readID) {
                        //total += 1;
                        lastID = readID;
                        IDsNum += 1;
                    }
                    PairKey key = {chr1, bin1, chr2, bin2, dist};
                    multiPairs[key] += 1;
                } else if (contactCount) { // uni file branch
                    total += contactCount;
                    
                    PairKey key = {chr1, bin1, chr2, bin2, dist};
                    uniPairsCountIntra[key] += contactCount;
              
                }
            }
        
            // Returns data for EM algorithm in the same order as Python yields.
            virtual EMData give_data_to_em() override {
                EMData data;
                data.resolution = resolution;
                data.total = total;
                data.IDsNum = IDsNum;
                data.itemsNum = itemsNum;
                std::cout << "Multi candidats deleted: " << MultiIgnored << "\n";
                std::cout << "Uni candidats deleted: " << UniIgnored << "\n";
                // For each pair in multiPairs, collect counts from uni data maps.
                for (const auto &entry : multiPairs) {
                    PairKey key = entry.first;
                    int multiCount = entry.second;
                    int uniCountIntra = 0;
                    if (uniPairsCountIntra.find(key) != uniPairsCountIntra.end()) {
                        uniCountIntra = uniPairsCountIntra[key];
                    }
                    data.pairsData[key] = {multiCount, uniCountIntra, 0};
                }
                data.binPairsNum = static_cast<int>(multiPairs.size());
                data.multiIDs = multiIDs;
                return data;
            }
        
            // Returns data for probability calculation in the same order as Python yields.
            virtual ProbCalcData give_data_to_prob_calc() override {
                ProbCalcData data;
                data.maxBin = maxBin;
                data.resolution = resolution;
                data.uniIntraContacts = uniPairsCountIntra;
                return data;
            }
        };
        