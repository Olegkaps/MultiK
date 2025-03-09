#include "collecter.h"


// Abstract Writer class
class Writer {
    public:
        std::vector<std::string> outFiles;
        Writer(const std::vector<std::string> &outFiles) : outFiles(outFiles) {}
        virtual ~Writer() {}
        virtual void operator()(DataCollecter &collecter) = 0;
    };
    
    // RNA_Writer inherits from Writer.
    class RNA_Writer : public Writer {
    public:
        std::string multi;
    
        RNA_Writer(const std::vector<std::string> &outFiles) : Writer(outFiles) {
            assert(outFiles.size() == 1 && ("RNA writer: 1 file names needed, given " + std::to_string(outFiles.size())).c_str());
            multi = outFiles[0];
        }
    
        virtual void operator()(DataCollecter &collecter) override {
            // s5_files: get data for probability calculation.
            ProbCalcData uniData = static_cast<StandardRNACollecter&>(collecter).give_data_to_prob_calc();
            int MaxOfBins = uniData.maxBin;
            int BinLength = uniData.resolution;
    
            // s6_files: get data for EM.
            EMData emData = static_cast<StandardRNACollecter&>(collecter).give_data_to_em();
            int resolution = emData.resolution;
            int total = emData.total;
            int IDsNum = emData.IDsNum;
            int binPairsNum = emData.binPairsNum;
            int itemsNum = emData.itemsNum;
            std::map<PairKey, std::array<int, 3>> binPairs = emData.pairsData;
            std::map<std::string, std::vector<PairKey>> IDs = emData.multiIDs;
            
            std::ofstream multiStream(multi);
            if (!multiStream) {
                std::cerr << "Error opening file for writing: " << multi << std::endl;
                exit(1);
            }
            multiStream << "Bin Length: " << resolution << "\n";
            multiStream << "Total reads: " << total << "\n";
            multiStream << "Num of reads: " << IDsNum << "\n";
            multiStream << "Num of bin pairs: " << binPairsNum << "\n";
            multiStream << "Num of items: " << itemsNum << "\n";
            multiStream << "#reads" << "\n";

            std::map<PairKey, int> pairsIndeces;
            std::vector<PairKey> allPairs;
            int i = 0;
            for (const auto &item : binPairs) {
                PairKey pair = item.first;
                pairsIndeces[pair] = i;
                allPairs.push_back(pair);
                i++;
            }
            
            // For each read ID in IDs.
            for (const auto &entry : IDs) {
                const std::string &ID = entry.first;
                const std::vector<PairKey> &lst = entry.second;
                multiStream << lst.size() << " " << ID << "\n";
                for (size_t j = 0; j < lst.size(); j++) {
                    std::string chr1 = lst[j].chr1;
                    std::string chr2 = lst[j].chr2;
                    int bin1 = lst[j].bin1;
                    int bin2 = lst[j].bin2;
                    PairKey key = {chr1, bin1, chr2, bin2};
                    int ind = pairsIndeces[key];
                    multiStream << chr1 << " " << bin1 << " " << chr2 << " " << bin2 << " " 
                                << std::abs(bin1 - bin2) << " " << ind << "\n";
                }
            }

            multiStream << "#bins" << "\n";
           
            for (const auto &pair : allPairs) {
                // Determine if intra-chromosomal: 1 if chr1 equals chr2, else 0.
                int intra = (pair.chr1 == pair.chr2) ? 0 : -1;
                int diff = std::abs(pair.bin2 - pair.bin1);
                int sumUni = binPairs[pair][1] + binPairs[pair][2];
                
                multiStream << binPairs[pair][0] << " " << intra << " " << diff << " " << sumUni << "\n";
            }

            multiStream.close();
        }
    };
    


    
    class contacts_Writer : public Writer {
        public:
            std::string multi;
            std::string uni;
        
            contacts_Writer(const std::vector<std::string> &outFiles) : Writer(outFiles) {
                assert(outFiles.size() == 2 && ("RNA writer: 2 file names needed, given " + std::to_string(outFiles.size())).c_str());
                multi = outFiles[0];
                uni = outFiles[1];
            }
        
            virtual void operator()(DataCollecter &collecter) override {
                // s5_files: get data for probability calculation.
                ProbCalcData uniData = static_cast<StandardRNACollecter&>(collecter).give_data_to_prob_calc();
                int MaxOfBins = uniData.maxBin;
                int BinLength = uniData.resolution;
                std::map<PairKey, int> uniIntraContacts = uniData.uniIntraContacts;
        
                // s6_files: get data for EM.
                EMData emData = static_cast<StandardRNACollecter&>(collecter).give_data_to_em();
                int resolution = emData.resolution;
                int total = emData.total;
                int IDsNum = emData.IDsNum;
                int binPairsNum = emData.binPairsNum;
                int itemsNum = emData.itemsNum;
                std::map<PairKey, std::array<int, 3>> binPairs = emData.pairsData;
                std::map<std::string, std::vector<PairKey>> IDs = emData.multiIDs;
                
                std::ofstream multiStream(multi);
                if (!multiStream) {
                    std::cerr << "Error opening file for writing: " << multi << std::endl;
                    exit(1);
                }
                multiStream << "Bin Length: " << resolution << "\n";
                multiStream << "Total reads: " << total << "\n";
                multiStream << "Num of reads: " << IDsNum << "\n";
                multiStream << "Num of bin pairs: " << binPairsNum << "\n";
                multiStream << "Num of items: " << itemsNum << "\n";
                multiStream << "#reads" << "\n";
                
                std::map<PairKey, int> pairsIndeces;
                std::vector<PairKey> allPairs;
                int i = 0;
                for (const auto &item : binPairs) {
                    PairKey pair = item.first;
                    pairsIndeces[pair] = i;
                    allPairs.push_back(pair);
                    i++;
                }

                // For each read ID in IDs.
                for (const auto &entry : IDs) {
                    const std::string &ID = entry.first;
                    const std::vector<PairKey> &lst = entry.second;
                    multiStream << lst.size() << " " << ID << "\n";
                    for (size_t j = 0; j < lst.size(); j++) {
                        std::string chr1 = lst[j].chr1;
                        std::string chr2 = lst[j].chr2;
                        int bin1 = lst[j].bin1;
                        int bin2 = lst[j].bin2;
                        int dist = lst[j].dist;
                        PairKey key = {chr1, bin1, chr2, bin2, dist};
                        int ind = pairsIndeces[key];
                        multiStream << chr1 << " " << bin1 << " " << chr2 << " " << bin2 << " " 
                                    << dist << " " << ind << "\n";
                    }
                }

                multiStream << "#bins" << "\n";

                
                for (const auto &pair : allPairs) {
                    // Determine if intra-chromosomal: 1 if chr1 equals chr2, else 0.
                    int intra = (pair.chr1 == pair.chr2) ? 0 : -1;
                    int diff = pair.dist;
                    int sumUni = binPairs[pair][1] + binPairs[pair][2];
                    multiStream << pair.chr2 << " " << pair.bin2 << " " << binPairs[pair][0] << " " << intra << " " << diff << " " << sumUni << "\n";
                }
                
                multiStream.close();

                std::ofstream uniStream(uni);
                for (const auto &entry : uniIntraContacts) {
                    PairKey key = entry.first;
                    int uniCount = entry.second;
                    std::string chr1 = key.chr1;
                    int bin1 = key.bin1;
                    std::string chr2 = key.chr2;
                    int bin2 = key.bin2;
                    int dist = key.dist;
                    uniStream << chr1 << "\t" << bin1 << "\t" << chr2 << "\t" << bin2 << "\t" 
                                    << dist << "\t" << uniCount << "\n";
                }
                uniStream.close();
            }
        };
        