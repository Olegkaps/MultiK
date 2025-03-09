#include <bits/stdc++.h>
#include "writter.h"



// Abstract Reader class
class Reader {
    public:
        std::vector<std::string> inFiles;
        int resolution;
        Reader(const std::vector<std::string> &inFiles, int resolution) : inFiles(inFiles), resolution(resolution) {}
        virtual ~Reader() {}
        virtual void operator()(DataCollecter &collector) = 0;
    };
    
    

    class mHiC_Reader : public Reader {
    public:
        std::string multiName;
        std::string uniName;
        mHiC_Reader(const std::vector<std::string> &inFiles, int resolution) : Reader(inFiles, resolution) {
            assert(inFiles.size() == 2 && ("mHiC-reader: 2 file names needed, given " + std::to_string(inFiles.size())).c_str());
            multiName = inFiles[0];
            uniName = inFiles[1];
        }
        
        virtual void operator()(DataCollecter &collector) override {
            std::ifstream multiStream(multiName);
            if (!multiStream) {
                std::cerr << "Error opening file for reading: " << multiName << std::endl;
                exit(1);
            }
            std::ifstream uniStream(uniName);
            if (!uniStream) {
                std::cerr << "Error opening file for reading: " << uniName << std::endl;
                exit(1);
            }
            std::string line;
            // Process multi file lines.
            while (std::getline(multiStream, line)) {
                if(line.empty())
                    continue;
                std::istringstream iss(line);
                std::string ID, chr1, chr2;
                int mid1, mid2;
                if (!(iss >> ID >> chr1 >> mid1 >> chr2 >> mid2))
                    continue;
                int bin1 = (mid1 - (resolution / 2)) / resolution;
                int bin2 = (mid2 - (resolution / 2)) / resolution;
                collector(bin1, bin2, chr1, chr2, ID);
            }
            // Process uni file lines.
            while (std::getline(uniStream, line)) {
                if(line.empty())
                    continue;
                std::istringstream iss(line);
                std::string chr1, chr2;
                int mid1, mid2, count;
                if (!(iss >> chr1 >> mid1 >> chr2 >> mid2 >> count))
                    continue;
                int bin1 = (mid1 - (resolution / 2)) / resolution;
                int bin2 = (mid2 - (resolution / 2)) / resolution;
                
                collector(bin1, bin2, chr1, chr2, "", count);
            }
            multiStream.close();
            uniStream.close();
        }
    };


    class contacts_Reader : public Reader {
        public:
            std::string ContactsName;
            std::string annotation;
            std::string chr_ids;
            int rna_ann = -1;

            std::map<std::string, std::string> chr_id2chr_name;

            std::map<std::string, std::vector<RNA_Annotation>> RNAs;
            
            
            contacts_Reader(const std::vector<std::string> &inFiles, int resolution, std::string annotation, std::string chr_ids) : Reader(inFiles, resolution) {
                assert(inFiles.size() == 1 && ("mHiC-reader: 1 file names needed, given " + std::to_string(inFiles.size())).c_str());
                ContactsName = inFiles[0];
                if (annotation != "") {
                    rna_ann = 0;
                    std::ifstream annotationStream(annotation);
                    std::ifstream chr_ids_stream(chr_ids);
                    std::string line;

                    if (!annotationStream) {
                        std::cerr << "Error opening file for reading: " << annotation << std::endl;
                        exit(1);
                    } 
                    if (!chr_ids_stream) {
                        std::cerr << "Error opening file for reading: " << chr_ids << std::endl;
                        exit(1);
                    }

                    while (std::getline(chr_ids_stream, line)) {
                        std::istringstream iss(line);
                        std::string chr_id, chr_name;
                        if (!(iss >> chr_id >> chr_name))
                            continue;
                        chr_id2chr_name[chr_id] = chr_name;
                    }
                    chr_ids_stream.close();
                    std::string ending = "RNA";

                    while (std::getline(annotationStream, line)) {
                        std::istringstream iss(line);
                        if (line.rfind("#", 0) == 0) {
                            continue;
                        }
                        std::string chr_id, db, type, n1, s, n2, meta;
                        int start, end;
                        if (!(iss >> chr_id >> db >> type >> start >> end >> n1 >> s >> n2 >> meta))
                            continue;
                        std::string chr_name = chr_id2chr_name[chr_id];
                        if (RNAs.find(chr_name) == RNAs.end()) {
                            std::vector<RNA_Annotation> a;
                            RNAs[chr_name] = a;
                        }

                        if (type.compare(type.length() - ending.length(), ending.length(), ending)) {
                            continue;
                        }
                        std::string name = meta.substr(3, meta.find(';')-3); // ID=<name>;
                        if (RNAs[chr_name].size() > 0 && RNAs[chr_name].back().a_start == start) {
                            if (RNAs[chr_name].back().a_end < end) {
                                RNAs[chr_name].pop_back();
                                RNAs[chr_name].push_back({name, start, end});
                            }
                        } else {
                            RNAs[chr_name].push_back({name, start, end});
                        }
                    }
                    annotationStream.close();

                    std::ofstream rna_binsStream(annotation + ".multik");
                    if (!rna_binsStream) {
                        std::cerr << "Error opening file for reading: " << annotation << std::endl;
                        exit(1);
                    }
                    std::cout << "RNA bins in file: " << annotation + ".multik" << std::endl;
                    for(const auto &entry : RNAs) {
                        std::string chr = entry.first;
                        std::sort(RNAs[chr].begin(), RNAs[chr].end(), compareRNA_Annotation);
                    }

                    for(const auto &entry : RNAs) {
                        std::string chr = entry.first;
                        std::vector<RNA_Annotation> a = entry.second;

                        rna_binsStream << chr << std::endl;
                        for (int i = 0; i < a.size(); i++) {
                            rna_binsStream << i << "\t" << a[i].a_name << "\n";
                        }
                    }

                    rna_binsStream.close();
                }
            }

            int generate_bin(int mid, int resolution, int useRNAbin, std::string chr){
                if (rna_ann == 0 && useRNAbin == 1) {
                    std::vector<RNA_Annotation> a = RNAs[chr];
                    int start, end, position;
                    start = 0;
                    end = a.size()-1;
                    while(start < end)
                    {
                        position  = (end + start) / 2;
                        if (a[position].a_start == mid)
                        {
                            break;
                        } else if (a[position].a_start < mid)
                        {
                            start = position + 1;
                        } else
                        {
                            end = position; 
                        }
                    }
                    if (a[position].a_end > mid) {
                        return position;
                    } else {
                        return -1;
                    }
                }
                return mid / resolution;
            }
            
            virtual void operator()(DataCollecter &collector) override {
                std::string line;
                

                // contacts
                std::ifstream contactsStream(ContactsName);
                if (!contactsStream) {
                    std::cerr << "Error opening file for reading: " << ContactsName << std::endl;
                    exit(1);
                }
                // Process file lines.
                std::getline(contactsStream, line);
                // SRR_ID	pairtype	r1_chr	r1_start	r1_end	r1_cigar	r2_chr	r2_start
                // r2_end	r2_cigar	r1_suppl_alignments	r2_suppl_alignments	NH_r1	NH_r2	NM_r1	NM_r2
                while (std::getline(contactsStream, line)) {
                    if(line.empty())
                        continue;
                    std::istringstream iss(line);
                    std::string SSR_ID, pairtype, chr1, chr2;
                    std::string cigar1, cigar2, NH1, NH2, NM1, NM2;
                    char alignments1[2000], alignments2[2000];
                    int start1, start2, end1, end2;
                    
                    if (!(iss >> SSR_ID >> pairtype >> chr1 >> start1 >> end1 >> cigar1 >> chr2 >> start2 >> end2 >> cigar2))
                        continue;

                    if (pairtype == "MU") {
                        if (!(iss >> alignments1))
                            continue;
                        strcpy(alignments2, "");
                    } else if (pairtype == "UM") {
                        if (!(iss >> alignments2))
                            continue;
                        strcpy(alignments1, "");
                    } else if (pairtype == "MM") {
                        if (!(iss >> alignments1 >> alignments2))
                            continue;
                    }


                    int mid1, mid2, dist;
                    if (pairtype == "UU"){
                        // Unique contact
                        mid1 = (end1+start1)/2;
                        mid2 = (end2+start2)/2;
                        dist = std::abs(mid1-mid2) / resolution;
                        int bin1 = generate_bin(mid1, resolution, 1, chr1); // diff args to GFF bin
                        int bin2 = generate_bin(mid2, resolution, 0, chr2); 
                        collector(bin1, bin2, chr1, chr2, "", 1, dist); // 
                    } else if (pairtype == "MU" || pairtype == "UM" || pairtype == "MM") {
                        // multi contact
                        std::map<Bin, int> bins1 = {};
                        std::map<Bin, int> bins2 = {};
                        
                        mid1 = (end1+start1)/2;
                        mid2 = (end2+start2)/2;

                        int bin1 = generate_bin(mid1, resolution, 1, chr1);
                        int bin2 = generate_bin(mid2, resolution, 0, chr2);
                        
                        Bin bin_1, bin_2;
                        bin_1.chr = chr1;
                        bin_2.chr = chr2;
                        bin_1.bin = bin1;
                        bin_2.bin = bin2;
                        bin_1.mid = mid1/resolution;
                        bin_2.mid = mid2/resolution;
                        
                        bins1[bin_1] += 1;
                        bins2[bin_2] += 1;
                        // suplementary alignments
                    if (alignments1 != "") {
                            char * aligns = std::strtok(alignments1, "(),;");
                            while (aligns != NULL) {   
                                bin_1.chr = aligns; // chr
                                aligns = std::strtok(NULL, "(),;"); // mid
                                int mid = std::atoi(aligns);
                                bin_1.bin = generate_bin(mid, resolution, 1, bin_1.chr);
                                bin_1.mid = mid/resolution;
                                bins1[bin_1] += 1;
                                aligns = std::strtok(NULL, "(),;"); // sigar
                                aligns = std::strtok(NULL, "(),;"); // n
                                aligns = std::strtok(NULL, "(),;"); // chr again
                            }
                        }
                        if (alignments2 != "") {
                            char * aligns = std::strtok(alignments2, "(),;");
                            while (aligns != NULL) {   
                                bin_2.chr = aligns; // chr
                                aligns = std::strtok(NULL, "(),;"); // mid
                                int mid = std::atoi(aligns);
                                bin_2.bin = generate_bin(mid, resolution, 0, bin_2.chr);
                                bin_2.mid = mid/resolution;
                                bins2[bin_2] += 1;
                                aligns = std::strtok(NULL, "(),;"); // sigar
                                aligns = std::strtok(NULL, "(),;"); // n
                                aligns = std::strtok(NULL, "(),;"); // chr again
                            }
                        }
                        if (bins1.size()*bins2.size() > 10) { //maybe more than 10
                            continue;
                        }
                        for(auto item1: bins1) {
                            Bin _bin1 = item1.first;
                            for(auto item2: bins2) {
                                Bin _bin2 = item2.first;
                                collector(_bin1.bin, _bin2.bin, _bin1.chr, _bin2.chr, SSR_ID, 0, std::abs(_bin1.mid - _bin2.mid));
                            }
                        }


                    }
                }
                contactsStream.close();
            }
        };
        