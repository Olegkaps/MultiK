from collections import defaultdict
from abc import abstractmethod
from typing import Type
import argparse


class DataCollecter:
    def __init__(self, resolution: int):
        self.uniPairsCountIntra = defaultdict(int)
        self.uniPairsCountInter = defaultdict(int)
        self.multiIDs = defaultdict(list)
        self.multiPairs = defaultdict(int)
        self.resolution = resolution
        self.maxBin = -1
        self.total = 0

        @abstractmethod
        def __call__(self, bin1: str, bin2: str, chr1: str, chr2: str, readID: str = None, contactCount: int = None):
            pass

        @abstractmethod
        def give_data_to_prob_calc(self):
            pass

        @abstractmethod
        def give_data_to_em(self):
            pass


class StandardRNACollecter(DataCollecter):
    def __init__(self, resolution: int):
        super().__init__(resolution)
        self.uniUp = []
        self.uniDown = []
        self.IDsNum = 0
        self.itemsNum = 0
        self.lastID = None

    def __call__(self, bin1: int, bin2: int, chr1: str, chr2: str, readID: str = None, contactCount: int = None):
        if readID:  # multi file
            self.itemsNum += 1
            self.multiIDs[readID].append([chr1, bin1, chr2, bin2])
            if self.lastID != readID:
                self.total += 1
                self.lastID = readID
                self.IDsNum += 1
            bin1, bin2 = min(bin1, bin2), max(bin1, bin2)
            self.multiPairs[(bin1, bin2)] += 1
        elif contactCount:  # uni file
            self.total += 1
            self.maxBin = max(self.maxBin, bin1, bin2)
            if bin1 - bin2 > 0:
                self.uniUp.append([chr1, bin1, chr2, bin2, contactCount])
            else:
                self.uniDown.append([chr1, bin1, chr2, bin2, contactCount])
            bin1, bin2 = min(bin1, bin2), max(bin1, bin2)
            if chr1 == chr2:
                self.uniPairsCountIntra[(bin1, bin2)] += contactCount
            else:
                self.uniPairsCountInter[(bin1, bin2)] += contactCount

    def give_data_to_em(self):
        self.pairsData = {}
        for pair in self.multiPairs:
            uniCountIntra = 0 if pair not in self.uniPairsCountIntra else self.uniPairsCountIntra[pair]
            uniCountInter = 0 if pair not in self.uniPairsCountInter else self.uniPairsCountInter[pair]
            self.pairsData[pair] = [self.multiPairs[pair], uniCountIntra, uniCountInter]
        self.binPairsNum = len(self.multiPairs)
        yield self.resolution
        yield self.total
        yield self.IDsNum
        yield self.binPairsNum
        yield self.itemsNum
        yield self.pairsData
        yield self.multiIDs

    def give_data_to_prob_calc(self):
        yield self.maxBin
        yield self.resolution
        yield self.uniUp
        yield self.uniDown


class Writer:
    def __init__(self, outFiles: list):
        self.outFiles = outFiles

    @abstractmethod
    def __call__(self, collecter):
        pass


class RNA_Writer(Writer):
    def __init__(self, outFiles: list):
        super().__init__(outFiles)
        assert len(self.outFiles) == 3, f"RNA writer: 3 file names needed, given {len(self.outFiles)}"
        self.upUni = self.outFiles[0]
        self.downUni = self.outFiles[1]
        self.multi = self.outFiles[2]

    def __call__(self, collecter):
        # s5_files
        uni_gen = collecter.give_data_to_prob_calc()
        MaxOfBins = next(uni_gen)
        BinLength = next(uni_gen)
        uniUpList = next(uni_gen)
        uniDownList = next(uni_gen)

        #upUni = open(self.upUni, "w")
        #print(f"Number of bins: {MaxOfBins}", file=upUni)
        #print(f"Length of a single bin: {BinLength}", file=upUni)
        #for binPair in uniUpList:
            #print(*binPair, sep="\t", file=upUni)
        #upUni.close()
        #downUni = open(self.downUni, "w")
        #print(f"Number of bins: {MaxOfBins}", file=downUni)
        #print(f"Length of a single bin: {BinLength}", file=downUni)
        #for binPair in uniDownList:
            #print(*binPair, sep="\t", file=downUni)
        #downUni.close()
        # frag files needed

        # s6_files
        mul_gen = collecter.give_data_to_em()
        resolution = next(mul_gen)
        total = next(mul_gen)
        IDsNum = next(mul_gen)
        binPairsNum = next(mul_gen)
        itemsNum = next(mul_gen)
        binPairs = next(mul_gen)
        IDs = next(mul_gen)

        multi = open(self.multi, "w")
        print(f"Bin Length: {resolution}", file=multi)
        print(f"Total reads: {total}", file=multi)
        print(f"Num of reads: {IDsNum}", file=multi)
        print(f"Num of bin pairs: {binPairsNum}", file=multi)
        print(f"Num of items: {itemsNum}", file=multi)
        print(f"#bins", file=multi)
        pairsIndeces = {}
        i = 0
        for pair in sorted(binPairs, key=lambda x: (x[0], x[1])):
            pairsIndeces[pair] = i
            i += 1
            print(f"{binPairs[pair][0]} {pair[0]} {pair[1]} {binPairs[pair][1]} {binPairs[pair][2]}", file=multi)
        print("#reads ", end="", file=multi)
        print(multi.tell() - 7, file=multi)
        for ID in IDs:
            lst = IDs[ID]
            print(len(lst), file=multi)
            for i in range(len(lst)):
                chr1, bin1, chr2, bin2 = lst[i]
                ind = pairsIndeces[(min(bin1, bin2), max(bin1, bin2))]
                print(f"{ID} {chr1} {bin1} {chr2} {bin2} {ind}", file=multi)
        multi.close()


class Reader:
    def __init__(self, inFiles: list, resolution: int):
        self.inFiles = inFiles
        self.resolution = resolution

    @abstractmethod
    def __call__(self, collector: Type[DataCollecter]):
        pass


class mHiC_Reader(Reader):
    def __init__(self, inFiles: list, resolution: int):
        super().__init__(inFiles, resolution)
        assert len(self.inFiles) == 2, f"mHiC-reader: 2 file names needed, given {len(self.inFiles)}"
        self.multiName = self.inFiles[0]
        self.uniName = self.inFiles[1]

    def __call__(self, collecter: Type[DataCollecter]):
        multi = open(self.multiName, "r")
        uni = open(self.uniName, "r")

        for line in multi:
            ID, chr1, mid1, chr2, mid2 = line.split()
            bin1 = (int(mid1) - self.resolution // 2) // self.resolution
            bin2 = (int(mid2) - self.resolution // 2) // self.resolution
            collecter(bin1, bin2, chr1, chr2, readID=ID)
        for line in uni:
            chr1, mid1, chr2, mid2, count = line.split()
            bin1 = (int(mid1) - self.resolution // 2) // self.resolution
            bin2 = (int(mid2) - self.resolution // 2) // self.resolution
            collecter(bin1, bin2, chr1, chr2, contactCount=int(count))

        multi.close()
        uni.close()


class Parser:
    def __init__(self, reader: Type[Reader], collecter: Type[DataCollecter], writer: Type[Writer]):
        self.reader = reader
        self.collecter = collecter
        self.writer = writer

    def __call__(self):
        self.reader(self.collecter)
        self.writer(self.collecter)
        print("ok")


class Test:
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Multi2Uni-Parser')
    parser.add_argument("--input_format", type=str, nargs="?", const="mHiC", choices=["mHiC"], help="Data format of input file(s): \
                       mHiC (Multi and Uni reads files, resolution)")
    parser.add_argument("--input_files", "-i", required=True, type=str, nargs='+',
                        help="Input files. If input_format=mHiC, then arrange is MULTI_FILE UNI_FILE")
    parser.add_argument("--output_format", type=str, nargs="?", const="2files", choices=["2files"],
                        help="Format of output files")
    parser.add_argument("--output_files", "-o", required=True, type=str, nargs='+',
                        help="Output files. <File to probality calculation> <2 files to Expectation Maximization algorithm>")
    parser.add_argument("--resolution", "-r", required=True, type=int, help="Resolution, lenght of a single bin")
    args = parser.parse_args()
    fileParser = Parser(mHiC_Reader(args.input_files, args.resolution), StandardRNACollecter(args.resolution),
                        RNA_Writer(args.output_files))
    fileParser()
