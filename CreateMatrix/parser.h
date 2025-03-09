#include "reader.h"


// Parser class which integrates Reader, DataCollecter, and Writer.
class Parser {
    public:
        Reader &reader;
        DataCollecter &collecter;
        Writer &writer;
        Parser(Reader &reader, DataCollecter &collecter, Writer &writer)
            : reader(reader), collecter(collecter), writer(writer) {}
    
        void operator()() {
            reader(collecter);
            writer(collecter);
            std::cout << "ok" << std::endl;
        }
    };
    