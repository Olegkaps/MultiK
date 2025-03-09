#include "parser.h"
#include "argumentParser.h"



int main(int argc, char* argv[]) {
    // Using a simple command-line parser.
    Args args = parseArgs(argc, argv);
    
    if (args.input_format == "mhic") {
        mHiC_Reader reader(args.input_files, args.resolution);
        StandardRNACollecter collecter(args.resolution);
        RNA_Writer writer(args.output_files);
        
        Parser fileParser(reader, collecter, writer);
        fileParser();
    } else if (args.input_format == "contacts") {
        contacts_Reader reader(args.input_files, args.resolution, args.annotation, args.chr_ids); //use different collecter and writter
        contactsRNACollecter collecter(args.resolution);
        contacts_Writer writer(args.output_files);
        
        Parser fileParser(reader, collecter, writer);
        fileParser();
    }

    return 0;
}
