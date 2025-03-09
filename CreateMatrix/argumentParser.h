
// Simple parser for command-line arguments.
struct Args {
    std::string input_format;
    std::vector<std::string> input_files;
    std::string output_format;
    std::vector<std::string> output_files;
    std::string annotation;
    std::string chr_ids;
    int resolution;
};

Args parseArgs(int argc, char* argv[]) {
    Args args;
    // Defaults for optional arguments.
    args.annotation = "";
    args.chr_ids = "";
    args.input_format = "mhic";
    args.output_format = "";
    args.resolution = 0;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--annotation" || arg == "-a") {
            if (i + 1 < argc) {
                args.annotation = argv[++i];
            }
        } else if (arg == "--chr_ids") {
            if (i + 1 < argc) {
                args.chr_ids = argv[++i];
            }
        } else if (arg == "--input_format") {
            if (i + 1 < argc) {
                args.input_format = argv[++i];
            }
        } else if (arg == "--input_files" || arg == "-i") {
            // Read all subsequent non-flag arguments until next flag.
            while (i + 1 < argc && argv[i+1][0] != '-') {
                args.input_files.push_back(argv[++i]);
            }
        } else if (arg == "--output_format") {
            if (i + 1 < argc) {
                args.output_format = argv[++i];
            }
        } else if (arg == "--output_files" || arg == "-o") {
            while (i + 1 < argc && argv[i+1][0] != '-') {
                args.output_files.push_back(argv[++i]);
            }
        } else if (arg == "--resolution" || arg == "-r") {
            if (i + 1 < argc) {
                args.resolution = std::stoi(argv[++i]);
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "CreateMatrix help" << std::endl;
            std::cout << "-i\t--input_files\tREQIURED: list of input files (depences on format type)" << std::endl;
            std::cout << "-o\t--output_files\tREQIURED: list of output files (depences on format type)" << std::endl;
            std::cout << "-a\t--annotation\tgff RNA annotation file (for contacts format only)" << std::endl;
            std::cout << "--chr_ids\tfile comparing chromosome id and chromosome name" << std::endl;
            std::cout << "-\t\t\t'ContactsFile.mk' file will be generated - mapping numbers to names" << std::endl;
            std::cout << "--input_format\tinput files format:" << std::endl;
            std::cout << "\t\tmhic" << std::endl;
            std::cout << "\t\t\t2 input files: MultiContactsFile UniContactsFile" << std::endl;
            std::cout << "\t\t\t1 output file: MatrixFile" << std::endl;
            std::cout << "\t\tcontacts" << std::endl;
            std::cout << "\t\t\t1 input file: ContactsFile" << std::endl;
            std::cout << "\t\t\t2 output files: MatrixFile UniContactsFile" << std::endl;
        }
    }
    // Basic check for required arguments.
    if (args.input_files.empty() || args.output_files.empty() || args.resolution == 0) {
        std::cerr << "Required arguments missing." << std::endl;
        exit(1);
    }
    if ((args.annotation != "") != (args.chr_ids != "")) {
        std::cerr << "You must use both --annotation and --chr_ids" << std::endl;
        exit(1);
    }
    if (args.input_format != "mhic" && args.input_format != "contacts") {
        std::cerr << "Invalid input file format" << std::endl;
    }
    return args;
}
