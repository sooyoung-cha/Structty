#include "Parameters.hpp"
#include <cmath>

bool is_nonnegative_number(const char* s) {
    if (s == nullptr || *s == '\0') return false;
    char* end = nullptr;
    errno = 0;
    double val = std::strtod(s, &end);
    if (end == s || *end != '\0') return false;
    if (errno == ERANGE) return false;
    if (val <= 0.0) return false;
    if (std::isnan(val)) return false;
    return true;
}

void print_help(){
    std::cout<<"-m, --mode:\n\t1. protein (default)\n\t2. chain\n\t3. rainbow"<<std::endl;
    std::cout<<"-d, --depth:\n\t1 .#@%*^-. (default)\n\t2. 7-character user input e.g. -d a134((%"<<std::endl;
    std::cout<<"-c, --chains:\n\tshow only the selected chains, see example/chainfile"<<std::endl;
    std::cout<<"-w, --width\n\t1. 3 (default)\n\t2. User input above 0, below 2000"<<std::endl;
    std::cout<<"-h, --height\n\t1. 3 (default)\n\t2. User input above 0, below 2000"<<std::endl;
    std::cout<<"-s, --structure:\n\tshow secondary structure (alpha helix, beta sheet)"<<std::endl;
    std::cout<<"-p, --predict:\n\tshow secondary structure with prediction if it is not described in the input file"<<std::endl;
    std::cout<<"-ut, --utmatrix:\n\trotate and translate, see example/utfile"<<std::endl;
}
Parameters::Parameters(int argc, char* argv[]) {
    arg_okay = true;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--help")) {
            print_help();
            std::exit(0);
        }
    }
    
    if (argc <= 1) {
        std::cerr << "Need input file dir" << std::endl;
        arg_okay = false;
        return;
    }

    for (int i = 1; i < argc; i++) {
        try {
            if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mode")) {
                if (i + 1 < argc) {
                    std::string val(argv[i + 1]);
                    std::transform(val.begin(), val.end(), val.begin(), ::tolower); // to lowercase
                    if (val == "chain" || val == "rainbow" || val == "protein") {
                        mode = val;
                        i++;
                    } else {
                        throw std::runtime_error("Error: Invalid value for --mode. Use 'protein', 'chain' or 'rainbow'.");
                    }
                } else {
                    throw std::runtime_error("Error: Missing value for -m / --mode.");
                }
            } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--depth")) {
                if (i + 1 < argc) {
                    std::string val(argv[i + 1]);
                    // std::transform(val.begin(), val.end(), val.begin(), ::tolower); // to lowercase
                    if (val != "") {
                        depthcharacter = val;
                        i++;
                    }
                } else {
                    throw std::runtime_error("Error: Missing value for -d / --depth.");
                }
            }
            else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--chains")) {
                if (i + 1 < argc) {
                    chainfile = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -c / --chains.");
                }
                // if (i + 1 < argc) {  
                //     if (argv[i+1] == nullptr || strlen(argv[i+1]) == 0) {  // if empty value
                //         throw std::runtime_error("Error: Chains argument is empty.");
                //     }
                //     while(in_file.size() - 1 != chains.size()){
                //         chains.push_back("-");
                //     }
                //     chains.push_back(argv[i+1]);  
                //     i++;
                // } else {
                //     throw std::runtime_error("Error: Missing argument for -c / --chains.");
                // }
            }
            else if (!strcmp(argv[i], "-w") || !strcmp(argv[i], "--width")) {
                if (i + 1 < argc) {
                    if (is_nonnegative_number(argv[i + 1])) {
                        width = std::stoi(argv[i + 1]);
                        ++i; 
                    } else {
                        throw std::runtime_error("Error: Parameter must be above 0");
                    }
                } else {
                    throw std::runtime_error("Error: Missing value for -w / --width.");
                }
            }
            else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--height")) {
                if (i + 1 < argc) {
                    if (is_nonnegative_number(argv[i + 1])) {
                        height = std::stoi(argv[i + 1]);
                        ++i; 
                    } else {
                        throw std::runtime_error("Error: Parameter must be above 0.");
                    }
                } else {
                    throw std::runtime_error("Error: Missing value for -h / --height.");
                }
            }
            else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--structure")) {
                show_structure = true;
            }
            else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--predict")) {
                predict_structure = true;
            }
            else if (!strcmp(argv[i], "-ut") || !strcmp(argv[i], "--utmatrix")) {
                if (i + 1 < argc) {
                    utmatrix = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -ut / --utmatrix.");
                }
            } else if (fs::exists(argv[i]) && fs::is_regular_file(argv[i]) && in_file.size() < 6){
                in_file.push_back(argv[i]);
            }
            else {
                throw std::runtime_error("Error: Unknown parameter: " + std::string(argv[i]));
            }
        }       
        catch (const std::exception& e) {
            std::cerr << "Wrong input parameters: " << e.what() << std::endl;
            std::cerr << "Error at argument: " << argv[i] << std::endl;
            arg_okay = false;
            return;
        }
    }
    while(in_file.size() != chains.size()){
        chains.push_back("-");
    }

    if (in_file.size() == 0){
        std::cerr << "Error: Need input file dir" << std::endl;
        arg_okay = false;
        return;
    }
    return;
}

void Parameters::print_args() {
    cout << "Input parameters >> " << endl;
    cout << "  in_file: " << endl;
    for (int i = 0; i < in_file.size(); i++) {
        std::cout << "\t" << in_file[i] << ": " << chains[i] << '\n'; 
    }
    cout << "  mode: " << mode << endl;
    cout << "  depthcharacter: " << depthcharacter << endl;
    cout << "  width: " << width << endl;
    cout << "  height: " << height << endl;
    cout << "  utmatrix: " << utmatrix << endl;
    cout << "  chainfile: " << chainfile << endl;
    cout << "  show_structure: " << show_structure << endl;

    cout << "\n";
    return;
}