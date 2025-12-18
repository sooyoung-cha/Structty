#pragma once
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

class Parameters{
    private:
        bool show_structure = false;
        bool predict_structure = false;
        bool arg_okay = true;
        int width = 3;
        int height = 3;
        vector<string> in_file;
        vector<string> chains;
        string utmatrix = "";
        string chainfile = "";
        string mode = "protein";
        string depthcharacter = "#@%*^-.";
    public:
        Parameters(int argc, char* argv[]);

        void print_args();

        bool is_valid_number(const std::string& str, int min, int max);

        // get, set
        vector<string>& get_in_file(){
            return in_file;
        }
        string get_in_file(int idx){
            if (idx < in_file.size()){
                return in_file[idx];
            }
            return "";
        }
        // vector<string>& get_chains(){
        //     return chains;
        // }
        // string get_chains(int idx){
        //     if (idx < chains.size()){
        //         return chains[idx];
        //     }
        //     return "";
        // }
        string get_chainfile(){
            return chainfile;
        }
        string get_utmatrix(){
            return utmatrix;
        }
        int get_width(){
            return width * 32 + 32;
        }
        int get_height(){
            return height * 8 + 24;
        }
        string get_mode(){
            return mode;
        }
        string get_depthcharacter(){
            return depthcharacter;
        }
        bool get_show_structure(){
            return show_structure;
        }
        bool check_arg_okay(){
            return arg_okay;
        }
};