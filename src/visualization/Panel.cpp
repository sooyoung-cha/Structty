#include "Panel.hpp"

Panel::Panel(int width) : panel_width(width) {}

void Panel::add_panel_info(const std::string& file_name, const std::map<char, int>& chain_info, const std::map<char, int>& chain_residue_info) {
    panel_atom_info[file_name] = chain_info;
    panel_residue_info[file_name] = chain_residue_info;
}

std::string Panel::get_panel_info() const {
    std::ostringstream oss;
    oss << "\n";
    oss << "*" << std::string(panel_width - 2, '=') << "*\n";

    oss << "\tW A S D : ^ < v >\t";
    oss << "R F : Zoom In/Out\t";
    oss << "X Y Z : Rotate X, Y, Z axis\t";
    oss << "C : Screenshot\t";
    oss << "Q : Quit\n";
    oss <<  std::string(panel_width, '-') << "\n";

    for (const auto& [file_name, chain_info] : panel_atom_info) {
        oss << file_name << "\n\t";
        int count = 0;
        for (const auto& [chainID, length] : chain_info) {
            if (panel_residue_info.find(file_name) == panel_residue_info.end()) {
                std::cerr << "[Error] file_name not found: " << file_name << "\n";
            }

            if (panel_residue_info.at(file_name).find(chainID) == panel_residue_info.at(file_name).end()) {
                std::cerr << "[Error] chainID not found: " 
                        << chainID 
                        << " in file " << file_name << "\n";
            }
            oss << chainID << ": " << panel_residue_info.at(file_name).at(chainID) << " (" << length << ")\t";
            count++;
            if (count % 3 == 0)
                oss << "\n\t";
        }
        if (count % 3 != 0)
            oss << "\n";

        oss << "\n";
    }
    oss << "*" << std::string(panel_width - 2, '=') << "*\n";

    return oss.str();
}