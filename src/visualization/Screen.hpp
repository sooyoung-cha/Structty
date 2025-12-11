#pragma once
#include "Protein.hpp"
#include "Atom.hpp"
#include "RenderPoint.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "Panel.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <ncurses.h>  
#include <cstdlib>

inline std::string get_home_dir() {
    if (const char* home = std::getenv("HOME")) return std::string(home);
    throw std::runtime_error("HOME is not set");
}

class Screen {
public:
    Screen(const int& width, const int& height, const bool& show_structure, const std::string& mode, const std::string& depthcharacter);
    ~Screen();
    bool handle_input();
    char get_pixel_char_from_depth(float z, float min_z, float max_z);
    void set_protein(const std::string& in_file, const std::string& target_chains, const bool& show_structure);
    void normalize_proteins(const std::string& utmatrix);
    void set_tmatrix();
    void set_utmatrix(const std::string& utmatrix, bool onlyU);
    void set_zoom_level(float zoom);
    void draw_screen();
    void assign_colors_to_points(std::vector<RenderPoint>& points, int protein_idx);
    void draw_line(std::vector<RenderPoint>& points,
                  int x1, int x2, 
                  int y1, int y2,
                  float z1, float z2, 
                  char chainID, char structure,
                  float min_z, float max_z);

private:
    float focal_offset = 10.0f;
    int structNum = -1;
    float camera_mul = 1.0;
    int screen_width, screen_height;
    bool screen_show_structure;
    bool yesUT = false;
    std::string screen_mode;
    std::string screen_depthcharacter;
    float aspect_ratio;
    std::vector<RenderPoint> screenPixels;   
    std::vector<Protein*> data;  
    std::vector<float> zoom_level;
    float ** vectorpointer;

    BoundingBox global_bb;
    Camera* camera;
    Panel* panel;

    std::unordered_map<char, int> chain_colors;

    void project();
    void project(std::vector<RenderPoint>& screenshotPixels, const int proj_width, const int proj_height);
    void clear_screen();
    void print_screen();
};