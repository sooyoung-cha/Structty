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

class Screen {
public:
    Screen(const int& width, const int& height, const bool& show_structure, const std::string& mode, const std::string& depthcharacter);
    ~Screen();
    bool handle_input();
    char get_pixel_char_from_depth(float z, float min_z, float max_z);
    void set_protein(const std::string& in_file, int ii, const bool& show_structure);
    void normalize_proteins(const std::string& utmatrix);
    void set_tmatrix();
    void set_utmatrix(const std::string& utmatrix, bool onlyU);
    void set_chainfile(const std::string& chainfile, int filesize);
    void set_zoom_level(float zoom);
    void draw_screen();
    void init_color_pairs();
    void assign_colors_to_points(std::vector<RenderPoint>& points, int protein_idx);
    void draw_line(std::vector<RenderPoint>& points,
                  int x1, int x2, 
                  int y1, int y2,
                  float z1, float z2, 
                  char chainID, char structure,
                  float min_z, float max_z);

private:
    int screen_width, screen_height;
    float aspect_ratio;
    bool screen_show_structure;
    bool yesUT = false;
    std::string screen_mode;
    std::string screen_depthcharacter;
    int structNum = -1;
    
    float focal_offset = 5.0f;
    float zoom_level;
    std::vector<float> pan_x;
    std::vector<float> pan_y;
    std::vector<std::string> chainVec; 
    float ** vectorpointer;
    std::vector<RenderPoint> screenPixels;   
    std::vector<Protein*> data;  

    BoundingBox global_bb;
    Camera* camera;
    Panel* panel;

    void project();
    void project(std::vector<RenderPoint>& screenshotPixels, const int proj_width, const int proj_height);
    void clear_screen();
    void print_screen(int panel_lines);
};