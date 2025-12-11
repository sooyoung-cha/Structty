#pragma once
#include "Protein.hpp"
#include "Atom.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "RenderPoint.hpp"
#include "lodepng.h"
#include <ctime>
#include <thread>
#include <chrono>  
#include <vector>
#include <iostream>
#include <filesystem> 
#include <unordered_map>

class Camera {
public:
    Camera(const std::string dir, const int width, const int height, const std::string mode);
    ~Camera() = default;

    void screenshot(const std::vector<RenderPoint>& screenPixels);    
    void renderPoint2image(const std::vector<RenderPoint>& screenPixels,
                           std::vector<RGBA>& screenImage);
    bool save_image(std::vector<RGBA>& screenImage);
    int get_alpha_from_depth(float z, float min_z, float max_z);

private:
    float focal_offset = 10.0f;

    std::string camera_dir;
    int camera_width, camera_height;
    int screenshot_idx;

    int height_duplicate = 2;

    std::string camera_mode = "default";
};