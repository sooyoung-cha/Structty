#pragma once
#include <limits>

struct RenderPoint {
    int   x = 0;
    int   y = 0;
    float depth = std::numeric_limits<float>::infinity();  
    char  pixel = ' ';                                 
    int   color_id = 0; 
    std::string  chainID = "";
    char  structure = 0;
};
