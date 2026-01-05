#include "Screen.hpp"

const float FOV = 90.0;
const float PI = 3.14159265359f;
const int MAX_STRUCT_NUM = 6;

Screen::Screen(const int& width, const int& height, const bool& show_structure, const std::string& mode, const std::string& depthcharacter) {
    screen_width = width;
    screen_height = height;
    screen_show_structure = show_structure;
    screen_mode = mode;
    screen_depthcharacter = depthcharacter;
    aspect_ratio = (float)screen_width / screen_height;
    zoom_level = 2.0f;

    start_color();
    use_default_colors();    
    init_color_pairs();
    
    camera = new Camera(width, height, mode);
    panel = new Panel(width, mode);
}

Screen::~Screen() {
    for (Protein* p : data) {
        delete p;
    }
    data.clear(); 
    
    delete camera;
    delete panel;
}

void Screen::init_color_pairs() {
    if (screen_mode == "protein") {
        for (int i = 0; i < data.size(); ++i) {
            init_pair(i+1, Palettes::UNRAINBOW[i], -1);
        }
    }

    else if (screen_mode == "chain") {
        int num_colors = sizeof(Palettes::UNRAINBOW) / sizeof(int);

        for (int i = 0; i < num_colors; ++i) {
            init_pair(i+1, Palettes::UNRAINBOW[i], -1);
        }
    }

    else if (screen_mode == "rainbow") {
        int num_colors = Palettes::RAINBOW.size();

        for (int i = 0; i < num_colors; ++i) {
            init_pair(i + 1, Palettes::RAINBOW[i], -1);
        }
    }
}

void Screen::set_protein(const std::string& in_file, int ii, const bool& show_structure) {
    Protein* protein = new Protein(in_file, chainVec.at(ii), show_structure);
    data.push_back(protein);
    pan_x.push_back(0.0f);
    pan_y.push_back(0.0f);
    
    init_color_pairs();
}

void Screen::set_tmatrix() {
    size_t filenum = data.size();
    vectorpointer = new float* [filenum];
    for (int i = 0; i < filenum; i++) {
        vectorpointer[i] = new float[3];
        vectorpointer[i][0] = 0;
        vectorpointer[i][1] = 0;
        vectorpointer[i][2] = 0;
    }
}
void Screen::set_chainfile(const std::string& chainfile, int filesize) {
    for (size_t i= 0; i < filesize; i++){
        chainVec.push_back("-");
    }
    if(chainfile == "") {
        return;
    }
    std::ifstream file(chainfile);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open chainfile\n";
        return;
    }
    // 0 1 2 3 4 5
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        int index;
        std::string chainlist;

        iss >> index;
        iss >> chainlist;
        if (index >= filesize) {
            std::cout << "Index in utmatrix file should be from 0 to yourfilenum - 1." << std::endl;
            continue;
        }
        // std::istringstream chainlistss(chainlist);
        chainVec[index]= chainlist;
    }
    file.close();
}

void Screen::set_utmatrix(const std::string& utmatrix, bool applyUT) {
    yesUT = !utmatrix.empty();

    const size_t filenum = data.size();
    float** matrixpointer = new float*[filenum];
    for (size_t i = 0; i < filenum; i++) {
        matrixpointer[i] = new float[9];
        for (int j = 0; j < 9; j++) {
            matrixpointer[i][j] = (j % 4 == 0) ? 1.f : 0.f; // identity
        }
    }

    if (utmatrix.empty()) {
        return;
    }

    std::ifstream file(utmatrix);
    if (!file.is_open()) {
        std::cerr << "Failed to open utmatrix file\n";
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        int index;
        std::string mat9Str;
        std::string mat3Str;

        iss >> index >> mat9Str >> mat3Str;
        if (index < 0 || index >= (int)filenum) continue;

        {
            std::istringstream mss(mat9Str);
            std::string val;
            int count = 0;
            while (std::getline(mss, val, ',') && count < 9) {
                matrixpointer[index][count++] = std::stof(val);
            }
        }

        {
            std::istringstream mss(mat3Str);
            std::string val;
            int count = 0;
            while (std::getline(mss, val, ',') && count < 3) {
                vectorpointer[index][count++] = std::stof(val);
            }
        }
    }

    if (applyUT) {
        for (size_t i = 0; i < filenum; i++) {
            data[i]->do_naive_rotation(matrixpointer[i]); // U
            data[i]->do_shift(vectorpointer[i]);          // T
        }
    }
    for (size_t i = 0; i < filenum; i++) delete[] matrixpointer[i];
    delete[] matrixpointer;
}

void Screen::normalize_proteins(const std::string& utmatrix) {
    const bool hasUT = !utmatrix.empty();
    for (size_t i = 0; i < data.size(); i++) {
        auto* p = data[i];
        p->load_data(vectorpointer[i], yesUT);
        panel->add_panel_info(p->get_file_name(),
                              p->get_chain_length(),
                              p->get_residue_count());
    }

    // 2) ut
    if (hasUT) {
        set_utmatrix(utmatrix, true);
    }

    // 3) bounding box calculation
    global_bb = BoundingBox();
    for (auto* p : data) {
        p->set_bounding_box();
        global_bb = global_bb + p->get_bounding_box();
    }

    // 4) scale calculation, every structure in the screen
    float max_ext = std::max(global_bb.max_x - global_bb.min_x,
                             global_bb.max_y - global_bb.min_y);
    max_ext = std::max(max_ext, global_bb.max_z - global_bb.min_z);
    float scale = (max_ext > 0.f) ? (2.0f / max_ext) : 1.0f;
    if (hasUT) {
        // relative position save

        float gx = 0.5f * (global_bb.min_x + global_bb.max_x);
        float gy = 0.5f * (global_bb.min_y + global_bb.max_y);
        float gz = 0.5f * (global_bb.min_z + global_bb.max_z);

        float global_shift[3] = { -gx, -gy, -gz };

        for (auto* p : data) {
            p->set_scale(scale);
            p->do_shift(global_shift); //same shift for every structure
            p->do_scale(scale);
        }

    } else {
        // If ut not input, go center

        for (auto* p : data) {
            float center_shift[3] = { -p->cx, -p->cy, -p->cz };

            p->set_scale(scale);
            p->do_shift(center_shift);
            p->do_scale(scale);
        }
    }
}



char Screen::get_pixel_char_from_depth(float z, float min_z, float max_z) {
    z -= focal_offset;
    float zn = (z - min_z) / (max_z - min_z);

    if (zn < 0.08f) return screen_depthcharacter[0];
    else if (zn < 0.18f) return screen_depthcharacter[1];
    else if (zn < 0.32f) return screen_depthcharacter[2];
    else if (zn < 0.50f) return screen_depthcharacter[3];
    else if (zn < 0.70f) return screen_depthcharacter[4];
    else if (zn < 0.85f) return screen_depthcharacter[5];
    else return screen_depthcharacter[6];
}

void Screen::draw_line(std::vector<RenderPoint>& points,
                      int x1, int x2, 
                      int y1, int y2,
                      float z1, float z2, 
                      std::string chainID, char structure,
                      float min_z, float max_z) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int dz = z2 - z1;
    int steps = std::max(abs(dx), abs(dy));
    if (steps == 0) steps = 1;

    float xIncrement = (float)dx / steps;
    float yIncrement = (float)dy / steps;
    float zIncrement = (float)dz / steps;

    float x = x1;
    float y = y1;
    float z = z1;

    for (int i = 0; i <= steps; ++i) {
        int ix = static_cast<int>(x);
        int iy = static_cast<int>(y);

        if (ix >= 0 && ix < screen_width && iy >= 0 && iy < screen_height) {
            points.push_back({ix, iy, z, get_pixel_char_from_depth(z, min_z, max_z), 0, chainID, structure});
        }

        x += xIncrement;
        y += yIncrement;
        z += zIncrement;
    }
}

void Screen::assign_colors_to_points(std::vector<RenderPoint>& points, int protein_idx) {
    if (screen_mode == "protein") {
        int num_colors = sizeof(Palettes::UNRAINBOW) / sizeof(int);
        int color_idx;

        for (auto& pt : points) {
            pt.color_id = (protein_idx % num_colors) + 1;
        }
    }

    else if (screen_mode == "chain") {
        int num_colors = sizeof(Palettes::UNRAINBOW) / sizeof(int);
        std::string cur_chain = points[0].chainID;
        int color_idx = 0;

        for (auto& pt : points) {
            std::string cID = pt.chainID;
            if (cID != cur_chain) {
                color_idx++;
                cur_chain = cID;
            }
            pt.color_id = (protein_idx * 10 + color_idx % num_colors) + 1;
        }
    }

    else if (screen_mode == "rainbow") {
        int num_points = points.size();
        int num_colors = Palettes::RAINBOW.size();

        for (int i = 0; i < num_points; i++) {
            int color_idx = (i * num_colors) / num_points;
            points[i].color_id = color_idx + 1;
        }
    }

    else {
        std::cerr << "Unknown mode: " << screen_mode << std::endl;
    }
}

void Screen::project() {
    float fovRads = 1.0 / tan((FOV / zoom_level) * 0.5 / 180.0 * PI);

    std::vector<RenderPoint> finalPoints;
    std::vector<RenderPoint> chainPoints;

    // project dots and connect them into line
    int protein_idx = 0;
    for (size_t ii = 0; ii < data.size(); ii++) {
        Protein* target = data[ii];
        chainPoints.clear();
        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.size() == 0) continue;

            int num_atoms = target->get_chain_length(chainID);
            int prevScreenX = -1, prevScreenY = -1;
            float prevZ = -1.0f;

            for (int i = 0; i < num_atoms; ++i) {
                float* position = chain_atoms[i].get_position();
                float x = position[0];
                float y = position[1];
                float z = position[2] + focal_offset;
                char structure = chain_atoms[i].get_structure();

                float projectedX = (x / z) * fovRads + pan_x[ii];
                float projectedY = (y / z) * fovRads + pan_y[ii];
                int screenX = (int)((projectedX + 1.0) * 0.5 * screen_width);
                int screenY = (int)((1.0 - projectedY) * 0.5 * screen_height);

                if (prevScreenX != -1 && prevScreenY != -1) {
                    draw_line(chainPoints, 
                             prevScreenX, screenX, 
                             prevScreenY, screenY, 
                             prevZ, z, 
                             chainID, structure, 
                             target->get_scaled_min_z(), target->get_scaled_max_z());
                }
                
                if (screenX >= 0 && screenX < screen_width && screenY >= 0 && screenY < screen_height) {
                    chainPoints.push_back({screenX, screenY, z, get_pixel_char_from_depth(z, target->get_scaled_min_z(), target->get_scaled_max_z()), 0, chainID, structure});
                }
                prevScreenX = screenX;
                prevScreenY = screenY;
                prevZ = z;
            }
        }
        assign_colors_to_points(chainPoints, protein_idx); 
        finalPoints.insert(finalPoints.end(), chainPoints.begin(), chainPoints.end());
        protein_idx++;

        for (const auto& pt : finalPoints) {
            int idx = pt.y * screen_width + pt.x;
            if (pt.depth < screenPixels[idx].depth) {
                screenPixels[idx].depth = pt.depth;
                screenPixels[idx].pixel = pt.pixel;
                screenPixels[idx].color_id = pt.color_id;
            }
        }
    }
}

void Screen::project(std::vector<RenderPoint>& projectPixels, const int proj_width, const int proj_height) {
    float fovRads = 1.0 / tan((FOV / zoom_level) * 0.5 / 180.0 * PI);

    std::vector<RenderPoint> finalPoints;
    std::vector<RenderPoint> chainPoints;

    // project dots and connect them into line
    int protein_idx = 0;
    for (size_t ii = 0; ii < data.size(); ii++) {
        Protein* target = data[ii];
        chainPoints.clear();
        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.size() == 0) continue;

            int num_atoms = target->get_chain_length(chainID);
            int prevScreenX = -1, prevScreenY = -1;
            float prevZ = -1.0f;

            for (int i = 0; i < num_atoms; ++i) {
                float* position = chain_atoms[i].get_position();
                float x = position[0];
                float y = position[1];
                float z = position[2] + focal_offset;
                char structure = chain_atoms[i].get_structure();

                float projectedX = (x / z) * fovRads + pan_x[ii];
                float projectedY = (y / z) * fovRads + pan_y[ii];
                int screenX = (int)((projectedX + 1.0) * 0.5 * proj_width);
                int screenY = (int)((1.0 - projectedY) * 0.5 * proj_height);

                if (prevScreenX != -1 && prevScreenY != -1) {
                    draw_line(chainPoints, 
                             prevScreenX, screenX, 
                             prevScreenY, screenY, 
                             prevZ, z, 
                             chainID, structure, 
                             target->get_scaled_min_z(), target->get_scaled_max_z());
                }
                
                if (screenX >= 0 && screenX < proj_width && screenY >= 0 && screenY < proj_height) {
                    chainPoints.push_back({screenX, screenY, z, get_pixel_char_from_depth(z, target->get_scaled_min_z(), target->get_scaled_max_z()), 0, chainID, structure});
                }
                prevScreenX = screenX;
                prevScreenY = screenY;
                prevZ = z;
            }
        }
        assign_colors_to_points(chainPoints, protein_idx); 
        finalPoints.insert(finalPoints.end(), chainPoints.begin(), chainPoints.end());
        protein_idx++;

        for (const auto& pt : finalPoints) {
            int idx = pt.y * proj_width + pt.x;
            if (pt.depth < projectPixels[idx].depth) {
                projectPixels[idx].depth = pt.depth;
                projectPixels[idx].pixel = pt.pixel;
                projectPixels[idx].color_id = pt.color_id;
            }
        }
    }
}
void Screen::clear_screen() {
    screenPixels.assign(screen_width * screen_height, RenderPoint());
}

void Screen::draw_screen() {
    clear_screen();
    project();

    int rows, cols;
    getmaxyx(stdscr, rows, cols);

    int panel_h = panel->get_height();
    if (panel_h > rows) panel_h = rows;

    // int margin = 3;
    int margin = 0;
    int offset = panel_h + margin;
    if (offset > rows) offset = rows;

    clear();

    print_screen(offset);

    int start_row = rows - panel_h;
    if (start_row < 0) start_row = 0;

    for (int r = start_row; r < rows; ++r) {
        move(r, 0);
        clrtoeol();
    }

    panel->draw_panel(start_row, 0, panel_h, cols);

    refresh();
}


void Screen::print_screen(int y_offset) {
    int rows, cols;
    getmaxyx(stdscr, rows, cols);

    for (int i = 0; i < screen_height; ++i) {
        int row = i - (y_offset/2)-3;
        if (row < 0) continue;
        if (row >= rows) break;

        int max_width = std::min(screen_width, cols);

        for (int j = 0; j < max_width; ++j) {
            int idx = i * screen_width + j;
            const RenderPoint& px = screenPixels[idx];

            if (px.color_id > 0) {
                attron(COLOR_PAIR(px.color_id));
                mvaddch(row, j, px.pixel);
                attroff(COLOR_PAIR(px.color_id));
            } else {
                mvaddch(row, j, px.pixel);
            }
        }
    }
}

void Screen::set_zoom_level(float zoom){
    if ((zoom_level + zoom > 0.5)&&(zoom_level + zoom < 15)){
        zoom_level += zoom;
    }
}

bool Screen::handle_input(){
    bool keep_show = true;

    auto pan_step_x = 2.0f * 4.0f / screen_width;
    auto pan_step_y = 2.0f * 2.0f / screen_height;

    auto apply_pan = [&](int idx, float dx, float dy){
        if (idx < 0 || idx >= (int)pan_x.size()) return;
        pan_x[idx] += dx;
        pan_y[idx] += dy;
    };

    int key = getch();
    switch(key){
        // select protein
        case 48:
        case 49:
        case 50:
        case 51:
        case 52:
        case 53:
        case 54:
            if (key == 48){
                structNum = -1;
            }
            else if (key - 48 <= data.size()) {
                structNum = key - 49;
            }
            break;
        // A, a (minus x-axis)
        case 65:
        case 97:
            if (structNum != -1) apply_pan(structNum, -pan_step_x, 0.0f);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, -pan_step_x, 0.0f);
            break;
        // D, d (plus x-axis)
        case 68:
        case 100:
            if (structNum != -1) apply_pan(structNum, +pan_step_x, 0.0f);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, +pan_step_x, 0.0f);
            break;
        // S, s (minus y-axis)
        case 83:
        case 115:
            if (structNum != -1) apply_pan(structNum, 0.0f, -pan_step_y);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, 0.0f, -pan_step_y);
            break;    
        // W, w (plus y-axis)
        case 87:
        case 119:
            if (structNum != -1) apply_pan(structNum, 0.0f, +pan_step_y);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, 0.0f, +pan_step_y);
            break;

        // X, x (rotate x-centered)
        case 88:
        case 120:
            if (structNum != -1) {
                data[structNum]->set_rotate(1, 0, 0);
            } else {
                for (int i = 0; i < data.size(); i++){
                    data[i]->set_rotate(1, 0, 0);
                }
            }
            break;  
        // Y, y (rotate y-centered)
        case 89:
        case 121:
            if (structNum != -1) {
                data[structNum]->set_rotate(0, 1, 0);
            } else {
                for (int i = 0; i < data.size(); i++){
                    data[i]->set_rotate(0, 1, 0);
                }
            }
            break;  
        // Z, z (rotate z-centered)
        case 90:
        case 122:
            if (structNum != -1) {
                data[structNum]->set_rotate(0, 0, 1);
            } else {
                for (int i = 0; i < data.size(); i++){
                    data[i]->set_rotate(0, 0, 1);
                }
            }
            break;  

        // F, f (zoom out)
        case 70:
        case 102:
            set_zoom_level(-0.3);
            break;   
        // R, R (zoom in)
        case 82:
        case 114:
            set_zoom_level(0.3);
            break;   

        // C, c (camera)
        case 67:
        case 99:     
        {     
            std::vector<RenderPoint> screenshotPixels;
            screenshotPixels.assign(screen_width * screen_height, RenderPoint());
            project(screenshotPixels, screen_width, screen_height);
            camera->screenshot(screenshotPixels);
            // camera->screenshot(screenPixels);
            break;
        }
        // Q, q
        case 81:
        case 113:
            keep_show = false;
            break;

        default:
            break;       
    }     

    return keep_show;
}