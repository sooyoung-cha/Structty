#include "Panel.hpp"

Panel::Panel(int width, const std::string& mode) : panel_width(width), panel_mode(mode) {}

void Panel::add_panel_info(const std::string& file_name, const std::map<char, int>& chain_info, const std::map<char, int>& chain_residue_info) {
    panel_atom_info[file_name] = chain_info;
    panel_residue_info[file_name] = chain_residue_info;
}

int Panel::get_height() const {
    int lines = 0;
    lines += 3; 

    for (const auto& [file_name, chain_info] : panel_atom_info) {
        lines += 1;

        int n = (int)chain_info.size();
        int chain_lines = (n == 0) ? 1 : ((n + 2) / 3); 
        lines += chain_lines;

        lines += 1;
    }

    lines += 1;
    return lines;
}

void Panel::draw_panel(int start_row, int start_col,
                       int max_rows, int max_cols) const
{
    int term_rows, term_cols;
    getmaxyx(stdscr, term_rows, term_cols);

    const int top    = std::max(0, start_row);
    const int left   = std::max(0, start_col);
    const int bottom = std::min(term_rows, top + std::max(0, max_rows));
    const int right  = std::min(term_cols, left + std::max(0, max_cols));

    if (top >= bottom || left >= right) return;

    auto in_rows = [&](int r) { return r >= top && r < bottom; };
    auto avail_w = [&]() { return right - left; };

    int r = top;
    const int c = left;

    const int num_colors = (int)(sizeof(Palettes::UNRAINBOW) / sizeof(int));

    auto draw_border_line = [&](int rr) {
        if (!in_rows(rr)) return;
        const int w = avail_w();
        if (w <= 0) return;

        move(rr, c);
        addch('*');

        int mid = std::max(0, w - 2);
        for (int i = 0; i < mid; ++i) addch('=');

        if (w >= 2) addch('*');
    };

    auto clear_line = [&](int rr) {
        if (!in_rows(rr)) return;
        move(rr, c);
        for (int i = 0; i < avail_w(); ++i) addch(' ');
        move(rr, c);
    };

    draw_border_line(r);
    ++r;
    if (r >= bottom) return;

    {
        clear_line(r);
        const std::string help =
            "\tW A S D : ^ < v >\tR F : Zoom In/Out\tX Y Z : Rotate X, Y, Z axis\tC : Screenshot\tQ : Quit";
        mvaddnstr(r, c, help.c_str(), avail_w());
        ++r;
        if (r >= bottom) return;
    }

    {
        if (!in_rows(r)) return;
        move(r, c);
        for (int i = 0; i < avail_w(); ++i) addch('-');
        ++r;
        if (r >= bottom) return;
    }

    // body
    for (const auto& [file_name, chain_info] : panel_atom_info) {
        if (r >= bottom) break;

        clear_line(r);
        mvaddnstr(r, c, file_name.c_str(), avail_w());
        ++r;
        if (r >= bottom) break;

        clear_line(r);
        mvaddch(r, c, '\t');

        int count = 0;
        for (const auto& [chainID, length] : chain_info) {
            if (r >= bottom) break;

            if (count > 0 && count % 3 == 0) {
                ++r;
                if (r >= bottom) break;
                clear_line(r);
                mvaddch(r, c, '\t');
            }

            int residue_cnt = 0;
            auto itF = panel_residue_info.find(file_name);
            if (itF != panel_residue_info.end()) {
                auto itC = itF->second.find(chainID);
                if (itC != itF->second.end()) residue_cnt = itC->second;
            }

            char buf[64];
            std::snprintf(buf, sizeof(buf), "%c: %d (%d)\t", chainID, residue_cnt, length);

            int x, y;
            getyx(stdscr, y, x);
            int remain = right - x;
            if (remain <= 0) { ++count; continue; }

            if (panel_mode == "chain") {
                int pair = Palettes::UNRAINBOW[count % num_colors];
                attron(COLOR_PAIR(pair));
                addnstr(buf, remain);
                attroff(COLOR_PAIR(pair));
            } else {
                addnstr(buf, remain);
            }

            ++count;
        }

        ++r;
        if (r >= bottom) break;

        clear_line(r);
        ++r;
    }

    if (r >= bottom) return;

    draw_border_line(r);
}
