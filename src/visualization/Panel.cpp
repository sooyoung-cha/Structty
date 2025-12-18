#include "Panel.hpp"

Panel::Panel(int width, const std::string& mode) : panel_width(width), panel_mode(mode) {}

void Panel::add_panel_info(const std::string& file_name, 
                           const std::map<char, int>& chain_info, 
                           const std::map<char, int>& chain_residue_info) {
    entries.push_back(Entry{
        file_name,
        chain_info,
        chain_residue_info
    });
}

int Panel::get_height() const {
    int lines = 0;
    lines += 3; 
    for (const auto& entry : entries) {
        lines += 1; 
        int n = (int)entry.chain_atom_info.size();
        int chain_lines = (n == 0) ? 1 : ((n + 2) / 3); // 3 per line
        lines += chain_lines;
        lines += 1; 
    }
    lines += 1; 
    return lines;
}

void Panel::draw_panel(int start_row, int start_col,
                       int max_rows, int max_cols) const {
    const int num_colors = (int)(sizeof(Palettes::UNRAINBOW) / sizeof(int));
    if (max_rows <= 0 || max_cols <= 0) return;

    // 출력 가능 영역(클리핑)
    const int top    = start_row;
    const int left   = start_col;
    const int bottom = start_row + max_rows; // exclusive
    const int right  = start_col + max_cols; // exclusive

    auto in_rows = [&](int rr){ return rr >= top && rr < bottom; };
    auto remain_cols = [&](int x){ return right - x; };

    auto clear_line = [&](int rr){
        if (!in_rows(rr)) return;
        move(rr, left);
        clrtoeol();
        move(rr, left);
    };

    auto put_str = [&](int& rr, int& x, const std::string& s){
        if (!in_rows(rr)) return;
        int rem = remain_cols(x);
        if (rem <= 0) return;
        addnstr(s.c_str(), rem);
        x += (int)s.size();
    };

    auto put_cstr = [&](int& rr, int& x, const char* s){
        if (!in_rows(rr)) return;
        int rem = remain_cols(x);
        if (rem <= 0) return;
        addnstr(s, rem);
    };

    int r = start_row;

    // Top border
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "*");
        int mid = std::max(0, std::min(panel_width - 2, max_cols - 2));
        put_str(r, x, std::string(mid, '='));
        if (max_cols >= 2) put_str(r, x, "*");
    }
    ++r;
    if (!in_rows(r)) return;

    // Help line
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "\tW A S D : ^ < v >\t");
        put_str(r, x, "R F : Zoom In/Out\t");
        put_str(r, x, "X Y Z : Rotate X, Y, Z axis\t");
        put_str(r, x, "C : Screenshot\t");
        put_str(r, x, "Q : Quit");
    }
    ++r;
    if (!in_rows(r)) return;

    // Separator
    clear_line(r);
    {
        move(r, left);
        int w = std::min(panel_width, max_cols);
        for (int i = 0; i < w; ++i) addch('-');
    }
    ++r;
    if (!in_rows(r)) return;

    // Body
    int file_idx = 0;
    for (const auto& entry : entries) {
        if (!in_rows(r)) break;

        const std::string& file_name = entry.file_name;
        const auto& chain_info       = entry.chain_atom_info;

        // protein 모드: 파일 단위 색
        int protein_pair = 0;
        if (panel_mode == "protein" && num_colors > 0) {
            protein_pair = (file_idx % num_colors) + 1; // 1..num_colors
        }

        // file name line
        clear_line(r);
        {
            int x = left;
            if (protein_pair > 0) attron(COLOR_PAIR(protein_pair));
            put_str(r, x, file_name);
            if (protein_pair > 0) attroff(COLOR_PAIR(protein_pair));
        }
        ++r;
        if (!in_rows(r)) break;

        // chain lines
        clear_line(r);
        move(r, left);
        addch('\t');
        int x = left + 1;

        int count = 0;
        for (const auto& [chainID, length] : chain_info) {
            if (!in_rows(r)) break;

            // 3개마다 줄바꿈 + 탭
            if (count > 0 && count % 3 == 0) {
                ++r;
                if (!in_rows(r)) break;
                clear_line(r);
                move(r, left);
                addch('\t');
                x = left + 1;
            }

            // residue count (Entry에서 안전 접근)
            int residue_cnt = 0;
            auto itC = entry.chain_residue_info.find(chainID);
            if (itC != entry.chain_residue_info.end()) residue_cnt = itC->second;

            char buf[64];
            std::snprintf(buf, sizeof(buf), "%c: %d (%d)\t", chainID, residue_cnt, length);

            // chain 모드: 체인 단위 색
            int chain_pair = 0;
            if (panel_mode == "chain" && num_colors > 0) {
                chain_pair = (count % num_colors) + 1; // 1..num_colors
            }

            // protein 모드면 protein_pair 우선
            int pair_to_use = (panel_mode == "protein") ? protein_pair : chain_pair;

            if (pair_to_use > 0) attron(COLOR_PAIR(pair_to_use));
            put_cstr(r, x, buf);
            if (pair_to_use > 0) attroff(COLOR_PAIR(pair_to_use));

            ++count;
        }

        // blank line
        ++r;
        if (!in_rows(r)) break;
        clear_line(r);
        ++r;

        ++file_idx;
    }

    if (!in_rows(r)) return;

    // Bottom border
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "*");
        int mid = std::max(0, std::min(panel_width - 2, max_cols - 2));
        put_str(r, x, std::string(mid, '='));
        if (max_cols >= 2) put_str(r, x, "*");
    }
}
