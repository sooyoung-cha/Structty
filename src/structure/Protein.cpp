#include "Protein.hpp"

static gemmi::Structure read_structure(const std::string& path) {
    gemmi::Structure st = gemmi::read_structure_file(path);
    st.remove_empty_chains();
    st.merge_chain_parts();
    return st;
}

static inline bool chain_ok(const std::string& target, char cid) {
    return (target == "-" || target.find(cid) != std::string::npos);
}

Protein::Protein(const std::string& in_file_, const std::string& target_chains_, const bool& show_structure_) {
    structureMaker = StructureMaker();

    in_file = in_file_;
    target_chains = target_chains_;
    show_structure = show_structure_;
    cx = cy = cz = scale = 0.0;
}

Protein::~Protein() {
}

std::map<char, std::vector<Atom>>& Protein::get_atoms() {
    return screen_atoms;  
}

std::map<char, int> Protein::get_residue_count() {
    return chain_res_count;
}

std::map<char, int> Protein::get_chain_length() {
    std::map<char, int> result;
    for (const auto& [chainID, atoms] : init_atoms) {
        result[chainID] = atoms.size();
    }
    return result;
}

int Protein::get_chain_length(char chainID) {
    if (screen_atoms.count(chainID)) {
        return screen_atoms[chainID].size();
    }
    return 0;
}

int Protein::get_length() {
    int total_atoms = 0;
    for (const auto& [chainID, atoms] : init_atoms) {
        total_atoms += atoms.size();
    }
    return total_atoms;
}

float Protein::get_scaled_min_z() { 
    return (bounding_box.min_z - cz) * scale; 
}

float Protein::get_scaled_max_z() { 
    return (bounding_box.max_z - cz) * scale; 
}

BoundingBox& Protein::get_bounding_box() {
     return bounding_box; 
}

void Protein::set_scale(float scale_) { 
    scale = scale_;
    cx = 0.5f * (bounding_box.min_x + bounding_box.max_x);
    cy = 0.5f * (bounding_box.min_y + bounding_box.max_y);
    cz = 0.5f * (bounding_box.min_z + bounding_box.max_z);
    ssPredictor.set_scale(1.0f/scale);
}    

bool Protein::is_ss_in_file(const std::string& in_file) {
    try {
        gemmi::Structure st = read_structure(in_file);
        return (!st.helices.empty() || !st.sheets.empty());
    } catch (...) {
        return false;
    }
}

void Protein::set_bounding_box() {
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            bounding_box.min_x = std::min(bounding_box.min_x, atom.x);
            bounding_box.min_y = std::min(bounding_box.min_y, atom.y);
            bounding_box.min_z = std::min(bounding_box.min_z, atom.z);
            bounding_box.max_x = std::max(bounding_box.max_x, atom.x);
            bounding_box.max_y = std::max(bounding_box.max_y, atom.y);
            bounding_box.max_z = std::max(bounding_box.max_z, atom.z);
        }
    }
}         

void Protein::count_seqres(const std::string& file) {
    // std::cout << "  count SEQRES\n";
    chain_res_count.clear();

    try {
        gemmi::Structure st = read_structure(file);

        for (const gemmi::Entity& ent : st.entities) {
            int len = (int)ent.full_sequence.size();
            if (len <= 0) continue;

            for (const std::string& cname : ent.subchains) {
                if (cname.empty()) continue;
                char cid = cname[0];
                chain_res_count[cid] = len;
            }
        }
    } catch (...) {
        std::cerr << "Error: gemmi SEQRES read failed\n";
    }
}

void Protein::load_init_atoms(const std::string& in_file, 
                              const std::string& target_chains,
                              const std::vector<std::tuple<char, int, char, int, char>>& ss_info, 
                              float * vectorpointers , bool yesUT) {
    // std::cout << "  load atoms\n";
    init_atoms.clear();

    gemmi::Structure st = read_structure(in_file);
    gemmi::Model& model = st.first_model();

    for (gemmi::Chain& chain : model.chains) {
        char cid = chain.name.empty() ? '?' : chain.name[0];
        if (!chain_ok(target_chains, cid))
            continue;

        for (gemmi::Residue& res : chain.residues) {
            const gemmi::Atom* ca = res.get_ca();
            if (!ca) continue;
            if (!res.seqid.num.has_value()) continue;

            int resn = (int)res.seqid.num;
            float x = (float)ca->pos.x;
            float y = (float)ca->pos.y;
            float z = (float)ca->pos.z;

            Atom a(x, y, z);

            for (auto& t : ss_info) {
                char sc; int s, ec; int e; char type;
                std::tie(sc, s, ec, e, type) = t;

                if (cid == sc && resn >= s && resn <= e) {
                    a.set_structure(type);
                    break;
                }
            }

            init_atoms[cid].push_back(a);
        }
    }
}

void Protein::load_init_atoms(const std::string& in_file, 
                              const std::string& target_chains, float * vectorpointers, bool yesUT) {
    // std::cout << "  load atoms\n";
    init_atoms.clear();

    gemmi::Structure st = read_structure(in_file);
    gemmi::Model& model = st.first_model();

    for (gemmi::Chain& chain : model.chains) {
        char cid = chain.name.empty() ? '?' : chain.name[0];
        if (!chain_ok(target_chains, cid))
            continue;

        for (gemmi::Residue& res : chain.residues) {
            const gemmi::Atom* ca = res.get_ca();
            if (!ca) continue;

            float x = (float)ca->pos.x;
            float y = (float)ca->pos.y;
            float z = (float)ca->pos.z;

            Atom a(x, y, z, 'x'); 
            init_atoms[cid].push_back(a);
        }
    }
}

void Protein::load_ss_info(const std::string& in_file,
                               const std::string& target_chains,
                               std::vector<std::tuple<char,int,char,int,char>>& ss_info)
{
    // std::cout << "  load SS info\n";
    ss_info.clear();

    gemmi::Structure st = read_structure(in_file);

    // Helix → H
    for (const gemmi::Helix& h : st.helices) {
        auto beg = h.start;
        auto end = h.end;

        if (!beg.res_id.seqid.num.has_value() ||
            !end.res_id.seqid.num.has_value())
            continue;

        char bc = beg.chain_name.empty() ? '?' : beg.chain_name[0];
        char ec = end.chain_name.empty() ? '?' : end.chain_name[0];
        if (!chain_ok(target_chains, bc)) continue;

        int bs = (int)beg.res_id.seqid.num;
        int es = (int)end.res_id.seqid.num;

        ss_info.emplace_back(bc, bs, ec, es, 'H');
    }

    // Sheet → S
    for (const gemmi::Sheet& sheet : st.sheets) {
        for (const gemmi::Sheet::Strand& s : sheet.strands) {
            auto beg = s.start;
            auto end = s.end;

            if (!beg.res_id.seqid.num.has_value() ||
                !end.res_id.seqid.num.has_value())
                continue;

            char bc = beg.chain_name.empty() ? '?' : beg.chain_name[0];
            char ec = end.chain_name.empty() ? '?' : end.chain_name[0];
            if (!chain_ok(target_chains, bc)) continue;

            int bs = (int)beg.res_id.seqid.num;
            int es = (int)end.res_id.seqid.num;

            ss_info.emplace_back(bc, bs, ec, es, 'S');
        }
    }
}

std::ostream& operator<<(std::ostream& os, const std::tuple<char, int, char, int, char>& t) {
    os << "("
       << std::get<0>(t) << ", "
       << std::get<1>(t) << ", "
       << std::get<2>(t) << ", "
       << std::get<3>(t) << ", "
       << std::get<4>(t) << ")";
    return os;
}

void Protein::load_data(float * vectorpointers, bool yesUT) {    
    // pdb
    if (in_file.find(".pdb") != std::string::npos || in_file.find(".cif") != std::string::npos) {
        if (show_structure){
            if (is_ss_in_file(in_file)){
                std::vector<std::tuple<char, int, char, int, char>> ss_info;
                load_ss_info(in_file, target_chains, ss_info);
                load_init_atoms(in_file, target_chains, ss_info, vectorpointers, yesUT);
            }
            else{
                load_init_atoms(in_file, target_chains, vectorpointers, yesUT);
                ssPredictor.run(init_atoms);
            }
        }
        else{
            load_init_atoms(in_file, target_chains, vectorpointers, yesUT);
        }
        
        if (init_atoms.empty()) {
            std::cerr << "Error: input PDB file is empty." << std::endl;
            return;
        }
        
        if (show_structure){
            structureMaker.calculate_ss_points(init_atoms, screen_atoms);
        }
        else{ screen_atoms = init_atoms; }
        count_seqres(in_file);
    }

    // others
    else{
        std::cerr << "Error: input file format is not supported." << std::endl;
        return;
    }
    std::cout << std::endl;
}

void Protein::set_rotate(int x_rotate, int y_rotate, int z_rotate){
    const float PI = 3.14159265359;
    // const float UNIT = 12;
    const float UNIT = 48;

    if (x_rotate != 0) {
        float values[9] = {1, 0, 0,
                   0, cos(x_rotate * PI / UNIT), -sin(x_rotate * PI / UNIT), 
                   0, sin(x_rotate * PI / UNIT), cos(x_rotate * PI / UNIT)};

        float* rotate_mat = new float[9];
        for (int i = 0; i < 9; i++) {
            rotate_mat[i] = values[i];
        }
        do_rotation(rotate_mat);
    }
    else if (y_rotate != 0) {
        float values[9] = {cos(y_rotate * PI / UNIT), 0, sin(y_rotate * PI / UNIT),
                   0, 1, 0, 
                   -sin(y_rotate * PI / UNIT), 0, cos(y_rotate * PI / UNIT)};

        float* rotate_mat = new float[9];
        for (int i = 0; i < 9; i++) {
            rotate_mat[i] = values[i];
        }
        do_rotation(rotate_mat);
    }
    else if (z_rotate != 0) {
        float values[9] = {cos(z_rotate * PI / UNIT), -sin(z_rotate * PI / UNIT), 0,
                   sin(z_rotate * PI / UNIT), cos(z_rotate * PI / UNIT), 0, 
                   0, 0, 1};

        float* rotate_mat = new float[9];
        for (int i = 0; i < 9; i++) {
            rotate_mat[i] = values[i];
        }
        do_rotation(rotate_mat);
    }

}


void Protein::set_shift(float shift_x, float shift_y, float shift_z) { 
    float* shift_mat = new float[3];
    shift_mat[0] = shift_x;
    shift_mat[1] = shift_y;
    shift_mat[2] = shift_z;
    do_shift(shift_mat);
}

void Protein::do_naive_rotation(float * rotate_mat) {
    float avgx = 0;
    float avgy = 0;
    float avgz = 0;
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            avgx = atom.x;
            avgy = atom.y;
            avgz = atom.z;
            atom.x  = avgx * rotate_mat[0] + avgy * rotate_mat[1]+ avgz * rotate_mat[2];
            atom.y  = avgx * rotate_mat[3] + avgy * rotate_mat[4]+ avgz * rotate_mat[5];
            atom.z  = avgx * rotate_mat[6] + avgy * rotate_mat[7]+ avgz * rotate_mat[8];
        }
    }
}
void Protein::do_rotation(float * rotate_mat) {
    float avgx = 0;
    float avgy = 0;
    float avgz = 0;
    int num = 0;
    // for (auto& [chainID, chain_atoms] : screen_atoms) {
    //     for (Atom& atom : chain_atoms) {
    //         avgx += atom.x;
    //         avgy += atom.y;
    //         avgz += atom.z;
    //         num += 1;   
    //     }
    // }
    // avgx /= num;
    // avgy /= num;
    // avgz /= num;

    float minx,maxx,miny,maxy,minz,maxz;
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            if (num == 0) {
                minx = atom.x;
                maxx = atom.x;
                miny = atom.y;
                maxy = atom.y;
                minz = atom.z;
                maxz = atom.z;
            } else {
                if (atom.x <minx) {
                    minx = atom.x;
                }
                if (atom.x > maxx) {
                    maxx = atom.x;
                }
                if (atom.y <miny) {
                    miny = atom.y;
                }
                if (atom.y > maxy) {
                    maxy = atom.y;
                }
                if (atom.z <minz) {
                    minz = atom.z;
                }
                if (atom.z > maxz) {
                    maxz = atom.z;
                }
            }
            num++;
        }
    }
    avgx = (minx + maxx) /2;
    avgy = (miny + maxy) /2;
    avgz = (minz + maxz) /2;
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            float x = atom.x;
            float y = atom.y;
            float z = atom.z;
            atom.x = avgx + rotate_mat[0] * (x - avgx) + rotate_mat[1] * (y - avgy) + rotate_mat[2] * (z - avgz);
            atom.y = rotate_mat[3] * (x - avgx) + avgy + rotate_mat[4] * (y - avgy) + rotate_mat[5] * (z - avgz);
            atom.z = rotate_mat[6] * (x - avgx) + rotate_mat[7] * (y - avgy) +  avgz + rotate_mat[8] * (z - avgz);
        }
    }
}

void Protein::do_shift(float* shift_mat) {
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            atom.x += shift_mat[0];
            atom.y += shift_mat[1];
            atom.z += shift_mat[2];
        }
    }
}


void Protein::do_scale(float scale) {
    for (auto& [chainID, chain_atoms] : screen_atoms) {
        for (Atom& atom : chain_atoms) {
            atom.x *= scale;
            atom.y *= scale;
            atom.z *= scale;
        }
    }
}