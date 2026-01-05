#include "StructureMaker.hpp"

const float PI = 3.14159265359f;

StructureMaker::StructureMaker() {}
StructureMaker::~StructureMaker() {}

std::vector<std::vector<Atom>> StructureMaker::extract_helix_segments(const Atom* atoms, int num_atoms) {
    std::vector<std::vector<Atom>> helices;
    std::vector<Atom> current;

    for (int i = 0; i < num_atoms; ++i) {
        if (atoms[i].structure == 'H') {
            current.push_back(atoms[i]);
        } else {
            if (current.size() >= 4) {
                helices.push_back(std::move(current));
            }
            current.clear();
        }
    }

    if (current.size() >= 4) {
        helices.push_back(std::move(current));
    }

    return helices;
}


void StructureMaker::compute_helix_axis(const std::vector<Atom>& helix, float (&center)[3], float (&axis)[3]) {
    size_t N = helix.size();
    if (N == 0) return;

    float sum[3] = {0.0f, 0.0f, 0.0f};

    // calculate center coord
    for (const Atom& atom : helix) {
        sum[0] += atom.x;
        sum[1] += atom.y;
        sum[2] += atom.z;
    }

    center[0] = sum[0] / N;
    center[1] = sum[1] / N;
    center[2] = sum[2] / N;

    // calculate covariance matrix
    float cov[3][3] = {0};

    for (const Atom& atom : helix) {
        float dx = atom.x - center[0];
        float dy = atom.y - center[1];
        float dz = atom.z - center[2];

        cov[0][0] += dx * dx;
        cov[0][1] += dx * dy;
        cov[0][2] += dx * dz;

        cov[1][0] += dy * dx;
        cov[1][1] += dy * dy;
        cov[1][2] += dy * dz;

        cov[2][0] += dz * dx;
        cov[2][1] += dz * dy;
        cov[2][2] += dz * dz;
    }

    // main axis change by power iteration
    float vec[3] = {1.0f, 0.0f, 0.0f};  // init vector
    for (int iter = 0; iter < 10; ++iter) {
        float x = cov[0][0]*vec[0] + cov[0][1]*vec[1] + cov[0][2]*vec[2];
        float y = cov[1][0]*vec[0] + cov[1][1]*vec[1] + cov[1][2]*vec[2];
        float z = cov[2][0]*vec[0] + cov[2][1]*vec[1] + cov[2][2]*vec[2];

        float norm = std::sqrt(x*x + y*y + z*z);
        if (norm == 0.0f) break;

        vec[0] = x / norm;
        vec[1] = y / norm;
        vec[2] = z / norm;
    }

    axis[0] = vec[0];
    axis[1] = vec[1];
    axis[2] = vec[2];
}



void StructureMaker::calculate_ss_points(std::map<std::string, std::vector<Atom>>& init_atoms,
                                         std::map<std::string, std::vector<Atom>>& ss_atoms) {
    // std::cout << "  apply structure\n";
    ss_atoms.clear();

    for (auto& [chainID, atoms] : init_atoms) {
        std::vector<Atom>& output = ss_atoms[chainID];
        size_t i = 0;
        while (i < atoms.size()) {
            char s = atoms[i].structure;

            if (s == 'H') {
                // helix start: find sequential H
                size_t start = i;
                while (i < atoms.size() && atoms[i].structure == 'H') ++i;
                size_t end = i;

                if (end - start >= 4) {
                    auto segment = std::vector<Atom>(atoms.begin() + start, atoms.begin() + end);

                    float center[3], axis[3];
                    compute_helix_axis(segment, center, axis);

                    float dx = segment.back().x - segment.front().x;
                    float dy = segment.back().y - segment.front().y;
                    float dz = segment.back().z - segment.front().z;
                    float length = std::sqrt(dx * dx + dy * dy + dz * dz);

                    const int steps = std::min<int>(circle_steps, (end - start));     

                    float up[3] = {0, 0, 1};
                    if (std::abs(axis[2]) > 0.99f) { up[0] = 1; up[2] = 0; }

                    float n1[3] = {
                        axis[1]*up[2] - axis[2]*up[1],
                        axis[2]*up[0] - axis[0]*up[2],
                        axis[0]*up[1] - axis[1]*up[0]
                    };
                    float n1_norm = std::sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
                    for (int k = 0; k < 3; ++k) n1[k] /= n1_norm;

                    float n2[3] = {
                        axis[1]*n1[2] - axis[2]*n1[1],
                        axis[2]*n1[0] - axis[0]*n1[2],
                        axis[0]*n1[1] - axis[1]*n1[0]
                    };

                    for (int s = 0; s <= steps; ++s) {
                        float t = static_cast<float>(s) / steps;
                        float base[3] = {
                            center[0] + axis[0] * (t - 0.5f) * length,
                            center[1] + axis[1] * (t - 0.5f) * length,
                            center[2] + axis[2] * (t - 0.5f) * length,
                        };

                        for (int a = 0; a < circle_steps; ++a) {
                            float theta = 2 * PI * a / circle_steps;
                            float dx = std::cos(theta);
                            float dy = std::sin(theta);

                            float px = base[0] + radius * (dx * n1[0] + dy * n2[0]);
                            float py = base[1] + radius * (dx * n1[1] + dy * n2[1]);
                            float pz = base[2] + radius * (dx * n1[2] + dy * n2[2]);

                            output.emplace_back(px, py, pz, 'H');
                        }
                    }
                } else {
                    // if too short, ignore
                    i = end;
                }
            }

            else if (s == 'S' && i + 1 < atoms.size() && atoms[i + 1].structure == 'S') {
                const Atom& p1 = atoms[i];
                const Atom& p2 = atoms[i + 1];

                float dx = p2.x - p1.x;
                float dy = p2.y - p1.y;
                float dz = p2.z - p1.z;
                float len = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (len == 0) { i++; continue; }

                float axis[3] = { dx / len, dy / len, dz / len };  // direction vector
                float up[3] = { 0.0f, 0.0f, 1.0f };
                if (std::abs(axis[2]) > 0.99f) {
                    up[0] = 1.0f; up[2] = 0.0f;  // if almost similar to z-axis, replace
                }

                // n1: vector perpendicular to axis
                float n1[3] = {
                    axis[1]*up[2] - axis[2]*up[1],
                    axis[2]*up[0] - axis[0]*up[2],
                    axis[0]*up[1] - axis[1]*up[0]
                };
                float n1_norm = std::sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
                for (int j = 0; j < 3; ++j) n1[j] /= n1_norm;

                // ribbon width vector, move to direction n1
                int line_steps = std::max<int>(2, static_cast<int>(len / 0.05f));

                for (int step = -width; step <= width; ++step) {
                    float offset[3] = {
                        n1[0] * step * 0.05f,
                        n1[1] * step * 0.05f,
                        n1[2] * step * 0.05f
                    };

                    float x1 = p1.x + offset[0];
                    float y1 = p1.y + offset[1];
                    float z1 = p1.z + offset[2];

                    float x2 = p2.x + offset[0];
                    float y2 = p2.y + offset[1];
                    float z2 = p2.z + offset[2];

                    for (int t = 0; t <= line_steps; ++t) {
                        float f = static_cast<float>(t) / line_steps;
                        float x = x1 + f * (x2 - x1);
                        float y = y1 + f * (y2 - y1);
                        float z = z1 + f * (z2 - z1);

                        output.emplace_back(x, y, z, 'S');
                    }
                }
                i++;  // sheet: pair
            }

            else {
                // no structure, just add
                output.push_back(atoms[i]);
                ++i;
            }
        }
    }
}
