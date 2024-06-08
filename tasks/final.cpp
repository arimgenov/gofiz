#include <cmath>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <string>

const int step = 10;
const int N = 1000, K_x = 100, K_y = 100;
const int P_x = 5, P_y = 5; // пространство скоростей
const int a = 10, c_x = 20, c_y = 50;
const int p0_x = 1, p0_y = 0; // начальные скорости
const double dv_x = 5.0 / P_x, dv_y = 5.0 / P_y;

const int wallStart = 50, wallEnd = 52; // координаты стены
const int holeStart = 48, holeEnd = 52;

// geometry
bool* skipMaskX;
bool* skipMaskY;

struct Segment {int coordinate, start, end;};

const Segment wallRightX[] = {
        {0,             0,               K_y},
        {wallStart + 1, 0,             holeStart + 1},
        {wallStart + 1, holeEnd + 1,     K_y},
        {wallEnd + 1,   0,               holeStart + 1},
        {wallEnd + 1,   holeEnd + 1,     K_y}
};
const Segment wallLeftX[] = {
        {wallStart, 0,             holeStart + 1},
        {wallStart, holeEnd + 1,   K_y},
        {wallEnd,   0,               holeStart + 1},
        {wallEnd,   holeEnd + 1, K_y},
//        {K_x - 1,   0, K_y}
};
const Segment wallUpY[] = {
        {0,             0,             wallEnd + 1},
        {holeStart + 1, wallStart + 1, wallEnd + 1},
        {holeEnd + 1, wallStart + 1, wallEnd + 1}
};
const Segment wallDownY[] = {
        {holeStart, wallStart + 1, wallEnd + 1},
        {holeEnd,   wallStart + 1, wallEnd + 1},
        {K_y - 1,   0,             wallEnd + 1}
};

int index(int k_x, int k_y, int p_x, int p_y) {
    return (p_y + P_y) * (2 * P_x + 1) * K_y * K_x
           + (p_x + P_x) * K_y * K_x
           + k_y * K_x
           + k_x;
}

void fillSkipMask() {
    for (Segment s: wallUpY) {
        int k_y = s.coordinate;
        for (int k_x = s.start; k_x < s.end; k_x++) {
            for (int p_y = 0; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    skipMaskY[index(k_x, k_y, p_x, p_y)] = true;
                }
            }
        }
    }

    for (Segment s: wallDownY) {
        int k_y = s.coordinate;
        for (int k_x = s.start; k_x < s.end; k_x++) {
            for (int p_y = -P_y; p_y < 0; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    skipMaskY[index(k_x, k_y, p_x, p_y)] = true;
                }
            }
        }
    }

    for (Segment s: wallRightX) {
        int k_x = s.coordinate;
        for (int k_y = s.start; k_y < s.end; k_y++) {
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = 0; p_x <= P_x; p_x++) {
                    skipMaskX[index(k_x, k_y, p_x, p_y)] = true;
                }
            }
        }
    }

    for (Segment s: wallLeftX) {
        int k_x = s.coordinate;
        for (int k_y = s.start; k_y < s.end; k_y++) {
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x < 0; p_x++) {
                    skipMaskX[index(k_x, k_y, p_x, p_y)] = true;
                }
            }
        }
    }
}

double initial(int k_x, int k_y, int p_x, int p_y) {
    if (k_x < wallStart) return std::exp(-((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
    else return 1e-9;
}

void set_initials(double data[]) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {
        for (int p_x = -P_x; p_x <= P_x; p_x++) {
            for (int k_y = 0; k_y < K_y; k_y++) {
                for (int k_x = 0; k_x < K_x; k_x++) {
                    data[index(k_x, k_y, p_x, p_y)] = initial(k_x, k_y, p_x, p_y);
                }
            }
        }
    }
}

double calculate_denom_x() {
    double denom = 0;
    for (int p_y = -P_y; p_y <= P_y; p_y++) for (int p_x = 0; p_x <= P_x; p_x++) denom += std::exp(- ((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.) * p_x;
    return denom;
}

double calculate_denom_y() {
    double denom = 0;
    for (int p_y = 0; p_y <= P_y; p_y++) for (int p_x = -P_x; p_x <= P_x; p_x++) denom += std::exp(- ((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.) * p_y;
    return denom;
}

double denom_x = calculate_denom_x();
double denom_y = calculate_denom_y();

void make_iteration_x(double next[], double data[]) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {
        for (int p_x = -P_x; p_x <= P_x; p_x++) {
            for (int k_y = 0; k_y < K_y; k_y++) {
                for (int k_x = 0; k_x < K_x; k_x++) {
                    if (skipMaskX[index(k_x, k_y, p_x, p_y)]) continue;
                    double g = p_x * 1.0 / P_x;
                    if (p_x > 0) {
                        next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] - g *
                                                                                            (data[index(k_x, k_y, p_x,
                                                                                                        p_y)] -
                                                                                             data[index(k_x - 1, k_y,
                                                                                                        p_x, p_y)]);
                    } else {
                        next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] - g *
                                                                                            (data[index(k_x + 1, k_y,
                                                                                                        p_x, p_y)] -
                                                                                             data[index(k_x, k_y, p_x,
                                                                                                        p_y)]);
                    }
                }
            }
        }
    }

    // Zero border condition
//    for (int p_y = -P_y; p_y <= P_y; p_y++) {
//        for (int p_x = -P_x; p_x <= P_x; p_x++) {
//            for (int k_y = 0; k_y < K_y; k_y++) {
//                if (p_x < 0) {
//                    next[index(K_x - 1, k_y, p_x, p_y)] = data[index(K_x - 1, k_y, p_x, p_y)];
//                }
//            }
//        }
//    }

    for (int p_y = -P_y; p_y <= P_y; p_y++) {
        for (int p_x = -P_x; p_x <= P_x; p_x++) {
            for (int k_y = 0; k_y < K_y; k_y++) {
                if (p_x < 0) {
                    next[index(K_x - 1, k_y, p_x, p_y)] = next[index(K_x - 1, k_y, -p_x, p_y)];
                }
            }
        }
    }

    // Mirror reflection
//    for (Segment s: wallRightX) {
//        int k_x = s.coordinate;
//        for (int k_y = s.start; k_y < s.end; k_y++) {
//            for (int p_y = -P_y; p_y <= P_y; p_y++) {
//                for (int p_x = 0; p_x <= P_x; p_x++) {
//                    next[index(k_x, k_y, p_x, p_y)] = next[index(k_x, k_y, -p_x, p_y)];
//                }
//            }
//        }
//    }
//    for (Segment s: wallLeftX) {
//        int k_x = s.coordinate;
//        for (int k_y = s.start; k_y < s.end; k_y++) {
//            for (int p_y = -P_y; p_y <= P_y; p_y++) {
//                for (int p_x = -P_x; p_x <= 0; p_x++) {
//                    next[index(k_x, k_y, p_x, p_y)] = next[index(k_x, k_y, -p_x, p_y)];
//                }
//            }
//        }
//    }

    // Diffuse reflection
    for (Segment s: wallRightX) {
        int k_x = s.coordinate;
        for (int k_y = s.start; k_y < s.end; k_y++) {
            double nom = 0;
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= 0; p_x++) {
                    nom += next[index(k_x, k_y, p_x, p_y)] * p_x;
                }
            }
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = 0; p_x <= P_x; p_x++) {
                    next[index(k_x, k_y, p_x, p_y)] = -nom / denom_x * std::exp(
                            -((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
                }
            }
        }
    }

    for (Segment s: wallLeftX) {
        int k_x = s.coordinate;
        for (int k_y = s.start; k_y < s.end; k_y++) {
            double nom = 0;
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = 0; p_x <= P_x; p_x++) {
                    nom += next[index(k_x, k_y, p_x, p_y)] * p_x;
                }
            }
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= 0; p_x++) {
                    next[index(k_x, k_y, p_x, p_y)] = nom / denom_x * std::exp(- ((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
                }
            }
        }
    }

}

void make_iteration_y(double next[], double data[]) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {
        for (int p_x = -P_x; p_x <= P_x; p_x++) {
            for (int k_y = 0; k_y < K_y; k_y++) {
                for (int k_x = 0; k_x < K_x; k_x++) {
                    if (skipMaskY[index(k_x, k_y, p_x, p_y)]) {
                        continue;
                    }
                    double g = p_y * 1.0 / P_y;
                    if (p_y > 0) {
                        next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] - g *
                                                                                            (data[index(k_x, k_y, p_x,
                                                                                                        p_y)] -
                                                                                             data[index(k_x, k_y - 1,
                                                                                                        p_x, p_y)]);
                    } else {
                        next[index(k_x, k_y, p_x, p_y)] = data[index(k_x, k_y, p_x, p_y)] - g *
                                                                                            (data[index(k_x, k_y + 1,
                                                                                                        p_x, p_y)] -
                                                                                             data[index(k_x, k_y, p_x,
                                                                                                        p_y)]);
                    }
                }
            }
        }
    }

    // Zero border condition
//    for (int p_y = -P_y; p_y <= P_y; p_y++) {
//        for (int p_x = -P_x; p_x <= P_x; p_x++) {
//            for (int k_x = wallEnd + 1; k_x < K_x; k_x++) {
//                if (p_y > 0) {
//                    next[index(k_x, 0, p_x, p_y)] = data[index(k_x, 0, p_x, p_y)];
//                } else {
//                    next[index(k_x, K_y - 1, p_x, p_y)] = data[index(k_x, K_y - 1, p_x, p_y)];
//                }
//            }
//        }
//    }

    for (int p_y = -P_y; p_y <= P_y; p_y++) {
        for (int p_x = -P_x; p_x <= P_x; p_x++) {
            for (int k_x = wallEnd + 1; k_x < K_x; k_x++) {
                if (p_y > 0) {
                    next[index(k_x, 0, p_x, p_y)] = next[index(k_x, 0, p_x, -p_y)];
                } else {
                    next[index(k_x, K_y - 1, p_x, p_y)] = next[index(k_x, K_y - 1, p_x, -p_y)];
                }
            }
        }
    }

    // Mirror reflection
//    for (Segment s: wallUpY) {
//        int k_y = s.coordinate;
//        for (int k_x = s.start; k_x < s.end; k_x++) {
//            for (int p_y = 0; p_y <= P_y; p_y++) {
//                for (int p_x = -P_x; p_x <= P_x; p_x++) {
//                    next[index(k_x, k_y, p_x, p_y)] = next[index(k_x, k_y, p_x, -p_y)];
//                }
//            }
//        }
//    }
//    for (Segment s: wallDownY) {
//        int k_y = s.coordinate;
//        for (int k_x = s.start; k_x < s.end; k_x++) {
//            for (int p_y = -P_y; p_y <= 0; p_y++) {
//                for (int p_x = -P_x; p_x <= P_x; p_x++) {
//                    next[index(k_x, k_y, p_x, p_y)] = next[index(k_x, k_y, p_x, -p_y)];
//                }
//            }
//        }
//    }

    // Diffuse reflection
    for (Segment s: wallUpY) {
        int k_y = s.coordinate;
        for (int k_x = s.start; k_x < s.end; k_x++) {
            double nom = 0;

            for (int p_y = -P_y; p_y <= 0; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    nom += next[index(k_x, k_y, p_x, p_y)] * p_y;
                }
            }
            for (int p_y = 0; p_y <= P_y; p_y++) {
                for (int p_x = P_x; p_x <= P_x; p_x++) {
                    next[index(k_x, k_y, p_x, p_y)] = -nom / denom_y * std::exp(
                            -((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
                }
            }
        }
    }

    for (Segment s: wallDownY) {
        int k_y = s.coordinate;
        for (int k_x = s.start; k_x < s.end; k_x++) {
            double nom = 0;
            for (int p_y = 0; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    nom += next[index(k_x, k_y, p_x, p_y)] * p_y;
                }
            }
            for (int p_y = -P_y; p_y <= 0; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    next[index(k_x, k_y, p_x, p_y)] = nom / denom_y * std::exp(- ((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
                }
            }

        }
    }
}

void write_total_count(char *filename, double *data, int i) {
    std::ofstream f;
    if (i == 0) {
        f.open(filename);
    } else {
        f.open(filename, std::ios_base::app);
    }

    double total_count = 0;
    for (int k_x = 0; k_x < K_x; k_x++) {
        for (int k_y = 0; k_y < K_y; k_y++) {

            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    total_count += data[index(k_x, k_y, p_x, p_y)];
                }
            }
        }
    }
    f << total_count << std::endl;
    f.close();
}

void write_to_file(char *filename, double *data) {
    std::ofstream f(filename);

    for (int k_x = 0; k_x < K_x; k_x++) {
        for (int k_y = 0; k_y < K_y; k_y++) {

            double concentration = 0;
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    concentration += data[index(k_x, k_y, p_x, p_y)];
                }
            }
            f << k_x << " " << k_y << " " << dv_x * dv_y * concentration << std::endl;
        }
    }

    f.close();
}


int main() {
    auto size = K_x * K_y * (2 * P_x + 1) * (2 * P_y + 1);
    auto data = new double[size];
    auto next = new double[size];
    skipMaskX = new bool [size]();
    skipMaskY = new bool [size]();
    fillSkipMask();

    set_initials(data);
    for (int i = 0; i < N; i++) {
        make_iteration_x(next, data);
        std::swap(next, data);
        make_iteration_y(next, data);
        std::swap(next, data);

        if (i % step == 0) {
            char filename[30];
            sprintf(filename, "data/final_gp/out_%03d.dat", i);
            write_to_file(filename, data);
            char totalFilename[40] = "data/final_count/total_count.dat";
            write_total_count(totalFilename, data, i);
        }
    }
}
