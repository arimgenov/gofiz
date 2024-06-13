#include <cmath>
#include <fstream>
#include <iostream>

const int step = 10;
const int N = 1500, K_x = 100, K_y = 100; // K - размер пространства
const int P_x = 30, P_y = P_x; // размеры скорости
const int a = 10, c_x = 50, c_y = 50; // c-положение пика
const int p0_x = 3, p0_y = 2; // начальные скорости

const double dv_x = 5.0/P_x, dv_y= dv_x;

const int ws = 49, we = 50; // координаты стены
const int rs = 49, re = 50;
const int hs = 39, he = 61;
// const int ts = 40, te = 60;

const int wrx[] = {0, 20, 20}; const int wrx_b[] = {0, 0, 60}; const int wrx_e[] = {K_y-1, 40, K_y-1}; 
const int wlx[] = {K_x-1, 19, 19}; const int wlx_b[] = {0, 0, 60}; const int wlx_e[] = {K_y-1, 40, K_y-1};
const int wuy[] = {0}; const int wuy_b[] = {0}; const int wuy_e[] = {K_x-1};
const int wdy[] = {K_y-1}; const int wdy_b[] = {0}; const int wdy_e[] = {K_x-1};

int index(int k_x, int k_y, int p_x, int p_y) {
    return (p_y + P_y) * (2 * P_x + 1) * K_x * K_y
    + (p_x + P_x) * K_x * K_y
    + k_y * K_x
    + k_x;
}

void set_initials1(double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        if (p_x == p0_x and p_y == p0_y) data[index(k_x,k_y,p_x,p_y)] = std::exp(-1.*((k_x-c_x)*(k_x-c_x)+(k_y-c_y)*(k_y-c_y))/(1.*a*a));
        // if (k_x < 60 and k_x > 40 and k_y < 60 and k_y > 40) data[index(k_x,k_y,p_x,p_y)] = std::exp(-1.*((k_x-c_x)*(k_x-c_x)+(k_y-c_y)*(k_y-c_y))/(1.*a*a));
        else data[index(k_x,k_y,p_x,p_y)] = 1e-9;
    }
    } } } }
}

void set_initials2(double* data, bool f) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        // if (k_x < ws) data[index(k_x,k_y,p_x,p_y)] = std::exp(-((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
        if (k_x < 60 and k_x > 40 and k_y < 60 and k_y > 40) data[index(k_x,k_y,p_x,p_y)] = std::exp(-((p_x * dv_x) * (p_x * dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
        else if (f) data[index(k_x,k_y,p_x,p_y)] = 1e-9;
    }
    } } } }
}

double calculate_denom_x() {
    double denom = 0;
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = 0; p_x <= P_x; p_x++) {
        denom += std::exp(-((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.) * p_x;
    } }
    return denom;
}

double calculate_denom_y() {
    double denom = 0;
    for (int p_y = 0; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {
        denom += std::exp(-((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.) * p_y;
    } }
    return denom;
}

double denom_x = calculate_denom_x();
double denom_y = calculate_denom_y();

void make_iter_x(double* next, double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        double gam = 1.* p_x /P_x;
        if (p_x > 0) {
            if (k_x == 0) continue; // граничное условие
            else next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y,p_x,p_y)] - data[index(k_x - 1,k_y,p_x,p_y)]);
        } else {
            if (k_x == K_x-1) continue;
            else next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x+1,k_y,p_x,p_y)] - data[index(k_x,k_y,p_x,p_y)]);
        }
    }
    } } } }
    // зеркальное граничное уловие
    // for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {
    // {
    //     if (p_x > 0) for(int k_x : wrx) next[index(k_x,k_y,p_x,p_y)] = next[index(k_x,k_y,-p_x,p_y)];
    //     else for (int k_x : wlx) next[index(k_x,k_y,p_x,p_y)] = next[index(k_x,k_y,-p_x,p_y)];
    // }
    // } } }
    // диффузное граничное условие
    for(int k_y = 0; k_y < K_y; k_y++) {
        for (int w=0; w<(sizeof(wrx)/sizeof(wrx[0])); w++) {
            int k_x = wrx[w];
            if (k_y<wrx_b[w] or k_y>wrx_e[w]) continue; // проверка на размеры стены
            double nom = 0;
            for (int p_y = -P_y; p_y <= P_y; p_y++) for (int p_x = -P_x; p_x < 0; p_x++) nom += next[index(k_x, k_y, p_x, p_y)] * p_x;
            for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = 0; p_x <= P_x; p_x++) {
                next[index(k_x, k_y, p_x, p_y)] = -nom / denom_x * std::exp(-((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
            } }
        }
        for (int w=0; w<(sizeof(wlx)/sizeof(wlx[0])); w++) {
            int k_x = wlx[w];
            if (k_y<wlx_b[w] or k_y>wlx_e[w]) continue; // проверка на размеры стены
            double nom = 0;
            for (int p_y = -P_y; p_y <= P_y; p_y++) for (int p_x = 0; p_x <= P_x; p_x++) nom += next[index(k_x, k_y, p_x, p_y)] * p_x;
            for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x < 0; p_x++) {
                next[index(k_x, k_y, p_x, p_y)] = nom / denom_x * std::exp(-((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
            } }
        }
    }
}

void make_iter_y(double* next, double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        double gam = 1.* p_y /P_y;
        if (p_y > 0) {
            if (k_y == 0) continue; // граничное условие
            else next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y,p_x,p_y)] - data[index(k_x,k_y-1,p_x,p_y)]);
        } else {
            if (k_y == K_y-1) continue;
            else next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y+1,p_x,p_y)] - data[index(k_x,k_y,p_x,p_y)]);
        }   
    }
    } } } }
    // зеркальное граничное уловие
    // for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_x = 0; k_x < K_x; k_x++) {
    // {
    //     if (p_y > 0) for(int k_y : wdy) next[index(k_x,k_y,p_x,p_y)] = next[index(k_x,k_y,p_x,-p_y)];
    //     else for (int k_y : wuy) next[index(k_x,k_y,p_x,p_y)] = next[index(k_x,k_y,p_x,-p_y)];
    // }
    // } } }
    // диффузное граничное условие
    for(int k_x = 0; k_x < K_x; k_x++) {
        for (int w=0; w<(sizeof(wuy)/sizeof(wuy[0])); w++) {
            int k_y = wuy[w];
            if (k_x<wuy_b[w] or k_x>wuy_e[w]) continue; // проверка на размеры стены
            double nom = 0;
            for (int p_y = -P_y; p_y < 0; p_y++) for (int p_x = -P_x; p_x <= P_x; p_x++) nom += next[index(k_x, k_y, p_x, p_y)] * p_y;
            for (int p_y = 0; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {
                next[index(k_x, k_y, p_x, p_y)] = -nom / denom_x * std::exp(-((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
            } }
        }
        for (int w=0; w<(sizeof(wdy)/sizeof(wdy[0])); w++) {
            int k_y = wdy[w];
            if (k_x<wdy_b[w] or k_x>wdy_e[w]) continue; // проверка на размеры стены
            double nom = 0;
            for (int p_y = 0; p_y <= P_y; p_y++) for (int p_x = -P_x; p_x <= P_x; p_x++) nom += next[index(k_x, k_y, p_x, p_y)] * p_y;
            for (int p_y = -P_y; p_y < 0; p_y++) {for (int p_x = -P_x; p_x < P_x; p_x++) {
                next[index(k_x, k_y, p_x, p_y)] = nom / denom_x * std::exp(-((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) / 2.);
            } }
        }
    }
}


void write_total_count(char *filename, double total, int i) {
    std::ofstream f;
    if (i == 0) f.open(filename);
    else f.open(filename, std::ios_base::app);

    f << total << std::endl;
    f.close();
}

void write_to_file(double* data, char* filename) {
    std::ofstream out(filename);
    
    for (int k_y = 0; k_y < K_y; k_y++) {
        for (int k_x = 0; k_x < K_x; k_x++) {
            double con = 0;
            
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    con += data[index(k_x,k_y,p_x,p_y)];
            } }

            out << k_x << " " << k_y << " " << dv_x * dv_y * con << std::endl;
    } }

    out.close();
}

int main() {
    auto size = K_x*K_y*(2*P_x+1)*(2*P_y+1);
    auto data = new double[size];
    auto next = new double[size];

    // set_initials2(data, true);
    set_initials1(data);

    double total = 0;
    // double total_next = 0;
    // for (int k_x = 0; k_x <= we; k_x++) for (int k_y = 0; k_y < K_y; k_y++) for (int p_y = -P_y; p_y <= P_y; p_y++) for (int p_x = -P_x; p_x <= P_x; p_x++) total += ((p_x*dv_x) * (p_x*dv_x) + (p_y * dv_y) * (p_y * dv_y)) * data[index(k_x, k_y, p_x, p_y)];

    for(int i=0; i<N; i++) {
        make_iter_x(next, data);
        std::swap(next, data);
        make_iter_y(next, data);
        std::swap(next, data);
        // set_initials2(data, false);

        if (i % step == 0) {
            char filename[50];
            sprintf(filename, "data/task3/out_%03d.dat", i); // запись в файл
            write_to_file(data, filename);

            total = 0;
            for (int k_x = 0; k_x <= we; k_x++) for (int k_y = 0; k_y < K_y; k_y++) for (int p_y = -P_y; p_y <= P_y; p_y++) for (int p_x = -P_x; p_x <= P_x; p_x++) total += (p_x*p_x * dv_x*dv_x + dv_y*dv_y*p_y*p_y) * data[index(k_x, k_y, p_x, p_y)];

            char totalFilename[50] = "data/final_count/total_count.dat";
            write_total_count(totalFilename, total, i);
        }
    }
    
    return 0;
}