#include <cmath>
#include <fstream>
#include <iostream>

const int step = 10;
const int N = 1500, K_x = 100, K_y = 100; // K - размер пространства
const int P_x = 13, P_y = 13; // размеры скорости
const int a = 20, c_x = 50, c_y = 50; // c-положение пика
const int p0_x = 2, p0_y = 0; // начальные скорости

const double dv_x = 5.0/P_x, dv_y= 5.0/P_y;

const int wrx[] = {0, 20, 20}; const int wrx_b[] = {0, 52, 0}; const int wrx_e[] = {K_y-1, K_y-1, 48}; 
const int wlx[] = {K_x-1, 19, 19}; const int wlx_b[] = {0, 52, 0}; const int wlx_e[] = {K_y-1, K_y-1, 48};
const int wuy[] = {0, 48}; const int wuy_b[] = {0, 19}; const int wuy_e[] = {K_x-1, 20};
const int wdy[] = {K_y-1, 52}; const int wdy_b[] = {0, 19}; const int wdy_e[] = {K_x-1, 20};

int index(int k_x, int k_y, int p_x, int p_y) {             
    return (p_y + P_y) * (2 * P_x + 1) * K_x * K_y
    + (p_x + P_x) * K_x * K_y
    + k_y * K_x
    + k_x;
}

void set_initials(double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        if (p_x == p0_x and p_y == p0_y) data[index(k_x,k_y,p_x,p_y)] = std::exp(-1. * ((k_x-c_x)*(k_x-c_x) + (k_y-c_y)*(k_y-c_y))/(1.*a*a));
        else data[index(k_x,k_y,p_x,p_y)] = 1e-9;
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
    
    set_initials(data);

    for(int i=0; i<N; i++) {
        make_iter_x(next, data);
        std::swap(next, data);
        make_iter_y(next, data);
        std::swap(next, data);

        if (i % step == 0) {
            char filename[30];
            sprintf(filename, "D:/VisCode/gofiz/data/out_%03d.dat", i); // запись в файл
            write_to_file(data, filename);
        }
    }
    
    return 0;
}